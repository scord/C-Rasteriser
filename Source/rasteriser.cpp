#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <glm/gtx/rotate_vector.hpp>
#include "glm/ext.hpp"
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::mat3;
using glm::vec2;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 700;
const int SCREEN_HEIGHT = 700;
SDL_Surface* screen;
int t;
float focalLength = SCREEN_HEIGHT/1.2;
vector<Triangle> triangles;
vec3 cameraPos(0,0,-4);
vec3 currentColor;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
float stencilBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 frameBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 lightPos(0,-0.5 ,-0.7);
vec3 lightPower = 14.0f*vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );

mat3 R = glm::mat3(0.1f);

const float yaw = 0.07;

vector<Triangle> shadowVolume;


struct Vertex
{
  vec2 texCo;
    vec3 position;
    vec3 normal;
    vec3 color;
    vec3 xtangent;
};

struct Pixel
{
    int x;
    int y;
    float zinv;
    vec3 color;
    vec3 pos3d;
    vec3 normal;
    vec2 texCo;
    vec3 xtangent;
};
vector<Pixel> framePixels;
/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();

void PixelShader( const Pixel& p, int lighting );
void VertexShader( const Vertex& v, Pixel& p );
void DrawPolygonEdges( const vector<vec3>& vertices );
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels );
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result, vec3** tex, vec3** bump);
void ComputePolygon( const vector<vec3>& vertices, vec3** tex, vec3** bump );
void ComputePolygonPixels( const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels , vec3** tex, vec3** bump);
void DrawStencil( const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels );

void savePos(){
  ofstream myfile;
  myfile.open ("pos.data");


  myfile << cameraPos.x;
  myfile  <<  " ";

  myfile << cameraPos.y;
  myfile  <<  " ";

  myfile << cameraPos.z;
  myfile  <<  " " << endl;

  myfile << lightPos.x;
  myfile  <<  " ";

  myfile << lightPos.y;
  myfile  <<  " ";

  myfile << lightPos.z;
  myfile  <<  " " << endl;

  for (int x = 0; x < 3;x++){
        for (int y = 0; y < 3;y++){
            myfile << R[x][y];
            myfile <<  " ";
        }
  }

  myfile << endl;

  myfile.close();
}

void readPos(){
  ifstream myfile;
  myfile.open ("pos.data");

  if (myfile.is_open()){

      string buf; // Have a buffer string

      string line;

      getline(myfile, line);

      stringstream ss(line); // Insert the string into a stream


      ss >> cameraPos.x;


      ss >> cameraPos.y;



      ss >> cameraPos.z;

      getline(myfile, line);

      ss.clear();
	  ss.str(line); // Insert the string into a stream


      ss >> lightPos.x;

      ss >> lightPos.y;

      ss >> lightPos.z;

      getline(myfile, line);


      ss.clear();
	  ss.str(line); // Insert the string into a stream

      for (int x = 0; x < 3;x++){
            for (int y = 0; y < 3;y++){
                      ss >> R[x][y];

            }
      }

      myfile.close();

  }
}


void ComputePolygon( const vector<Vertex>& vertices, vec3** tex, vec3** bump)
{
    int V = vertices.size();

    vector<Pixel> vertexPixels( V );

    for( int i=0; i<V; ++i )
	{
        VertexShader( vertices[i], vertexPixels[i] );
	}

    vector<Pixel> leftPixels;
    vector<Pixel> rightPixels;
    ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
    ComputePolygonPixels( leftPixels, rightPixels, tex, bump );
}

void DrawStencilBuffer( const vector<Vertex>& vertices )
{
    int V = vertices.size();

    vector<Pixel> vertexPixels( V );

    for( int i=0; i<V; ++i )
	{
        VertexShader( vertices[i], vertexPixels[i] );
	}

    vector<Pixel> leftPixels;
    vector<Pixel> rightPixels;
    ComputePolygonRows( vertexPixels, leftPixels, rightPixels);
    DrawStencil( leftPixels, rightPixels );
}

void ComputeFacing( Triangle& face, float& indicator )
{
    indicator = glm::dot(face.v0 - lightPos, face.normal);
}

void ComputePolygonPixels( const vector<Pixel>& leftPixels,
               const vector<Pixel>& rightPixels, vec3** tex, vec3** bump)
{

	for (int i = 0; i < leftPixels.size(); i++)
	{
		int cols = rightPixels[i].x - leftPixels[i].x + 1;
		vector<Pixel> line( rightPixels[i].x - leftPixels[i].x + 1 );

      Interpolate(leftPixels[i], rightPixels[i], line, tex, bump);

		for (int j = 0; j < cols; j++ )
		{
			Pixel p = line[j];
            framePixels.push_back(p);
			//PixelShader(p);
		}

	}
}


void DrawStencil( const vector<Pixel>& leftPixels,
               const vector<Pixel>& rightPixels )
{
	for (int i = 0; i < leftPixels.size(); i++)
	{
		int cols = rightPixels[i].x - leftPixels[i].x + 1;
		vector<Pixel> line( rightPixels[i].x - leftPixels[i].x + 1 );
		Interpolate(leftPixels[i], rightPixels[i], line, NULL, NULL);

		for (int j = 0; j < cols; j++ )
		{
			Pixel p = line[j];

            if (p.y < SCREEN_HEIGHT && p.y > 0 && p.x > 0 && p.x < SCREEN_WIDTH)
        	{

                if (depthBuffer[p.y][p.x] < p.zinv)
                {

                    stencilBuffer[p.y][p.x] += glm::dot(p.pos3d, p.normal);

                }

        	}
		}
	}
}

void PixelShader( const Pixel& p, int lighting )
{
	if (p.y < SCREEN_HEIGHT && p.y > 0 && p.x > 0 && p.x < SCREEN_WIDTH && depthBuffer[p.y][p.x] < p.zinv)
	{
        vec3 r = glm::normalize(lightPos - cameraPos -  p.pos3d);
    	float dist = glm::distance(lightPos - cameraPos, p.pos3d);

        vec3 illumination = p.color * indirectLightPowerPerArea;

        if (lighting == 1)
        {
    	    vec3 diff = lightPower / (4.0f*3.1416f*dist*dist);
            vec3 d = lightPos - cameraPos - p.pos3d;
            vec3 r = 2*glm::dot(d, p.normal)*p.normal - d;
            vec3 v = p.pos3d;
            vec3 h = (d + v) / glm::sqrt(glm::dot(d + v, d + v));
            float phong = glm::max(0.0f, glm::dot(r, - p.pos3d));
            float blinn = pow(glm::dot(p.normal, h),5);
    	    illumination += p.color * (0.1f*diff+0.9f*blinn);
        }
        frameBuffer[p.y][p.x] = illumination;
		depthBuffer[p.y][p.x] = p.zinv;
	}
}
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels )
{
// 1. Find max and min y-value of the polygon
//    and compute the number of rows it occupies.

	int maxY = -numeric_limits<int>::max();
	int minY = +numeric_limits<int>::max();

	for (int i = 0; i < vertexPixels.size(); ++i)
	{
		if (vertexPixels[i].y > maxY)
			maxY = vertexPixels[i].y;
		if (vertexPixels[i].y < minY)
			minY = vertexPixels[i].y;
	}

	int numRows = maxY - minY + 1;


// 2. Resize leftPixels and rightPixels
//    so that they have an element for each row.
	leftPixels.resize(numRows);
	rightPixels.resize(numRows);
// 3. Initialize the x-coordinates in leftPixels
//    to some really large value and the x-coordinates
//    in rightPixels to some really small value.

	for( int i=0; i<numRows; ++i )
	{
	    leftPixels[i].x  = +numeric_limits<int>::max();
	    rightPixels[i].x = -numeric_limits<int>::max();
	}



// 4. Loop through all edges of the polygon and use
//    linear interpolation to find the x-coordinate for
//    each row it occupies. Update the corresponding
//    values in rightPixels and leftPixels.

	ivec2 avec;
	ivec2 bvec;

	// Loop over all vertices and draw the edge from it to the next vertex:
	for( int i=0; i<vertexPixels.size(); ++i )
	{

		int i2 = (i+1)%vertexPixels.size(); // The next vertex

		Pixel a = vertexPixels[i];
		Pixel b = vertexPixels[i2];

		avec = ivec2(a.x, a.y);
		bvec = ivec2(b.x, b.y);

		int deltaX = glm::abs(a.x - b.x);
		int deltaY = glm::abs(a.y - b.y);

		int pixels = glm::max( deltaX, deltaY ) + 1;

		vector<Pixel> line( pixels );

		Interpolate( a, b, line , NULL, NULL);

		for (int j=0; j < pixels; ++j)
		{
			Pixel pixel = line[j];
			int y = pixel.y - minY;

            if ( y < leftPixels.size())
            {

    			if (leftPixels[y].x > pixel.x)
    				leftPixels[y] = pixel;
    			if (rightPixels[y].x < pixel.x)
    				rightPixels[y] = pixel;
            }

		}


	}

}

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	LoadTestModel(triangles);


	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

    Uint8* keystate = SDL_GetKeyState( 0 );
	if (keystate[SDLK_w])
	{
		lightPos.z += 0.1;
	}
	if (keystate[SDLK_s])
	{
		lightPos.z -= 0.1;
	}
	if (keystate[SDLK_a])
	{
		lightPos.x -= 0.1;
	}
	if (keystate[SDLK_d])
	{
		lightPos.x += 0.1;
	}
	if (keystate[SDLK_e])
	{
		lightPos.y += 0.1;
	}
	if (keystate[SDLK_q])
	{
		lightPos.y -= 0.1;
	}

  if( keystate[SDLK_r] )
{
cameraPos = vec3(0,0,-2);
lightPos = vec3(0, -0.5, -0.7 );
R = glm::mat3(0.1f);
// dcrease focal
}

if( keystate[SDLK_UP] )
{
cameraPos+= vec3( R[2][0], R[2][1], R[2][2] )*3.0f;
// Move camera forward
}
if( keystate[SDLK_DOWN] )
{
cameraPos-= vec3( R[2][0], R[2][1], R[2][2] )*3.0f;
// Move camera backward
}
if( keystate[SDLK_LEFT] )
{
R[2] = glm::rotateY(glm::vec3( R[2][0], R[2][1], R[2][2] ), -yaw);
R[0] = glm::rotateY(glm::vec3( R[0][0], R[0][1], R[0][2] ), -yaw);
// rotate camera to the left
}
if( keystate[SDLK_RIGHT] )
{
R[2] = glm::rotateY(glm::vec3( R[2][0], R[2][1], R[2][2] ), yaw);
R[0] = glm::rotateY(glm::vec3( R[0][0], R[0][1], R[0][2] ), yaw);
// totate camera to the right
}

}

void VertexShader( const Vertex& v, Pixel& p)
{
	vec3 v2 = v.position - cameraPos;
    vec3 forward = vec3(R[2][0], R[2][1], R[2][2]);
    v2 = glm::rotateY(v2, -atan2(forward.x, forward.z));
	p.x = focalLength * (v2.x / v2.z) + SCREEN_WIDTH/2.0f;
	p.y = focalLength * (v2.y / v2.z) + SCREEN_HEIGHT/2.0f;
	p.zinv = 1.0f/v2.z;
	p.pos3d = v2;
    forward = vec3(R[2][0], R[2][1], R[2][2]);
    p.normal = v.normal;
    p.xtangent = glm::rotateY(v.xtangent, atan2(forward.x, forward.z));
    p.color = v.color;
    p.texCo = v.texCo;
}

void Interpolate( Pixel a, Pixel b, vector<Pixel>& result , vec3** tex, vec3** bump)
{
    cout << a.color.r << endl;

    int N = result.size();
    float stepX = (float)(b.x-a.x) / float(max(N-1,1));
	float stepY = (float)(b.y-a.y) / float(max(N-1,1));
	float stepZ = (b.zinv-a.zinv) / float(max(N-1,1));
    vec3 stepN = (b.normal-a.normal) * 1.0f/float(max(N-1,1));
	vec3 stepI = vec3(b.pos3d*b.zinv-a.pos3d*a.zinv) / float(max(N-1,1));
    vec3 stepC = (b.color-a.color) *1.0f/float(max(N-1,1));
    vec2 stepT = (b.texCo*b.zinv - a.texCo*a.zinv)*1.0f/float(max(N-1,1));

    //cout << a.texCo.x << endl;

	Pixel current;
	float currentx = a.x;
	float currenty = a.y;
	float currentzinv = a.zinv;
    vec3 currentN = a.normal;

    current.color = a.color;
  current.pos3d = a.pos3d*a.zinv;
  current.texCo = a.texCo*a.zinv;//*a.zinv;
//cerr << "b " << b.texCo.x << " " << b.texCo.y << endl;
    for( int i=0; i<N; ++i )
    {
      //cerr << N << " " << (int)(glm::round(current.texCo.x)) << " " << (int)(glm::round(current.texCo.y)) <<  " " << stepT.x << " " << stepT.y << " " << i << endl;
        current.x = currentx;
		current.y = currenty;
		current.zinv = currentzinv;
current.xtangent = a.xtangent;

        if (tex != NULL){
        //  if ((int)(glm::round(current.texCo.x))  < 0 || (int)(glm::round(current.texCo.y))  < 0){
         // cerr << (int)(glm::round(current.texCo.x)) << " " << (int)(glm::round(current.texCo.y)) << " " << stepT.x << " " << stepT.y << " " << i << endl;
          //cerr << "b " << b.texCo.x << " " << b.texCo.y << endl;
        //}

          current.color = tex[(int)(current.texCo.x/current.zinv - 0.9)][(int)(current.texCo.y/current.zinv - 0.9)];

        } else {
          current.color += stepC;
        }

        if (bump != NULL){
            vec3 inv = bump[(int)(current.texCo.x/current.zinv - 0.9)][(int)(current.texCo.y/current.zinv - 0.9)];
            inv = vec3(inv.x, inv.y, inv.z);
            //current.normal = inv*2.0f - 1.0f - vec3(0,0,1) + currentN;

          current.normal = (bump[(int)(current.texCo.x/current.zinv - 0.9)][(int)(current.texCo.y/current.zinv - 0.9)]*2.0f - 1.0f);//glm::rotateY(bump[(int)(current.texCo.x/current.zinv - 0.9)][(int)(current.texCo.y/current.zinv - 0.9)]*2.0f - 1.0f, atan2(currentN.x, currentN.z));
          //vec3 tmp = glm::rotateY(currentN, atan2(currentN.x, currentN.z));
          //current.normal = glm::rotateX(current.normal, atan2(tmp.y, tmp.z));
          //current.normal = glm::rotateX(glm::rotateY(current.normal, atan2(current.xtangent.x, current.xtangent.z)), atan2(current.xtangent.y, current.xtangent.z));

        } else {
          current.normal = currentN;
      }

        result[i] = current;
        //result[i].texCo = result[i].texCo/result[i].zinv;
        result[i].pos3d = result[i].pos3d/result[i].zinv;
        result[i].texCo = result[i].texCo/result[i].zinv;
        result[i].color = result[i].color * a.color;
        currentx += stepX;
		currenty += stepY;
		currentzinv += stepZ;
		current.pos3d += stepI;
        current.texCo += stepT;

        currentN += stepN;
    }

}

void DrawLineSDL( SDL_Surface* surface, Pixel a, Pixel b, vec3 color )
{
    int cols = max(abs(a.y - b.y + 1), abs(a.x - b.x + 1));
    vector<Pixel> line(cols);
    Interpolate(a, b, line, NULL, NULL);
    for (int i = 0; i<cols; ++i)
    {
        PutPixelSDL(screen, line[i].x, line[i].y, color);
    }
}
void DrawEdges( const vector<Vertex>& vertices )
{
    int V = vertices.size();
    // Transform each vertex from 3D world position to 2D image position:
    vector<Pixel> projectedVertices( V );
    for( int i=0; i<V; ++i )
    {
        VertexShader( vertices[i], projectedVertices[i] );
    }
    // Loop over all vertices and draw the edge from it to the next vertex:
    for( int i=0; i<V; i+=1 )
    {
        int j = (i+1)%V; // The next vertex
        vec3 color( 1, 1, 0 );
        DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], color );
	}
}

int CheckEdge(Triangle triangle, Vertex v1, Vertex v2)
{
    for( int i=0; i<triangles.size(); ++i )
    {
        vector<Vertex> vertices(3);
        vertices[0].position = triangles[i].v0;
        vertices[1].position = triangles[i].v1;
        vertices[2].position = triangles[i].v2;

        for( int j=0; j<vertices.size(); ++j )
    	{
    		int j2 = (j+1)%vertices.size();

            if (vertices[j2].position == v1.position && vertices[j].position == v2.position)
            {
                float indicator1,indicator2 = 0;
                ComputeFacing(triangle, indicator1);
                ComputeFacing(triangles[i], indicator2);


                if ((indicator1 >= 0 && indicator2 < 0))
                {
                    return 1;
                }
            }

        }
    }
    return 0;
}

int ComputeSilhouetteEdges()
{
    vector<Vertex> edges;
    Vertex lightvert;
    lightvert.position = lightPos;
    int numEdges = 0;

    for( int i=0; i<triangles.size(); ++i )
    {
        vector<Vertex> vertices(3);
        vertices[0].position = triangles[i].v0;
        vertices[1].position = triangles[i].v1;
        vertices[2].position = triangles[i].v2;

        for( int j=0; j<vertices.size(); ++j )
    	{
    		int j2 = (j+1)%vertices.size();

            int sil = CheckEdge(triangles[i], vertices[j], vertices[j2]);
            if (sil)
            {
                vector<Vertex> v(3);

                v[0].position = lightPos + (vertices[j].position - lightPos)*1.0f;
                v[1].position = lightPos + (vertices[j2].position - lightPos)*1.0f;
                v[2].position = lightPos + (vertices[j2].position - lightPos)*5.0f;
                DrawEdges(v);
                edges.push_back(v[0]);
                edges.push_back(v[1]);
                edges.push_back(v[2]);
                Triangle triangle(v[0].position, v[1].position, v[2].position, vec3(1,1,1));
                triangle.ComputeNormal();
                shadowVolume.push_back(triangle);

                v[0].position = lightPos + (vertices[j].position - lightPos)*1.0f;
                v[1].position = lightPos + (vertices[j2].position - lightPos)*5.0f;
                v[2].position = lightPos + (vertices[j].position - lightPos)*5.0f;

                Triangle triangle2(v[0].position, v[1].position, v[2].position, vec3(1,1,1));
                triangle2.ComputeNormal();
                shadowVolume.push_back(triangle2);

                DrawEdges(v);
            }
            numEdges += sil;

        }


    }
    //DrawEdges(edges);
    return numEdges;
}



void Draw() {
    SDL_FillRect( screen, 0, 0 );
    if( SDL_MUSTLOCK(screen) )
        SDL_LockSurface(screen);

    shadowVolume.clear();
    framePixels.clear();

	for( int y=0; y<SCREEN_HEIGHT; ++y )
    {
    	for( int x=0; x<SCREEN_WIDTH; ++x )
        {
        	depthBuffer[y][x] = 0;
            stencilBuffer[y][x] = 0;
            frameBuffer[y][x] = vec3(0);
        }
    }
    int numEdges = ComputeSilhouetteEdges();
    for( int i=0; i<triangles.size(); ++i )
    {
        currentColor = triangles[i].color;
        vector<Vertex> vertices(3);
        vertices[0].position = triangles[i].v0;
        vertices[0].normal = triangles[i].normal;
        vertices[0].color = triangles[i].color;
        vertices[0].texCo = triangles[i].v0texCo;
        vertices[0].xtangent = triangles[i].xtangent;
        vertices[1].position = triangles[i].v1;
        vertices[1].normal = triangles[i].normal;
        vertices[1].color = triangles[i].color;
        vertices[1].texCo = triangles[i].v1texCo;
        vertices[1].xtangent = triangles[i].xtangent;
        vertices[2].position = triangles[i].v2;
        vertices[2].normal = triangles[i].normal;
        vertices[2].color = triangles[i].color;
        vertices[2].texCo = triangles[i].v2texCo;
        vertices[2].xtangent = triangles[i].xtangent;

        ComputePolygon( vertices , triangles[i].texture, triangles[i].bump);
    }

    for (int i=0; i<framePixels.size(); ++i)
    {
        PixelShader(framePixels[i], 0);
    }

    for ( int i=0; i<shadowVolume.size(); ++i)
    {

        vector<Vertex> vertices(3);
        vertices[0].position = shadowVolume[i].v0;
        vertices[0].normal = shadowVolume[i].normal;
        vertices[0].color = vec3(0);
        vertices[1].position = shadowVolume[i].v1;
        vertices[1].normal = shadowVolume[i].normal;
        vertices[1].color = vec3(0);
        vertices[2].position = shadowVolume[i].v2;
        vertices[2].normal = shadowVolume[i].normal;
        vertices[2].color = vec3(0);

        DrawStencilBuffer( vertices );
    }

    for( int y=0; y<SCREEN_HEIGHT; ++y )
    {
    	for( int x=0; x<SCREEN_WIDTH; ++x )
        {
        	depthBuffer[y][x] = 0;
        }
    }

    for (int i=0; i<framePixels.size(); ++i)
    {
        int x = framePixels[i].x;
        int y = framePixels[i].y;

        if (stencilBuffer[y][x] == 0)
            PixelShader(framePixels[i], 1);
    }

    for (int y=0; y<SCREEN_HEIGHT; ++y)
    {
        for (int x=0; x<SCREEN_WIDTH; ++x)
        {

            PutPixelSDL( screen, x, y, frameBuffer[y][x] );


        }
    }
    //numEdges = ComputeSilhouetteEdges();
    if ( SDL_MUSTLOCK(screen) )
        SDL_UnlockSurface(screen);
    SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

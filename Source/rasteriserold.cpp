#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::mat3;
using glm::vec2;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
float focalLength = SCREEN_HEIGHT/1.2;
vector<Triangle> triangles;
vec3 cameraPos(0,0,-2);
vec3 currentColor;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 lightPos(0,0,-1);
vec3 lightPower = 1.1f*vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );


struct Vertex
{
    vec3 position;
    vec3 normal;
    float reflectance;
};

struct Pixel
{
    int x;
    int y;
    float zinv;
    vec3 illumination;
};
/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();

void PixelShader( const Pixel& p );
void VertexShader( const Vertex& v, Pixel& p );
void DrawPolygonEdges( const vector<vec3>& vertices );
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels );
void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color );
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result );
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );
void DrawPolygon( const vector<vec3>& vertices );
void DrawPolygonRows( const vector<Pixel>& leftPixels,
               const vector<Pixel>& rightPixels );


void DrawPolygon( const vector<Vertex>& vertices )
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
    DrawPolygonRows( leftPixels, rightPixels );
}


void DrawPolygonRows( const vector<Pixel>& leftPixels,
               const vector<Pixel>& rightPixels )
{
	for (int i = 0; i < leftPixels.size(); i++)
	{
		int cols = rightPixels[i].x - leftPixels[i].x + 1;
		vector<Pixel> line( rightPixels[i].x - leftPixels[i].x + 1 );
		Interpolate(leftPixels[i], rightPixels[i], line);

		for (int j = 0; j < cols; j++ )
		{
			Pixel p = line[j];

			PixelShader(p);
		}

	}
}

void PixelShader( const Pixel& p )
{
	if (p.y < SCREEN_HEIGHT && p.y > 0 && p.x > 0 && p.x < SCREEN_WIDTH && depthBuffer[p.y][p.x] < p.zinv)
	{
		PutPixelSDL( screen, p.x, p.y, p.illumination );
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

		ivec2 delta = glm::abs( avec - bvec);
		int deltaX = glm::abs(a.x - b.x);
		int deltaY = glm::abs(a.y - b.y);

		int pixels = glm::max( deltaX, deltaY ) + 1;

		vector<Pixel> line( pixels );
		Interpolate( a, b, line );

		for (int j=0; j < pixels; ++j)
		{
			Pixel pixel = line[j];
			int y = pixel.y - minY;


			if (leftPixels[y].x > pixel.x)
				leftPixels[y] = pixel;
			if (rightPixels[y].x < pixel.x)
				rightPixels[y] = pixel;

		}


	}

}

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.

	LoadTestModel(triangles);


	//while( NoQuitMessageSDL() )
	//{
		Update();
		Draw();
	//}

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
}

void VertexShader( const Vertex& v, Pixel& p)
{
	vec3 v2 = v.position - cameraPos;
	p.x = focalLength * (v2.x / v2.z) + SCREEN_WIDTH/2.0f;
	p.y = focalLength * (v2.y / v2.z) + SCREEN_HEIGHT/2.0f;
	p.zinv = 1.0f/v2.z;
	vec3 r = glm::normalize(lightPos - v.position);
	float dist = glm::distance(lightPos, v.position);
	vec3 D = lightPower*max((float)glm::dot(r, v.normal), .0f) / (4.0f*3.1416f*dist*dist);
	p.illumination = currentColor * v.reflectance * (D + indirectLightPowerPerArea);

}

void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result )
{
    int N = result.size();
    vec2 step = vec2(b-a) / float(max(N-1,1));
    vec2 current( a );
    for( int i=0; i<N; ++i )
    {
        result[i] = current;
        current += step;
    }
}

void Interpolate( Pixel a, Pixel b, vector<Pixel>& result )
{
    int N = result.size();
    float stepX = (float)(b.x-a.x) / float(max(N-1,1));
	float stepY = (float)(b.y-a.y) / float(max(N-1,1));
	float stepZ = (b.zinv-a.zinv) / float(max(N-1,1));
	vec3 stepI = vec3(b.illumination-a.illumination) / float(max(N-1,1));
	cout << stepI.x << endl;
	Pixel current;
	float currentx = a.x;
	float currenty = a.y;
	float currentzinv = a.zinv;
	current.illumination = a.illumination;

    for( int i=0; i<N; ++i )
    {
		current.x = currentx;
		current.y = currenty;
		current.zinv = currentzinv;
        result[i] = current;
        currentx += stepX;
		currenty += stepY;
		currentzinv += stepZ;
		current.illumination += stepI;
    }
}

void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color )
{
	ivec2 delta = glm::abs( a - b );
	int pixels = glm::max( delta.x, delta.y ) + 1;


	vector<ivec2> line( pixels );
	Interpolate( a, b, line );

	for ( int i=0; i < pixels; i++)
	{
		vec3 color(1,1,1);
		PutPixelSDL( screen, line[i].x, line[i].y, color );
	}

}

/*void DrawPolygonEdges( const vector<vec3>& vertices )
{
    int V = vertices.size();
    // Transform each vertex from 3D world position to 2D image position:
    vector<Pixel> projectedVertices( V );
    for( int i=0; i<V; ++i )
    {
        VertexShader( vertices[i], projectedVertices[i] );
    }
    // Loop over all vertices and draw the edge from it to the next vertex:
    for( int i=0; i<V; ++i )
    {
        int j = (i+1)%V; // The next vertex
        vec3 color( 1, 1, 1 );
        DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], color );
	}
}
*/
void Draw() {
    SDL_FillRect( screen, 0, 0 );
    if( SDL_MUSTLOCK(screen) )
        SDL_LockSurface(screen);

	for( int y=0; y<SCREEN_HEIGHT; ++y )
    	for( int x=0; x<SCREEN_WIDTH; ++x )
        	depthBuffer[y][x] = 0;

    for( int i=0; i<triangles.size(); ++i )
    {
        currentColor = triangles[i].color;
        vector<Vertex> vertices(3);
        vertices[0].position = triangles[i].v0;
		vertices[0].normal = triangles[i].normal;
		vertices[0].reflectance = 1.0f;
        vertices[1].position = triangles[i].v1;
		vertices[1].normal = triangles[i].normal;
		vertices[1].reflectance = 1.0f;
        vertices[2].position = triangles[i].v2;
		vertices[2].normal = triangles[i].normal;
		vertices[2].reflectance = 1.0f;
        DrawPolygon( vertices );
    }
    if ( SDL_MUSTLOCK(screen) )
        SDL_UnlockSurface(screen);
    SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>
#include <SDL.h>

// Used to describe a triangular surface:
class Triangle
{
public:
	glm::vec3 v0;
	glm::vec3 v1;
	glm::vec3 v2;
	glm::vec3 normal;
	glm::vec3 color;
	glm::vec3** texture;
	glm::vec3** bump;
	glm::vec2 v0texCo;
	glm::vec2 v1texCo;
	glm::vec2 v2texCo;
	glm::vec3 xtangent;


	Triangle( glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec2 v0texCo, glm::vec2 v1texCo, glm::vec2 v2texCo, glm::vec3 xtangent, glm::vec3 color, glm::vec3** texture, glm::vec3** bump)
		: v0(v0), v1(v1), v2(v2),v0texCo(v0texCo), v1texCo(v1texCo), v2texCo(v2texCo), xtangent(xtangent), color(color), texture(texture), bump(bump)
	{
		ComputeNormal();
}

Triangle( glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 color)
	: v0(v0), v1(v1), v2(v2),color(color)
{
	ComputeNormal();
}


	void ComputeNormal()
	{
		glm::vec3 e1 = v1-v0;
		glm::vec3 e2 = v2-v0;
		normal = glm::normalize( glm::cross( e2, e1) );
	}
};

glm::vec3** getTexture(const char* filename)
{
	using glm::vec3;

	SDL_Surface* surface;
	surface=SDL_LoadBMP(filename);
	Uint32 *pixels = (Uint32 *)surface->pixels;
	vec3** texture;
	texture = new vec3 * [surface->w];



	for (int i = 0; i < surface->w; i++)
	{
		texture[i] = new vec3 [surface->h];
		assert(texture[i]);
	}

	for (int y = 0; y < surface->h; y++)
	{
		for (int x = 0; x < surface->w; x++)
		{
			SDL_Color col;

			Uint32 pixel = pixels[(y*surface->w)+x];

			SDL_GetRGB(pixel, surface->format, &col.r, &col.g, &col.b);
			Uint8 *colors = (Uint8*)&pixel;
			texture[x][y] = vec3((float)col.r/255.0f, (float)col.g/255.0f, (float)col.b/255.0f);
		}
	}
	return texture;
}

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& triangles )
{
	using glm::vec3;

	SDL_Init(SDL_INIT_VIDEO);

	// Defines colors:
	vec3 red(    0.75f, 0.15f, 0.15f );
	vec3 yellow( 0.75f, 0.75f, 0.15f );
	vec3 green(  0.15f, 0.75f, 0.15f );
	vec3 cyan(   0.15f, 0.75f, 0.75f );
	vec3 blue(   0.15f, 0.15f, 0.75f );
	vec3 purple( 0.75f, 0.15f, 0.75f );
	vec3 white(  0.75f, 0.75f, 0.75f );
	vec3 black( 0.0f, 0.0f, 0.0f);

	triangles.clear();
	triangles.reserve( 5*2*3 );

	// ---------------------------------------------------------------------------
	// Room

	float L = 555;			// Length of Cornell Box side.

	vec3 A(L,0,0);
	vec3 B(0,0,0);
	vec3 C(L,0,L);
	vec3 D(0,0,L);

	vec3 E(L,L,0);
	vec3 F(0,L,0);
	vec3 G(L,L,L);
	vec3 H(0,L,L);

	int height = 800;
	int width = 800;

	vec3** texture = getTexture("concrete.bmp");
	vec3** bump = getTexture("bump3.bmp");

	std::cout << sizeof(texture) << std::endl;
	// Floor:
	triangles.push_back( Triangle( C, B, A, glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height), glm::normalize(D-C),green, texture, bump) );
	triangles.push_back( Triangle( C, D, B, glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height),glm::normalize(D-C), green,texture, bump) );

	// Left wall
	triangles.push_back( Triangle( A, E, C, glm::vec2(width, 0), glm::vec2(width, height),glm::vec2(0, 0),glm::normalize(A-C), purple, texture, bump) );
	triangles.push_back( Triangle( C, E, G, glm::vec2(0, 0),  glm::vec2(width, height),glm::vec2(0, height),glm::normalize(A-C),purple, texture, bump) );

	// Right wall
	triangles.push_back( Triangle( F, B, D,  glm::vec2(width, height),glm::vec2(0, height), glm::vec2(0, 0),glm::normalize(H-D), yellow, texture, bump) );
	triangles.push_back( Triangle( H, F, D, glm::vec2(width, 0), glm::vec2(width, height), glm::vec2(0, 0),glm::normalize(H-D),yellow, texture, bump) );

	// Ceiling
	triangles.push_back( Triangle( E, F, G,glm::vec2(0, height),  glm::vec2(width, height), glm::vec2(0, 0),glm::normalize(H-G), cyan, texture, bump) );
	triangles.push_back( Triangle( F, H, G, glm::vec2(width, height), glm::vec2(width, 0), glm::vec2(0, 0),glm::normalize(H-G), cyan, texture, bump) );

	// Back wall
	triangles.push_back( Triangle( G, D, C, glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height),glm::normalize(H-G), white, texture, bump) );
	triangles.push_back( Triangle( G, H, D, glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height),glm::normalize(H-G),white, texture, bump) );

	// ---------------------------------------------------------------------------
	// Short block

	A = vec3(290,0,114);
	B = vec3(130,0, 65);
	C = vec3(240,0,272);
	D = vec3( 82,0,225);

	E = vec3(290,165,114);
	F = vec3(130,165, 65);
	G = vec3(240,165,272);
	H = vec3( 82,165,225);

	// Front
	triangles.push_back( Triangle(E,B,A, glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height),glm::normalize(F-E), red, texture, bump) );
	triangles.push_back( Triangle(E,F,B, glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height),glm::normalize(F-E), red, texture, bump) );

	// Front
	triangles.push_back( Triangle(F,D,B,glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height),glm::normalize(H-F),red, texture, bump) );
	triangles.push_back( Triangle(F,H,D,glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height),glm::normalize(H-F), red, texture, bump) );

	// BACK
	triangles.push_back( Triangle(H,C,D,glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height),glm::normalize(G-H),red, texture, bump) );
	triangles.push_back( Triangle(H,G,C,glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height),glm::normalize(G-H), red, texture, bump) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height),glm::normalize(A-E),red, texture, bump) );
	triangles.push_back( Triangle(E,A,C,glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height), glm::normalize(A-E), red, texture, bump) );

	// TOP
	triangles.push_back( Triangle(G,F,E,glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height),glm::normalize(H-G), red, texture, bump) );
	triangles.push_back( Triangle(G,H,F, glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height), glm::normalize(H-G), red, texture, bump) );

	// ---------------------------------------------------------------------------
	// Tall block

	A = vec3(423,0,247);
	B = vec3(265,0,296);
	C = vec3(472,0,406);
	D = vec3(314,0,456);

	E = vec3(423,330,247);
	F = vec3(265,330,296);
	G = vec3(472,330,406);
	H = vec3(314,330,456);

	// Front
	triangles.push_back( Triangle(E,B,A,glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height), glm::normalize(F-E), blue, texture, bump) );
	triangles.push_back( Triangle(E,F,B,glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height), glm::normalize(F-E), blue, texture, bump) );

	// Front
	triangles.push_back( Triangle(F,D,B,glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height), glm::normalize(H-F), blue, texture, bump) );
	triangles.push_back( Triangle(F,H,D,glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height), glm::normalize(H-F), blue, texture, bump) );

	// BACK
	triangles.push_back( Triangle(H,C,D,glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height), glm::normalize(G-H), blue, texture, bump) );
	triangles.push_back( Triangle(H,G,C,glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height), glm::normalize(G-H), blue, texture, bump) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height), glm::normalize(A-E), blue, texture, bump) );
	triangles.push_back( Triangle(E,A,C,glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height), glm::normalize(A-E), blue, texture, bump) );

	// TOP
	triangles.push_back( Triangle(G,F,E,glm::vec2(0, 0), glm::vec2(width, height), glm::vec2(0, height), glm::normalize(H-G), blue, texture, bump) );
	triangles.push_back( Triangle(G,H,F,glm::vec2(0, 0), glm::vec2(width, 0), glm::vec2(width, height), glm::normalize(H-G), blue, texture, bump) );


	// ----------------------------------------------
	// Scale to the volume [-1,1]^3

	for( size_t i=0; i<triangles.size(); ++i )
	{
		triangles[i].v0 *= 2/L;
		triangles[i].v1 *= 2/L;
		triangles[i].v2 *= 2/L;

		triangles[i].v0 -= vec3(1,1,1);
		triangles[i].v1 -= vec3(1,1,1);
		triangles[i].v2 -= vec3(1,1,1);

		triangles[i].v0.x *= -1;
		triangles[i].v1.x *= -1;
		triangles[i].v2.x *= -1;

		triangles[i].v0.y *= -1;
		triangles[i].v1.y *= -1;
		triangles[i].v2.y *= -1;

		triangles[i].ComputeNormal();
	}
}

#endif

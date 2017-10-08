#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>
#include <memory>
#include <thread>
#include <atomic>

#include <SDL.h>

#include "matrix.h"
#include "math.h"
#include "AAB.h"
#include "Sphere.h"
#include "Plane.h"
#include "Shape.h"
#include "Triangle.h"

using byte = uint8_t;

Vector3f rayTrace(Vector3f rayOrigin, Vector3f rayDirection, int depth);
void setPixel(SDL_Surface* surface);
bool drawPixel(SDL_Surface* surface, int i, int j, const Vector3f& col);

//need aspect ratio, as image might not be a square
constexpr int height = 640, width = 840;
constexpr float aspectRatio = width / float(height);
constexpr int numShapes = 5, maxDepth = 5;
static bool hardShadows = true, recurseReflect = false; 
static float FOV = 90.0f;

std::vector<Shape*> shapes(numShapes);
//area light from the box
static std::unique_ptr<AAB> areaLight = std::make_unique<AAB>(Vector3f(-2.5f, 20.f, -2.5f), Vector3f(9.f, 0.1f, 9.f));

//img will represent the view plane
static std::vector<std::vector<Vector3f> > img(height, std::vector<Vector3f>(width));

bool drawPixel(SDL_Surface* surface, int x, int y, const Vector3f& col) {
	if (y < 0 || y >= height || x < 0 || x >= width)
		return false;

	uint32_t colorSDL = SDL_MapRGB(surface->format, col.r, col.g, col.b);

	uint32_t* bufp;
	bufp = (uint32_t*)surface->pixels + y * surface->pitch / 4 + x;
	*bufp = colorSDL;
	bufp += surface->pitch / 4;

	return true;
}

void setPixel(SDL_Surface* surface) {
    float scale = tan(M_RAD(FOV) / 2); 
    Vector3f rayOrigin = Vector3f(0,0,0);
    std::vector<std::thread> threads;

#pragma omp parallel for schedule(dynamic,1) collapse(2)
//#pragma omp parallel for
    for(int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
            //normalize the pixel positions to the range [0,1], +0.5 so the ray passes the center of the pixel
            float pixNormX = (x + 0.5f)/width, pixNormY = (y+0.5f)/height;
            //remap coordinates to range [-1, 1] and reverse direction of y-axis and scale by FOV
            float camX = (2*pixNormX-1)*aspectRatio*scale, camY = (1 - 2*pixNormY)*scale;

            //since point lies on the image plane, which is 1 unit away from the camera's origin(0,0,0) the z-axis is set to -1
            Vector3f rayDirection = {camX, camY, -1};
            rayDirection = m_normalize(rayDirection);

			Vector3f color = rayTrace(rayOrigin, rayDirection, 0);
			img[y][x] = color;
            color = Vector3f(byte(Math::clamp(color.r*255, 1, 255)) 
                            ,byte(Math::clamp(color.g*255, 1, 255)) 
                            ,byte(Math::clamp(color.b*255, 1, 255)));
			drawPixel(surface, x, y, color);
		}
	}
}

inline int getShapeIntersection(const std::vector<Shape*>& shapes, const Vector3f& origin,const Vector3f& direction, float& t0, float& minT) {
    int shapeHit = -1;
    for(int k = 0; k < shapes.size(); ++k) {
        bool doesIntersect = shapes[k]->doesIntersect(origin, direction,t0);
        if(doesIntersect && t0 < minT) {
            minT = t0;
            shapeHit = k;
        }
    }
    return shapeHit;
}


Vector3f rayTrace(Vector3f rayOrigin, Vector3f rayDirection, int depth) {
        float minT = M_INFINITE, t0 = 0.0f;
        int shapeHit = getShapeIntersection(shapes, rayOrigin, rayDirection, t0, minT);
        Vector3f pixelColor, reflectColor;

        if(shapeHit != -1) {
            Vector3f p0 = rayOrigin + (minT * rayDirection); //point of intersection

            //Light properties
            Vector3f areaLightCenter = Vector3f(areaLight->position.x + (areaLight->size.x / 2),
                                                areaLight->position.y + (areaLight->size.y / 2),
                                                areaLight->position.z + (areaLight->size.z / 2));
            Vector3f lightIntensity = Vector3f(1.f); 

            Vector3f diffuseColor   = Vector3f(0.f);
            Vector3f specularColor  = Vector3f(0.f);
            int shininess = 0;

            //ambient lighting
            Vector3f ambient = shapes[shapeHit]->color * Vector3f(0.07f); //the "color" of the shadow

            Vector3f hitNormal  = m_normalize(shapes[shapeHit]->getNormal(p0, shininess, diffuseColor, specularColor));

            //Diffuse lighting 
            //Phong reflection model - https://en.wikipedia.org/wiki/Phong_reflection_model
            Vector3f lightRay   = m_normalize(areaLightCenter - p0); //this is the reflective ray point towards the light source
            Vector3f diffuse    = diffuseColor * lightIntensity * std::max(0.0f, m_dot(lightRay, hitNormal));

            //Specular lighting
            Vector3f reflection = m_normalize(2*m_dot(lightRay,hitNormal)*hitNormal - lightRay); 
            float maxCalc = std::max(0.0f, m_dot(reflection, m_normalize(rayOrigin-p0)));
            Vector3f specular = specularColor * lightIntensity * pow(maxCalc, shininess);

            if(recurseReflect) {
                if ((depth < maxDepth) && (shininess > 0)) {
                    Vector3f reflectionRayDirection = rayDirection - 2 * m_dot(rayDirection, hitNormal) * hitNormal;
                    Vector3f reflectionRayOrigin = p0 + (hitNormal * M_EPSILON);
                    reflectColor = reflectColor + 0.05f * rayTrace(reflectionRayOrigin, reflectionRayDirection, depth + 1);
                    return pixelColor = diffuse + specular + reflectColor;
                }
            }

            //emit a "shadow ray" from the interception and check if it hits a shape, it is a shadow if it hits
            int lightShapeHit = -1;
            if(hardShadows) {
                //normal is multiplied by epsilon so we dont get a shadow ray that hits itself
                lightShapeHit = getShapeIntersection(shapes, p0 + (M_EPSILON*hitNormal),lightRay, t0, minT);
                pixelColor = lightShapeHit != -1 ? ambient : diffuse + specular;
            }
            else {
                float sample = 9, totalRays = sample * sample, raysHit = 1.f;
                float softIncrement = areaLight->size.x / sample;
                for(float m = 0; m < areaLight->size.x; m+= softIncrement) {
                    for(float n = 0; n < areaLight->size.z; n+= softIncrement) {
                        float t0s = 0.f, minTs = M_INFINITE;

                        Vector3f areaLightPos = Vector3f(m, areaLight->position.y, n);
                        lightRay = m_normalize(areaLightPos - p0);

                        lightShapeHit = getShapeIntersection(shapes, p0 + (M_EPSILON*hitNormal),lightRay, t0s, minTs);
                        if(lightShapeHit != -1)
                            raysHit = raysHit - (1 /(float)totalRays);
                    }
                }
                pixelColor = Vector3f((raysHit) * (diffuse + specular));
            }
        }
        else {
            pixelColor = Vector3f(1.f); //ray did not hit hit anything, set color to white
        }
        return pixelColor;
}

void saveToFile(std::string fileName, const std::vector<std::vector<Vector3f> >& img) {
    std::ofstream file(fileName, std::ios::out | std::ios::binary);
    file << "P6\n" << width << " " << height << "\n255\n";
    for(int y = 0; y < height; ++y) {
        for(int x = 0; x < width; ++x) {
            Vector3f pixel = img[y][x];
            file << byte(Math::clamp(pixel.r*255, 1, 255)) 
                 << byte(Math::clamp(pixel.g*255, 1, 255)) 
                 << byte(Math::clamp(pixel.b*255, 1, 255));
        }
    }
    file.close();
}
void createShapes() {

    /*
	shapes.push_back(new Plane(Vector3f(0, -4, 0), Vector3f(0, 1, 0), Vector3f(0.3f)));
	//======================= SPHERES =======================
	shapes.push_back(new Sphere(Vector3f(0, 0, -20), 4, Vector3f(1, 0.32f, 0.36f)));
	shapes.push_back(new Sphere(Vector3f(5, -1, -15), 2, Vector3f(0.9f, 0.76f, 0.46f)));
	shapes.push_back(new Sphere(Vector3f(5, 0, -25), 3, Vector3f(0.65f, 0.77f, 0.97f)));
	shapes.push_back(new Sphere(Vector3f(-5.5f, 0, -15), 3, Vector3f(0.9f)));
    */
	shapes[0] = (new Plane(Vector3f(0, -4, 0), Vector3f(0, 1, 0), Vector3f(0.3f)));
	//======================= SPHERES =======================
	shapes[1]=(new Sphere(Vector3f(0, 0, -20), 4, Vector3f(1, 0.32f, 0.36f)));
	shapes[2]=(new Sphere(Vector3f(5, -1, -15), 2, Vector3f(0.9f, 0.76f, 0.46f)));
	shapes[3]=(new Sphere(Vector3f(5, 0, -25), 3, Vector3f(0.65f, 0.77f, 0.97f)));
	shapes[4] = (new Sphere(Vector3f(-5.5f, 0, -15), 3, Vector3f(0.9f)));
}

int main() {
    /*
    shapes[0]= new Plane(Vector3f(0,-4,0),Vector3f(0,1,0),Vector3f(0.2f,0.2,0.2)); // Dark grey floor
    shapes[1]= new Sphere(Vector3f(0,0,-20),4,Vector3f(1,0.32,0.36)); //Red
    shapes[2]= new Sphere(Vector3f(5,-1,-15),2,Vector3f(0.9,0.76,0.46)); //Yellow
    shapes[3]= new Sphere(Vector3f(5,0,-25),3,Vector3f(0.65,0.77,0.97)); //Light blue 
    shapes[4]= new Sphere(Vector3f(-5.5,0,-15),3,Vector3f(0.9,0.9,0.9)); //Light grey
    */

    SDL_Window* window = nullptr;
    SDL_Surface* surface = nullptr;

	window = SDL_CreateWindow("Raytracing", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height, SDL_WINDOW_SHOWN);
    surface = SDL_GetWindowSurface(window);

    std::cout << "Press s to toggle hard/soft shadows\n(Soft shadows will be laggy)\n\n" << "Press space to toggle between teapot and spheres\n\n" << "Press r to toggle the reflection\n\n" << "Use up and down to change the FOV\n" << std::endl;;
    bool done = false;
    SDL_Event event;
    while(!done) {
        //resets the surface
        SDL_FillRect(surface, NULL, SDL_MapRGB(surface->format, 0, 0, 0));

        createShapes();
		setPixel(surface);
		SDL_UpdateWindowSurface(window);

        SDL_PollEvent(&event);

        if(event.type == SDL_QUIT) done = true;
		else if (event.type == SDL_KEYDOWN) {
			switch (event.key.keysym.sym) {
            case SDLK_ESCAPE:
                done = true;
                break;
			case SDLK_UP:
				if (FOV > 20.0f) FOV -= 10.0f;
                std::cout << "Decreasing the FOV\n";
				break;

			case SDLK_DOWN:
				if (FOV < 150.0f) FOV += 10.0f;
                std::cout << "Increasing the FOV\n";
				break;

			case SDLK_s:
				hardShadows = !hardShadows;
                std::cout << "Switching the shadows, please wait...\n";
				break;

			case SDLK_r:
				recurseReflect = !recurseReflect;
                std::cout << "Switching the reflections, please wait...\n";
				break;

			default:
				break;
			}
		}

    }

    SDL_DestroyWindow(window);
    SDL_Quit();

    for(int i = 0; i < shapes.size(); ++i) {
        delete shapes[i];
    }

    //saveToFile("test.ppm", img);
    return 0;
}

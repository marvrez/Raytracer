#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdint>
#include <vector>

#include "matrix.h"
#include "math.h"
#include "Sphere.h"
#include "Plane.h"
#include "Shape.h"
#include "Triangle.h"

using byte = uint8_t;

//need aspect ratio, as image might not be a square
constexpr int height = 640, width = 840;
constexpr float aspectRatio = width / float(height);
constexpr int numShapes = 6;

int main() {
    Shape* shapes[numShapes];
    shapes[0]= new Plane(Vector3f(0,-4,0),Vector3f(0,1,0),Vector3f(0.2,0.2,0.2)); // Dark grey floor
    shapes[1]= new Sphere(Vector3f(0,0,-20),4,Vector3f(1,0.32,0.36)); //Red
    shapes[2]= new Sphere(Vector3f(5,-1,-15),2,Vector3f(0.9,0.76,0.46)); //Yellow
    shapes[3]= new Sphere(Vector3f(5,0,-25),3,Vector3f(0.65,0.77,0.97)); //Light blue 
    shapes[4]= new Sphere(Vector3f(-5.5,0,-15),3,Vector3f(0.9,0.9,0.9)); //Light grey
    shapes[5]= new Triangle(Vector3f(0,3,-5), Vector3f(-1.9,1.2,-5), Vector3f(1.9, 1.2, -5), Vector3f(0.8,0,1)); //Purple

    //img will represent the view plane
    std::vector<std::vector<Vector3f> > img(height, std::vector<Vector3f>(width));

    //create field of view for cam at 30 degrees(standard for video games)
    float scale = tan(M_RAD(90.f) / 2); //30 degrees field of view
    Vector3f rayOrigin = Vector3f(0,0,0);

    for(int i = 0; i < height; ++i) {
        for(int j = 0; j < width; ++j) {
            //normalize the pixel positions to the range [0,1], +0.5 so the ray passes the center of the pixel
            float pixNormX = (j + 0.5f)/width, pixNormY = (i+ 0.5f)/height;
            //remap coordinates to range [-1, 1] and reverse direction of y-axis and scale by FOV
            //Points are now in the camera space
            float camX = (2*pixNormX-1)*aspectRatio*scale, camY = (1 - 2*pixNormY)*scale;

            //since point lies on the image plane, which is 1 unit away from the camera's origin(0,0,0) the z-axis is set to -1
            Vector3f rayDirection = {camX, camY, -1};
            rayDirection = m_normalize(rayDirection);

            float minT = M_INFINITE, t0 = 0.0f;
            int shapeHit = -1;

            for(int k = 0; k < numShapes; ++k) {
                bool doesIntersect = shapes[k]->doesIntersect(rayOrigin, rayDirection,t0);
                if(doesIntersect && t0 < minT) {
                    minT = t0;
                    shapeHit = k;
                }
                img[i][j] = shapeHit == -1 ? Vector3f(1) : shapes[shapeHit]->color;
            }
        }
    }

    std::ofstream file("./test.ppm", std::ios::out | std::ios::binary);
    file << "P6\n" << width << " " << height << "\n255\n";
    for(int y = 0; y < height; ++y) {
        for(int x = 0; x < width; ++x) {
            Vector3f vec = img[y][x];
            file << byte(Math::clamp(vec.x*255, 1, 255)) 
                 << byte(Math::clamp(vec.y*255, 1, 255)) 
                 << byte(Math::clamp(vec.z*255, 1, 255));
        }
    }
    file.close();

    return 0;
}

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdint>
#include <vector>

#include "matrix.h"
#include "math.h"
#include "Sphere.h"

using byte = uint8_t;

//need aspect ratio, as image might not be a square
constexpr int height = 640, width = 840;
constexpr float aspectRatio = width / float(height);
constexpr int numSpheres = 5;

int main() {
    Sphere sphereList[numSpheres];
    sphereList[0]= Sphere(Vector3f(0,0,-20),4,Vector3f(1.f, 0.32f, 0.36f)); //Red sphere
    sphereList[1]= Sphere(Vector3f(5,-1,-15),2,Vector3f(0.9f,0.76f,0.46f)); //Yellow sphere
    sphereList[2]= Sphere(Vector3f(5,0,-25),3,Vector3f(0.65f,0.77f,0.97f)); //Light blue
    sphereList[3]= Sphere(Vector3f(-5.5,0,-15),3,Vector3f(0.80f,0.80f,0.80f)); //Light grey
    sphereList[4]= Sphere(Vector3f(0,-10004,-20),10000,Vector3f(0.2f,0.2f,0.2f)); //Dark grey(using a large sphere as floor for now)

    //img will represent the view plane
    std::vector<std::vector<Vector3f> > img(height, std::vector<Vector3f>(width));

    //create field of view for cam at 30 degrees(standard for video games)
    float scale = tan(M_RAD(30.f) / 2); //30 degrees field of view
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
            int sphereHit = -1;

            for(int k = 0; k < numSpheres; ++k) {
                bool doesIntersect = sphereList[k].doesIntersect(t0, rayOrigin, rayDirection);
                //save nearest sphere
                if(doesIntersect && t0 < minT) {
                    minT = t0;
                    sphereHit = k;
                }
                img[i][j] = sphereHit == -1 ? Vector3f(1) : sphereList[sphereHit].color;
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

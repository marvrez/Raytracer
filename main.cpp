#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdint>
#include <vector>

#include "matrix.h"
#include "math.h"

using byte = uint8_t;

//need aspect ratio, as image might not be a square
constexpr int height = 640, width = 840, aspectRatio = width / height;

int main() {
    //img will represent the view plane
    Vector3f** img = new Vector3f*[height];
    for(int i = 0; i < height; ++i) 
        img[i] = new Vector3f[width];

    for(int i = 0; i < height; ++i) {
        for(int j = 0; j < width; ++j) {
            //normalize the pixel positions to the range [0,1], +0.5 so the ray passes the center of the pixel
            float pixNormX = (j + 0.5f)/width, pixNormY = (i+ 0.5f)/height;

            //remap coordinates to range [-1, 1] and reverse direction of y-axis
            float pixRemapX = (2*pixNormX-1) * aspectRatio, pixRemapY = 1 - 2*pixNormY;

            //create field of view for cam at 30 degrees(standard for video games)
            float camX = pixRemapX * tan(M_RAD(30.f)/2), camY = pixRemapY * tan(M_RAD(30.f) / 2);

            //since point lies on the image plane, which is 1 unit away from the camera's origin(0,0,0) the z-axis is set to -1
            Vector3f rayDirection = {camX, camY, 1};
            img[i][j] = m_normalize(rayDirection);
        }
    }

    std::ofstream file("./test.ppm", std::ios::out | std::ios::binary);
    file << "P6\n" << width << " " << height << "\n255\n";
    for(int y = 0; y < height; ++y) {
        for(int x = 0; x < width; ++x) {
            Vector3f vec = img[y][x];
            file << byte(std::min(1.f, vec.x * 255.f)) 
                 << byte(std::min(1.f, vec.y * 255.f)) 
                 << byte(std::min(1.f, vec.z * 255.f));
        }
    }
    file.close();

    return 0;
}

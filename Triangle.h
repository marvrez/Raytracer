#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "matrix.h"
#include "Shape.h"

//based off barycentric coordinates
struct Triangle : public Shape {
    Triangle();
    Triangle(Vector3f a, Vector3f b, Vector3f c, Vector3f color);

    bool doesIntersect(Vector3f rayOrigin, Vector3f rayDirection, float& t);

    Vector3f a,b,c;
};

#endif

#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "matrix.h"
#include "Shape.h"

//based off barycentric coordinates
struct Triangle : public Shape {
    Triangle();
    Triangle(Vector3f a, Vector3f b, Vector3f c,
             Vector3f u, Vector3f v, Vector3f w,
             Vector3f color);

    virtual bool doesIntersect(Vector3f rayOrigin, Vector3f rayDirection, float& t);
    virtual Vector3f getNormal(Vector3f p0, int& shininess,
                               Vector3f& diffuseColor, Vector3f& specularColor);

    Vector3f a,b,c;
    Vector3f u,v,w; //normal vectors to a,b and c respectively
};

#endif

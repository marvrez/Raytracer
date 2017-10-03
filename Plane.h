#ifndef PLANE_H
#define PLANE_H

#include "matrix.h"
#include "Shape.h"

struct Plane : public Shape {
    Plane();
    Plane(Vector3f position, Vector3f normal, Vector3f color);

    virtual bool doesIntersect(Vector3f rayOrigin, Vector3f rayDirection, float& t);

    virtual Vector3f getNormal(Vector3f p0, int& shininess,
                               Vector3f& diffuseColor, Vector3f& specularColor);

    Vector3f normal;
    Vector3f position;
};

#endif

#ifndef SPHERE_H
#define SPHERE_H

#include "math.h"
#include "Shape.h"

struct Sphere : public Shape {
    Sphere();
    Sphere(Vector3f centre, float radius, Vector3f color);

    virtual bool doesIntersect(Vector3f rayOrigin, Vector3f rayDirection, float& t);

    virtual Vector3f getNormal(Vector3f p0, int& shininess,
                               Vector3f& diffuseColor, Vector3f& specularColor);

    Vector3f centre;
    float radius, radius2; //cache radius^2 so we save some time calculating it over and over again
};


#endif

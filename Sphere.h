#ifndef SPHERE_H
#define SPHERE_H

#include "math.h"

struct Sphere {
    Sphere();
    Sphere(Vector3f centre, float radius, Vector3f color);

    bool doesIntersect(float& t, Vector3f rayOrigin, Vector3f rayDirection);

    Vector3f centre;
    Vector3f color;
    float radius, radius2; //cache radius^2 so we save some time calculating it over and over again
};


#endif

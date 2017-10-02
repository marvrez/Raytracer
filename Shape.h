#ifndef SHAPE_H
#define SHAPE_H

#include "matrix.h"

struct Shape {
    Shape();
    Shape(Vector3f position, Vector3f color);

    virtual bool doesIntersect(Vector3f rayOrigin, Vector3f rayDirection, float& t) = 0;

    Vector3f position;
    Vector3f color;
};

#endif

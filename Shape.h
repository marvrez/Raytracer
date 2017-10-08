#ifndef SHAPE_H
#define SHAPE_H

#include "matrix.h"

struct Shape {
    Shape();
    Shape(Vector3f color);
    virtual ~Shape() = default;

    virtual bool doesIntersect(Vector3f rayOrigin, Vector3f rayDirection, float& t) = 0;

    //Diffuse color => Color that object reflects when illuminated by white light.
    //Specular color => color of highlights on the shiny surface. 
    //(|Specular col.| == |Diffuse col.|) => less shiny surface.
    virtual Vector3f getNormal(Vector3f p0, int& shininess,
                               Vector3f& diffuseColor, Vector3f& specularColor) = 0;

    Vector3f color;
};

#endif

#include "Plane.h"

#include "math.h"

Plane::Plane() 
    : Shape(), normal(0)
{
}

Plane::Plane(Vector3f position, Vector3f normal, Vector3f color) 
    : Shape(color), position(position), normal(normal)
{
}

bool Plane::doesIntersect(Vector3f rayOrigin, Vector3f rayDirection, float& t) {
    float denom = m_dot(rayDirection, normal);
    if(abs(denom) < M_EPSILON) return false; //ray is paralell with plane
    float result = m_dot(position - rayOrigin, normal) / denom;
    return (t = result) >= 0;
}

Vector3f Plane::getNormal(Vector3f p0, int& shininess,
                          Vector3f& diffuseColor, Vector3f& specularColor)
{
    shininess = 0;
    diffuseColor = Vector3f(0.8f, 0.8f, 0.8f);
    specularColor = this->color;
    return this->normal;
}

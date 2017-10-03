#include "Triangle.h"

Triangle::Triangle() 
    : Shape(), a(0), b(0), c(0)
{
}

Triangle::Triangle(Vector3f a, Vector3f b, Vector3f c, Vector3f color) 
    : Shape(color), a(a), b(b), c(c)
{
}

//https://en.wikipedia.org/wiki/Möller-Trumbore_intersecton_algorithm
bool Triangle::doesIntersect(Vector3f rayOrigin, Vector3f rayDirection, float& t) {
    Vector3f e1 = b-a, e2 = c-a;
    float u = m_dot(rayOrigin-a, m_cross(rayDirection, e2)) / m_dot(e1, m_cross(rayDirection, e2));
    float v = m_dot(rayDirection, m_cross(rayOrigin-a, e1)) / m_dot(e1, m_cross(rayDirection, e2));

    if(u < 0 || u > 1) return false;
    else if (v < 0 || u+v > 1) return false;
    else {
        t = m_dot(e2, m_cross(rayOrigin-a, e1)) / m_dot(e1, m_cross(rayDirection,e2));
        return true;
    }
}

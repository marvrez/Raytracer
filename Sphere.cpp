#include "Sphere.h"

Sphere::Sphere() 
    : Shape(), centre(0), radius(0), radius2(0)
{
}

Sphere::Sphere(Vector3f centre, float radius, Vector3f color) 
    : Shape(color),centre(centre), radius(radius), radius2(radius*radius)
{
}

bool Sphere::doesIntersect(Vector3f rayOrigin, Vector3f rayDirection, float& t) {
    //geometric solution for finding the intersection of a ray and sphere
    Vector3f L = this->centre - rayOrigin;
    float tca = m_dot(L, rayDirection); 
    if(tca < 0) return false; //L and rayDirection points in opposite directions
    float d2 = m_dot(L,L) - tca*tca;
    if(d2 > this->radius2) return false; //ray misses the sphere
    float thc = sqrt(radius2-d2);
    t = tca - thc; //save point of intersection
    return true;
}

Vector3f Sphere::getNormal(Vector3f p0, int& shininess,
                           Vector3f& diffuseColor, Vector3f& specularColor)
{
    shininess = 128; //adds a glow effect to the spheres
    diffuseColor = this->color;
    specularColor = Vector3f(0.7f, 0.7f, 0.7f);
    return (p0 - this->centre);
}

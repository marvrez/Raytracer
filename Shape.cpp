#include "Shape.h"

Shape::Shape() 
    : position(0), color(0)
{
}

Shape::Shape(Vector3f position, Vector3f color) 
    : position(position), color(color)
{
}

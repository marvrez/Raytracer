#ifndef AAB_H
#define AAB_H

#include "matrix.h"

//AAB == Axis-aligned-box
//We will be using this class to create sof shadows by essentially trying to shoot shadow rays
//towards an AAB rather than a point as we did for hard shadows
struct AAB {
    AAB();
    AAB(Vector3f position, Vector3f size);
    
    Vector3f position, size;
};

#endif

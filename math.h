#ifndef MATH_H
#define MATH_H

#include "matrix.h"

#include <cmath>
#include <complex>
#include <float.h>

namespace Math {
//Defines
#define M_PI            3.1415926535897932384626433832795f
#define M_RAD(val)        ((val)*0.0174532925199432957692369076848f)
#define M_DEG(val)        ((val)*57.295779513082320876798154814105f)
#define M_LARGE_EPSILON   1e-2f
#define M_EPSILON         1e-4f
#define M_TINY_EPSILON    1e-5f
#define M_INFINITE      3.4e38f
#define M_EULER_NUMBER  2.71828182845904523536f
#define M_SQRT2         1.4142135623730950488016887242097f
#define M_SQRT3         1.7320508075688772935274463415059f

//(kinda optimized) math functions
inline float fSinCosRecurse (float ret, float rad, const float rad2, float radpow, float fact, const int prec, int cur) {
    radpow *= -rad2;
    fact   *= cur++;
    fact   *= cur++;

    ret += radpow / fact;

    if(cur < prec)
        return fSinCosRecurse(ret, rad, rad2, radpow, fact, prec, cur);
    else
        return ret;
}

inline float fSin (float rad, int precision = 4) {
    while(rad >  M_PI) rad -= (M_PI*2.0f);
    while(rad < -M_PI) rad += (M_PI*2.0f);

    const float rad2 = (rad*rad);
    float ret = rad;

    return fSinCosRecurse(ret, rad, rad2, rad, 1.0f, (precision << 1) + 1, 2);
}

inline float fCos (float rad, int precision = 4) {
    while(rad >  M_PI) rad -= (M_PI*2.0f);
    while(rad < -M_PI) rad += (M_PI*2.0f);

    const float rad2 = (rad*rad);
    float ret = 1.0f;

    return fSinCosRecurse(ret, rad, rad2, 1.0f, 1.0f, (precision << 1) + 1, 1);
}

inline bool fEqual(float f1, float f2, float precision = M_LARGE_EPSILON) {
    return fabsf(f1-f2) <= precision;
}

inline int factorial(int n) {
    int res = 1;
    for (int i = 1; i <= n; ++i)
        res *= i;
    return res;
}

inline int digitSum(int n) {
    return n ? (n%10) + digitSum(n/10) : 0;
}

inline unsigned long long combination(unsigned long long n, unsigned long long k) {
    if (k > n) return 0;

    if (k > n/2) k = n-k;

    unsigned long long r = 1;
    for (unsigned long long d = 0; d < k; ++d) {
        r *= (n - d);
        r /= (d + 1);
    }
    return r;
}

inline std::pair<std::complex<float>, std::complex<float> > solveQuadratic(float a, float b, float c) {
    std::complex<float> x1 = 0, x2 = 0;
    float determinant = b*b - (4 * a*c);

    if (determinant > 0) {
        x1 = (-b + sqrt(determinant)) / (2*a);
        x2 = (-b - sqrt(determinant)) / (2*a);
    }
    else if (determinant == 0) {
        x1 = (-b + sqrt(determinant)) / (2*a);
        x2 = (-b + sqrt(determinant)) / (2*a);
    }
    else {
        x1.real(-b / (2*a));
        x1.imag(sqrt(-determinant) / (2*a));

        x2.real(-b / (2*a));
        x2.imag(-sqrt(-determinant) / (2*a));
    }

    return std::make_pair(x1, x2);
}

inline bool doIntersect(Vector2f p1, Vector2f q1, Vector2f p2, Vector2f q2) {
    static auto det = [](const Vector2f& u, const Vector2f& v) {return u.x*v.y - u.y*v.x;};
    return (det(p2-p1, q1-p1)*det(p2-p1, q2-p1) < 0) && (det(q2-q1, p1-q1)*det(q2-q1, p2-q1) < 0);
}

inline Vector2f getIntersection(Vector2f p1, Vector2f q1, Vector2f p2, Vector2f q2) {
    if(doIntersect(p1,q1,p2,q2)) {
        if(p2.x - p1.x != 0 && q2.x - q1.x != 0) {
            float a1 = (p2.y - p1.y)/(p2.x - p1.x);
            float a2 = (q2.y - q1.y)/(q2.x - q1.x);
            float b1 = p1.y - a1*p1.x;
            float b2 = q1.y - a2*q1.x;
            float Ix = (b2-b1)/(a1-a2);
            float Iy = a1*Ix+b1;
            return Vector2f(Ix,Iy);
        }

        if(p2.x - p1.x == 0) {
            float a2 = (q2.y - q1.y)/(q2.x - q1.x);
            float b1 = p1.x;
            float b2 = q1.y - a2*q1.x;
            float Ix = b1;
            float Iy = a2*b1+b2;
            return Vector2f(Ix,Iy);
        }

        if(q2.x - q1.x == 0) {
            float a1 = (p2.y - p1.y)/(p2.x - p1.x);
            float b1 = p1.y - a1*p1.x;
            float b2 = q1.x;
            float Ix = b2;
            float Iy = a1*b2+b1;
            return Vector2f(Ix,Iy);
        }
    }

    return Vector2f();
}

inline float gaussianFunction(float maxVal, float wideness, float x) {
    return maxVal*exp(-pow(x,2)/(2*pow(wideness,2)));
}

inline float clamp(float val, float lower, float upper) {
    return std::max(lower, std::min(val, upper));
}

inline float fastInvSqrt(float x) {
    //https://en.wikipedia.org/wiki/Fast_inverse_square_root
    float xhalf = 0.5f * x;
    int i = *(int*)&x;          // Integer representation of float
    i = 0x5f3759df - (i >> 1);  // Initial guess
    x = *(float*)&i;            // Converted back to floating point
    x = x*(1.5f-(xhalf*x*x));   // One round of Newton-Raphson's method
    return x;
}



} //namespace Math;


#endif // MATH_H

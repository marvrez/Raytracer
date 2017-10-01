#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <ostream>
#include <stdexcept>

template <int n, typename T = float>
struct Vector {
    T data[n];
    T& operator[](int i) { return data[i]; }

    const T& operator[](int i) const { return data[i]; }
};

template <int rows, int cols, typename T = float>
struct Matrix {
    T data[rows*cols];

    T &operator[](int i) { return data[i]; }
    T& operator()(int row, int col) { return data[row + col*rows]; }

    const T& operator()(int row, int col) const { return data[row + col*rows]; }
    const T& operator[](int i) const { return data[i]; }

    float det() const;
};

//Most used vector types
using Vector2f = Vector<2>;
using Vector3f = Vector<3>;
using Vector4f = Vector<4>;

using Vector2i = Vector<2, int>;
using Vector3i = Vector<3, int >;
using Vector4i = Vector<4, int>;

using Vector2u = Vector<2, unsigned>;
using Vector3u = Vector<3, unsigned>;
using Vector4u = Vector<4, unsigned>;

//typical matrix types
using Matrix2f = Matrix<2,2>;
using Matrix3f = Matrix<3,3>;
using Matrix4f = Matrix<4,4>;

using Matrix2i = Matrix<2,2, int>;
using Matrix3i = Matrix<3,3, int >;
using Matrix4i = Matrix<4,4, int>;

using Matrix2u = Matrix<2,2, unsigned>;
using Matrix3u = Matrix<3,3, unsigned>;
using Matrix4u = Matrix<4,4, unsigned>;

///////////////// Vector Specializations /////////////////
// Specializations for n = 2, 3, 4 so that we can access each
// component of a vector through .a notation, as well as access
// truncated vectors.

template <typename T>
struct Vector<2, T> {
    union {
        T data[2];
        struct { T x, y; };
    };
    T& operator[](int i) { return data[i]; }
    const T& operator[](int i) const { return data[i]; }

    Vector() = default;
    explicit Vector(T num);
    Vector(T x, T y);
};

template <typename T>
struct Vector<3, T> {
    union {
        T data[3];
        struct { T x, y, z; };
        struct { T r, g, b; };
        Vector<2, T> xy;
    };
    T& operator[](int i) { return data[i]; }
    const T& operator[](int i) const { return data[i]; }

    Vector() = default;
    explicit Vector(T num);
    Vector(T x, T y, T z);
};

template <typename T>
struct Vector<4, T> {
    union {
        T data[4];
        struct { T x, y, z, w; };
        struct { T r, g, b, a; };
        Vector<3, T> xyz;
        Vector<3, T> rgb;
        Vector<2, T> xy;
    };
    T& operator[](int i) { return data[i]; }
    const T& operator[](int i) const { return data[i]; }

    Vector() = default;
    explicit Vector(T num);
    Vector(T x, T y, T z, T w);
};

///////////////// Matrix Specializations /////////////////
// The data layout for each matrix is interpreted as column
// order. For example
// | a11 a12 |
// | a21 a22 |
// is stored in memory as [a11, a21, a12, a22].
// The rows of the matrix can be accessed via .a1 or .a2.

template <typename T>
struct Matrix<2, 2, T> {
    union {
        T data[4];
        T mat[2][2];
        struct { T a11, a21, a12, a22; };
        struct { Vector<2,T> a1, a2; };
    };
    const T& operator[](int i) const { return data[i]; }
    T& operator[](int i) { return data[i]; }

    float det() const;
    Matrix inverse();

    Matrix() = default;
    Matrix(int i, float theta = 0); //rotate
    Matrix(T a0, T a1, T a2, T a3);
};

template <typename T>
struct Matrix<3, 3, T> {
    union {
        T data[9];
        T mat[3][3];
        struct { T a11, a21, a31, a12, a22, a32, a13, a23, a33; };
        struct { Vector<3,T> a1, a2, a3; };
    };
    T& operator[](int i) { return data[i]; }
    const T& operator[](int i) const { return data[i]; }

    float det() const;

    Matrix() = default;
    Matrix(int i, float theta = 0);
    Matrix(int i, float arg1, float arg2);
    Matrix(int i, const Vector<3,T>& v, float arg1);
    Matrix(T a0, T a1, T a2,
           T a3, T a4, T a5,
           T a6, T a7, T a8);
};

template <typename T>
struct Matrix<4, 4, T> {
    union {
        T data[16];
        T mat[4][4];
        struct { T a11, a21, a31, a41, a12, a22, a32, a42, a13, a23, a33, a43, a14, a24, a34, a44; };
        struct { Vector<4,T> a1, a2, a3, a4; };
    };
    const T& operator[](int i) const { return data[i]; }
    T& operator[](int i) { return data[i]; }

    float det() const;

    Matrix() = default;
    Matrix(int i, float theta = 0);
    Matrix(int i, float arg1, float arg2);
    Matrix(int i, const Vector<3,T>& p);
    Matrix(T a0, T a1, T a2, T a3,
           T a4, T a5, T a6, T a7,
           T a8, T a9, T a10, T a11,
           T a12, T a13, T a14, T a15);
};

//////////// Matrix and vector constructors ////////////

#define templ_t template<typename T>

templ_t Vector<2,T>::Vector(T num) : x(num), y(num)  { }
templ_t Vector<3,T>::Vector(T num) : x(num), y(num), z(num) { }
templ_t Vector<4,T>::Vector(T num) : x(num), y(num), z(num), w(num) { }

templ_t Vector<2,T>::Vector(T x, T y) : x(x), y(y) { }
templ_t Vector<3,T>::Vector(T x, T y, T z) : x(x), y(y), z(z) { }
templ_t Vector<4,T>::Vector(T x, T y, T z, T w) : x(x), y(y), z(z), w(w) { }

//TODO: use enum instead?
#define MATRIX_IDENTITY    0
#define MATRIX_ROTATION    1
#define MATRIX_TRANSLATION 2
#define MATRIX_ROTATION_X  3
#define MATRIX_ROTATION_Y  4
#define MATRIX_ROTATION_Z  5

//Clockwise rotations for Matrix2
#define MATRIX2_ROTATE_90  6
#define MATRIX2_ROTATE_180 7
#define MATRIX2_ROTATE_270 8

templ_t Matrix<2,2,T>::Matrix(int i, float theta) {
    switch(i) {
    case MATRIX_ROTATION: {
        float st = sin(theta), ct = cos(theta);
        data[0] = ct, data[1] = -st,
        data[2] = st, data[3] = ct;
    }
    break;
    case MATRIX2_ROTATE_90: {
        data[0] = 0; data[1] = -1;
        data[2] = 1; data[3] = 0;
    }
    break;
    case MATRIX2_ROTATE_180: {
        data[0] = -1; data[1] = 0;
        data[2] =  0;  data[3] = -1;
    }
    break;
    case MATRIX2_ROTATE_270: {
        data[0] = 0;  data[1] = 1;
        data[2] = -1; data[3] = 0;
    }
    break;
    case MATRIX_IDENTITY:
    default: {
        data[0] = 1; data[1] = 0,
        data[2] = 0; data[3] = 1;
    }
    break;
    }
}

templ_t Matrix<2,2,T>::Matrix(T a0, T a1, T a2, T a3) {
    data[0] = a0; data[1] = a1;
    data[2] = a2; data[3] = a3;
}

templ_t float Matrix<2,2,T>::det() const {
    return data[0]*data[3] - data[1]*data[2];
}

templ_t Matrix<2,2,T> Matrix<2,2,T>::inverse() {
    return (1 / this->det()) * Matrix<2,2,T>(data[3], -data[1],
                                            -data[2],  data[0]);
}

templ_t Matrix<3,3,T>::Matrix(int i, float theta) {
    switch(i){
        float ct, st;
        case MATRIX_ROTATION:
        case MATRIX_ROTATION_X: {
            st = sin(theta), ct = cos(theta);
            data[0] = 1; data[1] = 0;  data[2] = 0;
            data[3] = 0; data[4] = ct; data[5] = -st;
            data[6] = 0; data[7] = st; data[8] = ct;
        }
        break;
        case MATRIX_ROTATION_Y: {
            st = sin(theta), ct = cos(theta);
            data[0] = ct;  data[1] = 0;  data[2] = st;
            data[3] = 0;   data[4] = 1;  data[5] = 0;
            data[6] = -st; data[7] = 0;  data[8] = ct;
        }
        break;
        case MATRIX_ROTATION_Z: {
            st = sin(theta), ct = cos(theta);
            data[0] = ct; data[1] = -st;data[2] = 0;
            data[3] = st; data[4] = ct; data[5] = 0;
            data[6] = 0;  data[7] = 0;  data[8] = 1;
        }
        break;
        case MATRIX_IDENTITY:
        default:
            for(int i = 0; i < 9; i++)
                data[i] = i%4 == 0;
            break;
    }
}

templ_t Matrix<3,3,T>::Matrix(int i, float arg1, float arg2) {
    switch(i) {
        case MATRIX_ROTATION: {
            float st = sin(arg1), ct = cos(arg1),
            sp = sin(arg2), cp = cos(arg2);
            data[0] = ct;    data[1] = -st;   data[2] = 0;
            data[3] = st*cp; data[4] = ct*cp; data[5] = -sp;
            data[6] = st*sp; data[7] = ct*sp; data[8] = cp;
        }
        break;
    }
}

templ_t Matrix<3,3,T>::Matrix(int i, const Vector<3,T>& v, float arg1) {
    switch(i) {
        case MATRIX_ROTATION: {
            float len = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
            Matrix<3,3,T> what(0,    -v.z,  v.y,
                               v.z,     0, -v.x,
                              -v.y,   v.x,    0);
            Matrix<3,3,T> u = Matrix<3,3,T>(MATRIX_IDENTITY) + what/(len)*sin(arg1) + what*what/pow(len, 2)*(1 - cos(arg1));

            for(int i = 0 ; i< 9; ++i)
                data[i] = u.data[i];
        }
        break;
    }
}

templ_t Matrix<3,3,T>::Matrix(T a0, T a1, T a2,
                              T a3, T a4, T a5,
                              T a6, T a7, T a8)
{
    data[0] = a0; data[1] = a1; data[2] = a2;
    data[3] = a3; data[4] = a4; data[5] = a5;
    data[6] = a6; data[7] = a7; data[8] = a8;
}

templ_t float Matrix<3,3,T>::det() const {
    return data[0]*(data[4]*data[8] - data[5]*data[7]) -
           data[1]*(data[3]*data[8] - data[5]*data[6]) +
           data[2]*(data[3]*data[7] - data[4]*data[6]);
}

templ_t Matrix<4,4,T>::Matrix(int i, float theta) {
    switch(i){
        float st, ct;
        case MATRIX_ROTATION:
        case MATRIX_ROTATION_X: {
            st = sin(theta), ct = cos(theta);
            data[0] = 1; data[1] = 0; data[2] = 0; data[3] = 0;
            data[4] = 0; data[5] = ct; data[6] = -st; data[7] = 0;
            data[8] = 0; data[9] = st; data[10] = ct; data[11] = 0;
            data[12] = 0; data[13] = 0; data[14] = 0; data[15] = 1;
        }
        break;
        case MATRIX_ROTATION_Y: {
            st = sin(theta), ct = cos(theta);
            data[0] = ct; data[1] = 0; data[2] = st; data[3] = 0;
            data[4] = 0; data[5] = 1; data[6] = 0; data[7] = 0;
            data[8] = -st; data[9] = 0; data[10] = ct; data[11] = 0;
            data[12] = 0; data[13] = 0; data[14] = 0; data[15] = 1;
        }
        break;
        case MATRIX_ROTATION_Z: {
            st = sin(theta), ct = cos(theta);
            data[0] = ct; data[1] = -st; data[2] = 0; data[3] = 0;
            data[4] = st; data[5] = ct; data[6] = 0; data[7] = 0;
            data[8] = 0; data[9] = 0; data[10] = 1; data[11] = 0;
            data[12] = 0; data[13] = 0; data[14] = 0; data[15] = 1;
        }
        break;
        case MATRIX_IDENTITY:
        default:
            for(int i = 0; i < 16; ++i)
                data[i] = i%5 == 0;
        break;
    }
}

templ_t Matrix<4,4,T>::Matrix(int i, float arg1, float arg2) {
    switch(i){
        case MATRIX_ROTATION: {
            float st = sin(arg1), ct = cos(arg1),
            sp = sin(arg2), cp = cos(arg2);
            data[0] = ct; 	 data[1] = -st;   data[2] = 0;   data[3] = 0;
            data[4] = st*cp; data[5] = ct*cp; data[6] = -sp; data[7] = 0;
            data[8] = st*sp; data[9] = ct*sp; data[10] = cp; data[11] = 0;
            data[12] = 0; 	 data[13] = 0; 	  data[14] = 0;  data[15] = 1;
        }
        break;
    }
}

templ_t Matrix<4,4,T>::Matrix(int i, const Vector<3,T>& p) {
    switch(i){
        case MATRIX_TRANSLATION: {
            for(int i = 0; i< 16; ++i)
                data[i] = i%5 ==0;
            data[3] = p.x;
            data[7] = p.y;
            data[11] = p.z;
        }
        break;
    }
}

templ_t Matrix<4,4,T>::Matrix(T a0, T a1, T a2, T a3,
                              T a4, T a5, T a6, T a7,
                              T a8, T a9, T a10, T a11,
                              T a12, T a13, T a14, T a15)
{
    data[0] = a0;   data[1] = a1;   data[2] = a2;   data[3] = a3;
    data[4] = a4;   data[5] = a5;   data[6] = a6;   data[7] = a7;
    data[8] = a8;   data[9] = a9;  data[10] = a10; data[11] = a11;
    data[12] = a12; data[13] = a13; data[14] = a14; data[15] = a15;
}

templ_t float Matrix<4,4,T>::det() const {
    return data[0] * Matrix<3,3,T>(data[5],data[6],data[7],data[9],data[10],data[11],data[13],data[14],data[15]).det()
          -data[1] * Matrix<3,3,T>(data[4],data[6],data[7],data[8],data[10],data[11],data[12],data[14],data[15]).det()
          +data[2] * Matrix<3,3,T>(data[4],data[5],data[7],data[8],data[9], data[11],data[12],data[13],data[15]).det()
          -data[3] * Matrix<3,3,T>(data[4],data[5],data[6],data[8],data[9], data[10],data[12],data[13],data[14]).det();
}
#undef templ_t

///////////////// Matrix functions /////////////////

template<int r, int c, typename T>
float Matrix<r,c,T>::det() const {
    if(r != c) return 0; //matrix has no determinant

    if(r == 2 && c == 2) return data[0]*data[3] - data[1]*data[2];

    const auto minor = [&](int a, int b) {
        Matrix<r-1,c-1,T> mat;
        int offset_i = 0, offset_j = 0;
        for(int i = 0; i < r - 1; ++i) {
            if(i == b) offset_i++;
            for(int j = 0; j < c - 1; ++j) {
                if(j == a) offset_j++;
                mat[j + i*(c-1)] = (*this)(j + offset_j, i + offset_i);
            }
            offset_j = 0;
        }
        return mat;
    };

    //TODO: MAKE LAMBDAS _minor and _det!!!!!!
    float sum = 0;
    int factor = 1;
    for(int i = 0; i < c; i++) {
        sum += factor * data[i] * minor(i, 0).det();
        factor *= -1;
    }
    return sum;
}

template <int r, int c, typename T>
T* m_element(Matrix<r, c, T>* m, int row, int column) {
    return m->data + row + (column * r);
}

template <int r, int c, typename T>
Vector<r, T> m_getcolumn(const Matrix<r, c, T>& m, int column) {
    Vector<r, T> result;
    for (int i = 0; i < r; ++i)
        result.data[i] = m.data[i + (column * r)];
    return result;
}

template <int ra, int ca, int cb, typename T>
Matrix<ra, cb, T> operator*(Matrix<ra, ca> a, Matrix<ca, cb, T> b) {
    Matrix<ra, cb, T> result = {};
    float *entry = (float*)result.data;
    for (int col = 0; col < cb; col++)
    for (int row = 0; row < ra; row++) {
        for (int i = 0; i < ca; ++i) {
            float x = a.data[i * ra + row];
            float y = b.data[i + col * ca];
            *entry += x * y;
        }
        entry++;
    }
    return result;
}

template <int r, int c>
Matrix<r, c> operator+(const Matrix<r, c>& a, const Matrix<r, c>& b) {
    Matrix<r, c> result;
    for (int i = 0; i < r*c; ++i) result.data[i] = a.data[i] + b.data[i];
    return result;
}

template <int r, int c>
Matrix<r, c> operator-(const Matrix<r, c>& a, const Matrix<r, c>& b) {
    Matrix<r, c> result;
    for (int i = 0; i < r*c; ++i) result.data[i] = a.data[i] - b.data[i];
    return result;
}

template <int r, int c>
Matrix<r, c> operator*(const Matrix<r, c>& a, float s) {
    Matrix<r, c> result;
    for (int i = 0; i < r*c; ++i) result.data[i] = a.data[i] * s;
    return result;
}
template <int r, int c>
Matrix<r, c> operator*(float s, const Matrix<r, c>& a) {
    return a*s;
}

template <int r, int c>
Matrix<r, c> operator-(const Matrix<r, c>& a) {
    Matrix<r, c> result;
    for (int i = 0; i < r*c; ++i) result.data[i] = -a.data[i];
    return result;
}


template <int r, int c, typename T>
Matrix<c, r, T> m_transpose(Matrix<r, c, T> m) {
    Matrix<c, r, T> result = {};
    for (int row = 0; row < r; row++)
        for (int col = 0; col < c; col++) {
            T* a = m_element(&result, col, row);
            T* b = m_element(&m, row, col);
            *a = *b;
        }
    return result;
}

//TODO: m_invert() ?

///////////////// Vector functions /////////////////
// 	   a := vector of dimension N and type T
//     b := another vector of dimension N and type T
//     s := scalar of same type T
//
// The following operators are defined
//
//     -A        :: Component-wise negation of A
//     A+B = B+A :: B added component-wise to A
//     A-B       :: B subtracted component-wise from A
//     A*s = s*A :: s multiplied each component of A
//     A*B = B*A :: A multiplied component-wise by B
//     A/B       :: A divided component-wise by B
//     dot(A, B) :: The inner product of A and B
#define vector_template template <int n, typename T>
#define vec_t Vector<n, T>

vector_template vec_t  operator+(vec_t a, vec_t b)   { vec_t r; for (int i=0; i<n; ++i) r.data[i] = a.data[i]+b.data[i]; return r; }
vector_template vec_t  operator-(vec_t a, vec_t b)   { vec_t r; for (int i=0; i<n; ++i) r.data[i] = a.data[i]-b.data[i]; return r; }
vector_template vec_t  operator-(vec_t a)            { vec_t r; for (int i=0; i<n; ++i) r.data[i] = -a.data[i];          return r; }
vector_template vec_t  operator*(vec_t v, float s)   { vec_t r; for (int i=0; i<n; ++i) r.data[i] = v.data[i]*s;         return r; }
vector_template vec_t  operator*(float s, vec_t v)   { vec_t r; for (int i=0; i<n; ++i) r.data[i] = v.data[i]*s;         return r; }
vector_template vec_t  operator*(vec_t a, vec_t b)   { vec_t r; for (int i=0; i<n; ++i) r.data[i] = a.data[i]*b.data[i]; return r; }
vector_template vec_t  operator/(vec_t a, vec_t b)   { vec_t r; for (int i=0; i<n; ++i) r.data[i] = a.data[i]/b.data[i]; return r; }
vector_template vec_t  operator/(vec_t v, float s)   { vec_t r; for (int i=0; i<n; ++i) r.data[i] = v.data[i]/s;         return r; }

vector_template vec_t& operator+=(vec_t& a, vec_t b)  { for (int i = 0; i < n; ++i) a.data[i] += b.data[i]; return a; }
vector_template vec_t& operator-=(vec_t& a, vec_t b)  { for (int i = 0; i < n; ++i) a.data[i] -= b.data[i]; return a; }
vector_template vec_t& operator*=(vec_t& v, float s)  { for (int i = 0; i < n; ++i) v.data[i] *= s;         return v; }
vector_template vec_t& operator*=(vec_t& a, vec_t b)  { for (int i = 0; i < n; ++i) a.data[i] *= b.data[i]; return a; }
vector_template vec_t& operator/=(vec_t& a, vec_t b)  { for (int i = 0; i < n; ++i) a.data[i] /= b.data[i]; return a; }
vector_template vec_t& operator/=(vec_t& v, float s)  { for (int i = 0; i < n; ++i) v.data[i] /= s;         return v; }

vector_template
inline T dot(vec_t a, vec_t b) {
    T result = 0;
    for (int i = 0; i < n; ++i)
        result += a.data[i] * b.data[i];
    return result;
}
#undef vector_template
#undef vec_t

template<typename T>
inline Vector<3,T> m_cross(const Vector<3,T>& a, const Vector<3,T>& b) {
    return Vector<3,T>{a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

// Two-dimensional "cross-product" is the Z component of the cross-product between two 3D vectors where Z=0.
// Alternatively it can be thought as the dotproduct of a vector and another vector but rotated 90 degrees.
template<typename T>
inline T m_cross(const Vector<2,T>& a, const Vector<2,T>& b) {
    return a.x*b.y-a.y*b.x;
}

///////////// Vector matrix functions /////////////
template <int r, int c, typename T>
Vector<r, T> operator*(const Matrix<r, c, T>& m, const Vector<c, T>& x) {
    Vector<r, T> result = {};
    for (int col = 0; col < c; col++)
        for (int row = 0; row < r; row++)
            result.data[row] += m.data[row + col * r] * x.data[col];
    return result;
}

template<typename T>
Vector<2,T> operator*(const Matrix<2,2,T>& m, const Vector<2,T>& b) { return m.a1*b.x + m.a2*b.y; }
template<typename T>
Vector<3,T> operator*(const Matrix<3,3,T>& m, const Vector<3,T>& b) { return m.a1*b.x + m.a2*b.y + m.a3*b.z; }
template<typename T>
Vector<4,T> operator*(const Matrix<4,4,T>& m, const Vector<4,T>& b) { return m.a1*b.x + m.a2*b.y + m.a3*b.z + m.a4*b.w; }

template <int n, typename T>
inline float m_length(const Vector<n, T>& v) {
    float result = sqrt(m_dot(v, v));
    return result;
}

template <int n, typename T>
inline Vector<n, T> m_normalize(const Vector<n,T>& v) {
    Vector<n, T> result = v / m_length(v);
    return result;
}

///////////// Util functions for vector and matrix /////////////
template<int n, typename T>
std::ostream& operator<<(std::ostream& os, const Vector<n, T>& vec) {
    os << "[";
    for(int i = 0; i < n; ++i) {
        os << vec.data[i];
        if(i != n-1) os << ", ";
    }
    os << "]";
    return os;
}

template<int r, int c, typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<r,c, T>& mat) {
    for(int i = 0; i < r; ++i) {
        for(int j = 0; j < c; ++j) {
            int idx = j + i*r;
            if(j == 0)  os << "[" << mat[idx] << ", ";
            else if(j == r-1) os << mat[idx] << "]\n";
            else os << mat[idx] << ", ";
        }
    }
    return os;
}

#endif // MATRIX_H

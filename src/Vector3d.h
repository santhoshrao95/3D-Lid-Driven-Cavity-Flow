#pragma once

#include <cmath>

struct Vector3d
{
    double x, y, z;

    Vector3d() {}
    Vector3d(const double x, const double y, const double z) : x(x), y(y), z(z) {}
};

Vector3d negate(const Vector3d& v);
Vector3d operator - (const Vector3d& v);

Vector3d add(const Vector3d& a, const Vector3d& b);
Vector3d operator + (const Vector3d& a, const Vector3d& b);

Vector3d sub(const Vector3d& a, const Vector3d& b);
Vector3d operator - (const Vector3d& a, const Vector3d& b);

double dot(const Vector3d& a, const Vector3d& b);
double norm(const Vector3d& v);

Vector3d div(const Vector3d& v, const double d);
Vector3d operator / (const Vector3d& v, const double d);

Vector3d mul(const double s, const Vector3d& v);
Vector3d operator * (const double s, const Vector3d& v);

Vector3d normalized(const Vector3d& v);
Vector3d withMaxNorm(const Vector3d& v, const double maxNorm);
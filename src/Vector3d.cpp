#include "Vector3d.h"

Vector3d negate(const Vector3d& v)
{
    return Vector3d(-v.x, -v.y, -v.z);
}
Vector3d operator - (const Vector3d& v)
{
    return negate(v);
}

Vector3d add(const Vector3d& a, const Vector3d& b)
{
    return Vector3d(a.x + b.x, a.y + b.y, a.z + b.z);
}
Vector3d operator + (const Vector3d& a, const Vector3d& b)
{
    return add(a, b);
}

Vector3d sub(const Vector3d& a, const Vector3d& b)
{
    return Vector3d(a.x - b.x, a.y - b.y, a.z - b.z);
}
Vector3d operator - (const Vector3d& a, const Vector3d& b)
{
    return sub(a, b);
}

double dot(const Vector3d& a, const Vector3d& b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}
double norm(const Vector3d& v)
{
    return sqrt(dot(v, v));
}
Vector3d div(const Vector3d& v, const double d)
{
    return Vector3d(v.x / d, v.y / d, v.z / d);
}
Vector3d operator / (const Vector3d& v, const double d)
{
    return div(v, d);
}

Vector3d mul(const double s, const Vector3d& v)
{
    return Vector3d(s * v.x, s * v.y, s * v.z);
}
Vector3d operator * (const double s, const Vector3d& v)
{
    return mul(s, v);
}

Vector3d normalized(const Vector3d& v)
{
    const auto vNorm = norm(v);
    return (vNorm > 0) ? div(v, vNorm) : Vector3d(0, 0, 0);
}
Vector3d withMaxNorm(const Vector3d& v, const double maxNorm)
{
    const auto vNorm = norm(v);
    return (vNorm <= maxNorm) ? v : maxNorm * normalized(v);
}
#include "FiniteDifference.h"

double gradient1stOrderBackwardDiff(const std::vector<double>& f, const size_t i, const double dx)
{
    return (f[i] - f[i - 1]) / dx;
}
double gradient2ndOrderCentralDiff(const std::vector<double>& f, const size_t i, const double dx)
{
    return (f[i + 1] - (2 * f[i]) + f[i - 1]) / (dx * dx);
}
Vector3d gradient1stOrderCentralDiff(//modified
    const std::vector<std::vector<std::vector<double>>>& f,
    const size_t xIndex, const size_t yIndex, const size_t zIndex,
    const double dx, const double dy, const double dz
)
{
    const auto dfdx = (f[xIndex + 1][yIndex][zIndex] - f[xIndex - 1][yIndex][zIndex]) / (2 * dx);
    const auto dfdy = (f[xIndex][yIndex + 1][zIndex] - f[xIndex][yIndex - 1][zIndex]) / (2 * dy);
    const auto dfdz = (f[xIndex][yIndex][zIndex + 1] - f[xIndex][yIndex][zIndex - 1]) / (2 * dz);

    return Vector3d(dfdx, dfdy, dfdz);
}
double divergence1stOrderBackwardDiff( //modified
    const std::vector<std::vector<std::vector<Vector3d>>>& v,
    const size_t xIndex, const size_t yIndex, const size_t zIndex,
    const double dx, const double dy, const double dz
)
{
    const auto dvxdx = (v[xIndex][yIndex][zIndex].x - v[xIndex - 1][yIndex][zIndex].x) / dx;
    const auto dvydy = (v[xIndex][yIndex][zIndex].y - v[xIndex][yIndex - 1][zIndex].y) / dy;
    const auto dvzdz = (v[xIndex][yIndex][zIndex].z - v[xIndex][yIndex][zIndex - 1].z) / dz;

    return dvxdx + dvydy + dvzdz;
}
double laplacian2ndOrderCentralDiff( //modified
    const std::vector<std::vector<std::vector<double>>>& f,
    const size_t xIndex, const size_t yIndex, const size_t zIndex,
    const double dx, const double dy, const double dz
)
{
    const auto d2fdx2 = (f[xIndex + 1][yIndex][zIndex] - (2 * f[xIndex][yIndex][zIndex]) + f[xIndex - 1][yIndex][zIndex])
        / (dx * dx);
    const auto d2fdy2 = (f[xIndex][yIndex + 1][zIndex] - (2 * f[xIndex][yIndex][zIndex]) + f[xIndex][yIndex - 1][zIndex])
        / (dy * dy);
    const auto d2fdz2 = (f[xIndex][yIndex][zIndex + 1] - (2 * f[xIndex][yIndex][zIndex]) + f[xIndex][yIndex][zIndex - 1])
        / (dy * dy);


    return d2fdx2 + d2fdy2 + d2fdz2;
}
Vector3d laplacian3rdOrderCentralDiff( // modified
    const std::vector<std::vector<std::vector<Vector3d>>>& v,
    const size_t xIndex, const size_t yIndex, const size_t zIndex,
    const double dx, const double dy, const double dz
)
{
    const auto d2vxdx = (v[xIndex + 1][yIndex][zIndex].x - (2 * v[xIndex][yIndex][zIndex].x) + v[xIndex - 1][yIndex][zIndex].x)
        / (dx * dx);
    const auto d2vxdy = (v[xIndex][yIndex + 1][zIndex].x - (2 * v[xIndex][yIndex][zIndex].x) + v[xIndex][yIndex - 1][zIndex].x)
        / (dy * dy);
    const auto d2vxdz = (v[xIndex][yIndex][zIndex + 1].x - (2 * v[xIndex][yIndex][zIndex].x) + v[xIndex][yIndex][zIndex - 1].x)
        / (dz * dz);

    const auto d2vydx = (v[xIndex + 1][yIndex][zIndex].y - (2 * v[xIndex][yIndex][zIndex].y) + v[xIndex - 1][yIndex][zIndex].y)
        / (dx * dx);
    const auto d2vydy = (v[xIndex][yIndex + 1][zIndex].y - (2 * v[xIndex][yIndex][zIndex].y) + v[xIndex][yIndex - 1][zIndex].y)
        / (dy * dy);
    const auto d2vydz = (v[xIndex][yIndex][zIndex + 1].y - (2 * v[xIndex][yIndex][zIndex].y) + v[xIndex][yIndex][zIndex - 1].y)
        / (dz * dz);

    const auto d2vzdx = (v[xIndex + 1][yIndex][zIndex].z - (2 * v[xIndex][yIndex][zIndex].z) + v[xIndex - 1][yIndex][zIndex].z)
        / (dx * dx);
    const auto d2vzdy = (v[xIndex][yIndex + 1][zIndex].z - (2 * v[xIndex][yIndex][zIndex].z) + v[xIndex][yIndex - 1][zIndex].z)
        / (dy * dy);
    const auto d2vzdz = (v[xIndex][yIndex][zIndex + 1].z - (2 * v[xIndex][yIndex][zIndex].z) + v[xIndex][yIndex][zIndex - 1].z)
        / (dz * dz);
    return Vector3d(d2vxdx + d2vxdy + d2vxdz, d2vydx + d2vydy + d2vydz, d2vzdx + d2vzdy + d2vzdz);
}

/*std::vector<std::vector<double>> iteratePoissonsEquation(
    const std::vector<std::vector<double>>& p, const std::vector<std::vector<double>>& b,
    const size_t numX, const size_t numY, const double dx, const double dy)
{
    const auto dx2 = dx * dx;
    const auto dy2 = dy * dy;

    auto newP = p;
    for(size_t i = 1; i < (numX - 1); i++)
    {
        for(size_t j = 1; j < (numY - 1); j++)
        {
            const auto term1Numerator =
                (dy2 * (p[i + 1][j] + p[i - 1][j])) +
                (dx2 * (p[i][j + 1] + p[i][j - 1]));
            const auto termDenominator = 2 * (dx2 + dy2);
            const auto term1 = term1Numerator / termDenominator;
            const auto term2 = ((dx2 * dy2) / termDenominator) * b[i][j];

            newP[i][j] = term1 - term2;
        }
    }

    return newP;
}

std::vector<std::vector<double>> getBForIncompressibleNavierStokes(
    const std::vector<std::vector<Vector3d>>& v, const double rho,
    const size_t numX, const size_t numY, const double dx, const double dy, const double dt)
{
    std::vector<std::vector<double>> b(numX, std::vector<double>(numY));

    for(size_t i = 1; i < numX - 1; i++)
    {
        for(size_t j = 1; j < numY - 1; j++)
        {
            const auto term1 = (1.0 / dt) * divergence1stOrderBackwardDiff(v, i, j, dx, dy);
            const auto term2 = pow((v[i + 1][j].x - v[i - 1][j].x) / (2 * dx), 2);
            const auto term3 = 2 *
                ((v[i][j + 1].x - v[i][j - 1].x) / (2 * dy)) *
                ((v[i + 1][j].y - v[i - 1][j].y) / (2 * dx));
            const auto term4 = pow((v[i][j + 1].y - v[i][j - 1].y) / (2 * dy), 2);

            b[i][j] = rho * (term1 - term2 - term3 - term4);
        }
    }

    return b;
}*/
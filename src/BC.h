#pragma once

#include<vector>

void applyFlowVelocityBC(std::vector<std::vector<std::vector<long double>>> &u,
                        std::vector<std::vector<std::vector<long double>>> &v,
                        std::vector<std::vector<std::vector<long double>>> &w,
                        std::vector<std::vector<std::vector<long double>>> &p,
                        size_t nx, size_t ny, size_t nz,
                        size_t nxp2, size_t nyp2, size_t nzp2);

void applyPressureBC(std::vector<std::vector<std::vector<long double>>> &phi,
                    size_t nx, size_t ny, size_t nz,
                    size_t nxp2, size_t nyp2, size_t nzp2);
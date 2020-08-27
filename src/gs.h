#pragma once
//#ifndef NSSolverRK3_CN_H
//#define NSSolverRK3_CN_H
#include "BC.h"
#include <vector>
#include<math.h>
#include<iostream>


void gs(std::vector<std::vector<std::vector<long double>>> &phi,
        std::vector<std::vector<std::vector<long double>>> &rhsn,
        size_t nx, size_t ny, size_t nz,
        size_t nxp2, size_t nyp2, size_t nzp2,
        long double dx, long double dy, long double dz);
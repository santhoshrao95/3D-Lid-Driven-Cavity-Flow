#pragma once
//#ifndef NSSolverRK3_CN_H
//#define NSSolverRK3_CN_H
#include "BC.h"
#include "gs.h"
#include<fstream>
#include <vector>

//#endif
void NSSolverRK3_CN(std::vector<std::vector<std::vector<long double>>> &u,
                    std::vector<std::vector<std::vector<long double>>> &v,
                    std::vector<std::vector<std::vector<long double>>> &w,
                    std::vector<std::vector<std::vector<long double>>> &p,
                    size_t nx, size_t ny, size_t nz,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    long double dx, long double dy, long double dz, long double Time,const long double initialTimestep, long double Re);

long double calculateDelt(std::vector<std::vector<std::vector<long double>>> &u,
                    std::vector<std::vector<std::vector<long double>>> &v,
                    std::vector<std::vector<std::vector<long double>>> &w,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    const long double dx, const long double dy, const long double dz);

void rksubstep(long double a1, long double a2, long double a3, long double dt_rk, long double delT,
                    std::vector<std::vector<std::vector<long double>>> &u,
                    std::vector<std::vector<std::vector<long double>>> &v,
                    std::vector<std::vector<std::vector<long double>>> &w,
                    std::vector<std::vector<std::vector<long double>>> &p,
                    size_t nx, size_t ny, size_t nz,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    long double dx, long double dy, long double dz, long double Re);

void uRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &uTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re);

void vRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &vTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re);

void wRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &wTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re);

void var_advance(std::vector<std::vector<std::vector<long double>>> &v1,
                    std::vector<std::vector<std::vector<long double>>> &v2,
                    std::vector<std::vector<std::vector<long double>>> &v3,
                    const long double f1, const long double f2, const long double f3,
                     const long double nxp2, const long double nyp2, const long double nzp2);
void psource(std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        std::vector<std::vector<std::vector<long double>>> &v1Temp,
                        std::vector<std::vector<std::vector<long double>>> &v2Temp,
                        std::vector<std::vector<std::vector<long double>>> &v3Temp,
                        std::vector<std::vector<std::vector<long double>>> &src,
                        size_t nx, size_t ny, size_t nz,
                        size_t nxp2, size_t nyp2, size_t nzp2,
                        const long double dx, const long double dy, const long double dz, const long double rkdt);
void pgrad(std::vector<std::vector<std::vector<long double>>> &pn,
                        std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        const long double dx, const long double dy, const long double dz,
                        size_t nxp2, size_t nyp2, size_t nzp2, const long double dtrk);

void print(const long double lx, const long double ly, const long double lz,
            const long double dx, const long double dy, const long double dz,
            size_t nx, size_t ny, size_t nz,
            size_t nxp2, size_t nyp2, size_t nzp2, 
            std::vector<std::vector<std::vector<long double>>> &u,
            std::vector<std::vector<std::vector<long double>>> &v,
            std::vector<std::vector<std::vector<long double>>> &w,
            std::vector<std::vector<std::vector<long double>>> &p,
            size_t count, const long Time, const long delT);

void printlevel(size_t count, const long double totaltime, const long double delT);


 
#pragma once

//#include"BC.h"
#include"NSSolverRK3_CN.h"
#include<vector>

class CavityFlow{
    public:
    CavityFlow();
    CavityFlow(const long double xL,const long double yL,const long double zL,
                                         const long double timeLength, const long double initialTimestep, const long double Re,
                                         size_t nx, size_t ny, size_t nz);
    //void applyFlowVelocityBoundarycondition();
    void run();



    private:
    const long double Lx, Ly, Lz, timeLength,initialTimestep, Re;
    size_t nx, ny, nz;
    size_t nxp2, nyp2, nzp2;
    const long double dx = Lx / (nxp2-1), dy = Ly / (nyp2-1), dz = Lz / (nzp2-1);
    std::vector<std::vector<std::vector<long double>>> u;
    std::vector<std::vector<std::vector<long double>>> v;
    std::vector<std::vector<std::vector<long double>>> w;
    std::vector<std::vector<std::vector<long double>>> p;

};
#include"CavityFlow.h"
#include<vector>



//default constructor
CavityFlow::CavityFlow():Lx(1.0),Ly(1.0),Lz(1.0),timeLength(10.0), initialTimestep(0.001),
                                           Re(200.0),nx(32),ny(32),nz(32),
                                           nxp2(nx+2),nyp2(ny+2),nzp2(nz+2)
{
    u = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    v = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    w = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    p = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

}

//parameterized constructor
CavityFlow::CavityFlow(const long double xL,const long double yL,const long double zL,
                                         const long double timeLength, const long double initialTimestep, const long double Re,
                                         size_t nx, size_t ny, size_t nz):
                                         Lx(xL),Ly(yL),Lz(zL),timeLength(timeLength),initialTimestep(initialTimestep),Re(Re),
                                         nx(nx),ny(ny),nz(nz),nxp2(nx+2),nyp2(ny+2),nzp2(nz+2)
{
    u = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    v = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    w = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    p = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

}

/*void CavityFlow::applyFlowVelocityBoundarycondition()
{
    //x-min boundary condition

    for (size_t j = 0; j < nyp2; j++)
    {
        for (size_t k = 0; k < nzp2; k++)
        {
            u[0][j][k] = 0.0;
            v[0][j][k] = - v[1][j][k];
            w[0][j][k] = - w[1][j][k];
            p[0][j][k] = p[1][j][k];
        }    
    }
    //x-max boundary condition

    for (size_t j = 0; j < nyp2; j++)
    {
        for (size_t k = 0; k < nzp2; k++)
        {
            u[nx][j][k] = 0.0;
            v[nx][j][k] = - v[nx-1][j][k];
            w[nx][j][k] = - w[nx-1][j][k];
            p[nx+1][j][k] = p[nx][j][k];
        }    
    }
    //z-min boundary condition

    for (size_t i = 0; i < nxp2; i++)
    {
        for (size_t j = 0; j < nzp2; j++)
        {
            w[i][j][0] = 0.0;
            u[i][j][0] = - u[i][j][1];
            v[i][j][0] = - v[i][j][1];
            p[i][j][0] = p[i][j][1];
        }    
    }    
    //z-max boundary condition

    for (size_t i = 0; i < nxp2; i++)
    {
        for (size_t j = 0; j < nzp2; j++)
        {
            w[i][j][nz] = 0.0;
            u[i][j][nz] = 2.0 - u[i][j][nz-1];
            v[i][j][nz] = - v[i][j][nz-1];
            p[i][j][nz+1] = p[i][j][nz];
        }    
    }
    //y boundary condition //periodic

    for (size_t i = 0; i < nxp2; i++)
    {
        for (size_t k = 0; k < nzp2; k++)
        {
            u[i][0][k] = u[i][ny][k];
            u[i][ny+1][k] = u[i][1][k];

            v[i][0][k] = v[i][ny-1][k];
            v[i][ny][k] = v[i][1][k];

            w[i][0][k] = w[i][ny][k];
            w[i][ny+1][k] = w[i][1][k];

            p[i][0][k] = p[i][ny][k];
            p[i][ny+1][k] = p[i][1][k];
        }    
    }
}*/

void CavityFlow::run()
{
    NSSolverRK3_CN(u,v,w,p,nx,ny,nz,nxp2,nyp2,nzp2,dx,dy,dz,timeLength,initialTimestep, Re);
}
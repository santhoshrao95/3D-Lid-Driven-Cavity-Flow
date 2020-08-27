#include"BC.h"

void applyFlowVelocityBC(std::vector<std::vector<std::vector<long double>>> &u,
                        std::vector<std::vector<std::vector<long double>>> &v,
                        std::vector<std::vector<std::vector<long double>>> &w,
                        std::vector<std::vector<std::vector<long double>>> &p,
                        size_t nx, size_t ny, size_t nz,
                        size_t nxp2, size_t nyp2, size_t nzp2)
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
            v[nx+1][j][k] = - v[nx][j][k];
            w[nx+1][j][k] = - w[nx][j][k];
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
            u[i][j][nz+1] = 2.0 - u[i][j][nz];
            v[i][j][nz+1] = - v[i][j][nz];
            p[i][j][nz+1] = - p[i][j][nz];
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
}


void applyPressureBC(std::vector<std::vector<std::vector<long double>>> &phi,
                    size_t nx, size_t ny, size_t nz,
                    size_t nxp2, size_t nyp2, size_t nzp2)
{
    
    for (size_t k = 0; k < nzp2; k++)
    {
        for (size_t j = 0; j < nyp2; j++)
        {
            for (size_t i = 0; i < nxp2; i++)
            {
                phi[0][j][k] =  phi[1][j][k];
                phi[nx+1][j][k] = phi[nx][j][k];
                phi[i][j][0] =   phi[i][j][1];
                phi[i][j][nz+1] =  - phi[i][j][nz];
                //periodic BC
                phi[i][0][k] = phi[i][ny][k];
                phi[i][ny+1][k] = phi[i][1][k];

            }
            
        }
        
    }
    

}
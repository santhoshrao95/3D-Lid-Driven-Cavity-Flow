#include"gs.h"


void gs(std::vector<std::vector<std::vector<long double>>> &phi,
        std::vector<std::vector<std::vector<long double>>> &rhsn,
        size_t nx, size_t ny, size_t nz,
        size_t nxp2, size_t nyp2, size_t nzp2,
        long double dx, long double dy, long double dz)
        {
        long double ap = (-2.0/(dx*dx))+(-2.0/(dy*dy))+(-2.0/(dz*dz));
        long double ax = (1.0/(dx*dx)), ay = (1.0/(dy*dy)), az = (1.0/(dz*dz));
        std::vector<std::vector<std::vector<long double>>> phiold;
        phiold = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
        long double l2norm = 0.0;
        phiold = phi;
        do
        {
        for (size_t k = 1; k <(nzp2-1); k++)
        {
            for (size_t j = 1; j < (nyp2-1); j++)
            {
                for (size_t i = 1; i < (nxp2-1); i++)
                {
                    phi[i][j][k] = (rhsn[i][j][k] - ((ax*(phi[i][j][k+1]+phi[i][j][k-1])) 
                                               + (ay*(phi[i+1][j][k]+phi[i-1][j][k])) 
                                                + (az*(phi[i][j+1][k]+phi[i][j-1][k]))))/ap;
                }
            }
        }        

        double err = 0.0;
        for (int k = 1; k < (nzp2-1); k++) // calculate err
        {
            for (int j = 1; j < (nyp2-1); j++)
            {
                for (int i = 1; i < (nxp2-1); i++)
                {
                    err = err+ pow((phi[i][j][k]-phiold[i][j][k]),2.0);
                }

            }

        }    

        l2norm = sqrt(err);

        std::cout << "l2norm is " << l2norm << std::endl;

        phiold = phi;
        applyPressureBC(phi,nx,ny,nz,nxp2,nyp2,nzp2);

        } while (l2norm>0.0001);
        

        applyPressureBC(phi,nx,ny,nz,nxp2,nyp2,nzp2);
        }


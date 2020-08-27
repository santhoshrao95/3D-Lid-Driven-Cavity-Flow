//#pragma once

#include <vector>
#include <algorithm>
#include "NSSolverRK3_CN.h"

void printlevel(size_t count, const long double totaltime, const long double delT)
{
    std::ofstream output("whereiam", std::ios::app);
    output<<"I have entered count: " << count << " " << totaltime << " " << delT << std::endl;
    output.close();
}
void print(const long double lx, const long double ly, const long double lz,
            const long double dx, const long double dy, const long double dz,
            size_t nx, size_t ny, size_t nz,
            size_t nxp2, size_t nyp2, size_t nzp2, 
            std::vector<std::vector<std::vector<long double>>> &u,
            std::vector<std::vector<std::vector<long double>>> &v,
            std::vector<std::vector<std::vector<long double>>> &w,
            std::vector<std::vector<std::vector<long double>>> &p,
            size_t count, const long double Time, const long double delT)
{
  const double x0 = 0.0, y0 = 0.0, z0 = 0.0;
  std::vector<double> xe(nxp2, 0.0),ye(nxp2, 0.0),ze(nxp2, 0.0),xc(nxp2, 0.0),yc(nxp2, 0.0),zc(nxp2, 0.0);
    for (size_t i = 0; i < nxp2; i++)
    {
        xe[i] = x0 + lx/double(nxp2-1) * double(i+1) - 0.5 * lx / double(nxp2-1);
        ye[i] = x0 + ly/double(nyp2-1) * double(i+1) - 0.5 * ly / double(nyp2-1);
        ze[i] = z0 + lz/double(nzp2-1) * double(i+1) - 0.5 * lz / double(nzp2-1);
    }
    
    for (size_t i = 1; i < nxp2; i++)
    {
        xc[i] = (xe[i]+xe[i-1]) / 2.0;
        yc[i] = (ye[i]+ye[i-1]) / 2.0;
        zc[i] = (ze[i]+ze[i-1]) / 2.0;
    }
    
    xc[0] = xe[0] - 0.5 * (xe[1]-xe[0]);
    yc[0] = ye[0] - 0.5 * (ye[1]-ye[0]);
    zc[0] = ze[0] - 0.5 * (ze[1]-ze[0]);

    std::string init("output");
    std::string add(std::to_string(count));
    //init.append(add);
    //init.append("_");
    //init.append(std::to_string(Time));
    init.append(".csv.");
    init.append(add);


    std::ofstream output(init, std::ios::trunc);
    output << "x " <<"y " << "z " << "u " <<"v " << "w " << "p " <<std::endl; 
    for (size_t k = 0; k < nzp2; k++)
    {
        for (size_t j = 0; j < nyp2; j++)
        {
            for (size_t i = 0; i < nxp2; i++)
            {
                output << xc[i] << " " << yc[j] << " " << zc[k] << " "  << u[i][j][k] << " " << v[i][j][k] << " " << w[i][j][k] << " " << p[i][j][k] << " " << std::endl;  
            }
            
        }
        
    }
    output.close();

}

long double calculateDelt(std::vector<std::vector<std::vector<long double>>> &u,
                    std::vector<std::vector<std::vector<long double>>> &v,
                    std::vector<std::vector<std::vector<long double>>> &w,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    const long double dx, const long double dy, const long double dz)
{
    const long double cmax = 1.2; 
    long double umax = 0.0, vmax = 0.0, wmax = 0.0;
  for (size_t k = 1; k < nzp2-1; k++)
  {
    for (size_t j = 1; j < nyp2-1; j++)
    {
      for (size_t i = 1; i < nxp2-2; i++)
      {
              if (umax<abs(u[i][j][k]))
                umax = abs(u[i][j][k]);
              /*if (vmax<v[i][j][k])
                vmax = v[i][j][k];
              if (wmax<w[i][j][k])
                wmax = w[i][j][k];*/
            }            
        }       
    }
  for (size_t k = 1; k < nzp2-1; k++)
  {
    for (size_t j = 1; j < nyp2-2; j++)
    {
      for (size_t i = 1; i < nxp2-1; i++)
      {
              /*if (umax<u[i][j][k])
                umax = u[i][j][k];*/
              if (vmax<abs(v[i][j][k]))
                vmax = abs(v[i][j][k]);
              /*if (wmax<w[i][j][k])
                wmax = w[i][j][k];*/
            }            
        }       
    }
  for (size_t k = 1; k < nzp2-2; k++)
  {
    for (size_t j = 1; j < nyp2-1; j++)
    {
      for (size_t i = 1; i < nxp2-1; i++)
      {
              /*if (umax<u[i][j][k])
                umax = u[i][j][k];
              if (vmax<v[i][j][k])
                vmax = v[i][j][k];*/
              if (wmax<abs(w[i][j][k]))
                wmax = abs(w[i][j][k]);
            }            
        }       
    }
    long double deltx = cmax * dx / umax;
    long double delty = cmax * dy / vmax;
    long double deltz = cmax * dz / wmax;

    return std::min(std::min(deltx, delty), deltz);
}

void uRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &uTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re)
{
  long double uaver1 = 0.0, uaver2 = 0.0, vaver1 = 0.0, vaver2 = 0.0, waver1 = 0.0, waver2 = 0.0;
  long double rdx = 1.0/dx , rdy = 1.0/dy, rdz = 1.0/dz, rRe = 1.0/Re;
  //uTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  
  for (size_t k = 1; k < nzp2-1; k++)
  {
    for (size_t j = 1; j < nyp2-1; j++)
    {
      for (size_t i = 1; i < nxp2-2; i++)
      {
        uaver2 = 0.5 * ( un[i+1][j][k] + un[i][j][k] );
        uaver1 = 0.5 * ( un[i][j][k] + un[i-1][j][k] );
        uTemp[i][j][k] = uTemp[i][j][k] + ( uaver1*uaver1 - uaver2*uaver2) * rdx;


        uaver2 = 0.5 * ( un[i][j][k+1] + un[i][j][k] );
        uaver1 = 0.5 * ( un[i][j][k] + un[i][j][k-1] );
        waver2 = 0.5 * ( wn[i+1][j][k] + wn[i][j][k] );
        waver1 = 0.5 * ( wn[i+1][j][k-1] + wn[i][j][k-1] );
        uTemp[i][j][k] = uTemp[i][j][k] + ( uaver1*waver1 - uaver2*waver2) * rdz;

        uaver2 = 0.5 * ( un[i][j+1][k] + un[i][j][k] );
        uaver1 = 0.5 * ( un[i][j][k] + un[i][j-1][k] );
        waver2 = 0.5 * ( vn[i+1][j][k] + vn[i][j][k] );
        waver1 = 0.5 * ( vn[i+1][j-1][k] + vn[i][j-1][k] );
        uTemp[i][j][k] = uTemp[i][j][k] + ( uaver1*vaver1 - uaver2*vaver2) * rdy;

        uTemp[i][j][k] = uTemp[i][j][k] + rRe*( ( un[i+1][j][k] - un[i][j][k] )*rdx
                                              - ( un[i][j][k] - un[i-1][j][k] )*rdx )*rdx               
                                        + rRe*( ( un[i][j+1][k] - un[i][j][k] )*rdy  
                                              - ( un[i][j][k] - un[i][j-1][k] )*rdy )*rdy               
                                        + rRe*( ( un[i][j][k+1] - un[i][j][k] )*rdy   
                                              - ( un[i][j][k] - un[i][j][k-1] )*rdy )*rdy ; 
      }
      
    }
    
  }

}

void vRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &vTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re)
{
  long double uaver1 = 0.0, uaver2 = 0.0, vaver1 = 0.0, vaver2 = 0.0, waver1 = 0.0, waver2 = 0.0;
  long double rdx = 1.0/dx , rdy = 1.0/dy, rdz = 1.0/dz, rRe = 1.0/Re;
  // vTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  
  for (size_t k = 1; k < nzp2-1; k++)
  {
    for (size_t j = 1; j < nyp2-2; j++)
    {
      for (size_t i = 1; i < nxp2-1; i++)
      {
        vaver2 = 0.5 * ( vn[i][j+1][k] + vn[i][j][k] );
        vaver1 = 0.5 * ( vn[i][j][k] + vn[i][j-1][k] );
        vTemp[i][j][k] = vTemp[i][j][k] + ( vaver1*vaver1 - vaver2*vaver2) * rdy;

        vaver2 = 0.5 * ( vn[i+1][j][k] + vn[i][j][k] );
        vaver1 = 0.5 * ( vn[i][j][k] + vn[i-1][j][k] );
        uaver2 = 0.5 * ( un[i][j+1][k] + un[i][j][k] );
        uaver1 = 0.5 * ( un[i-1][j+1][k] + un[i-1][j][k] );
        vTemp[i][j][k] = vTemp[i][j][k] + ( uaver1*vaver1 - uaver2*vaver2) * rdx;

        vaver2 = 0.5 * ( vn[i][j][k+1] + vn[i][j][k] );
        vaver1 = 0.5 * ( vn[i][j][k] + vn[i][j][k-1] );
        waver2 = 0.5 * ( wn[i][j+1][k] + wn[i][j][k] );
        waver1 = 0.5 * ( wn[i][j+1][k-1] + wn[i][j][k-1] );
        vTemp[i][j][k] = vTemp[i][j][k] + ( vaver1*waver1 - vaver2*waver2) * rdz;

        vTemp[i][j][k] = vTemp[i][j][k] + rRe*( ( vn[i+1][j][k] - vn[i][j][k] )*rdx
                                              - ( vn[i][j][k] - vn[i-1][j][k] )*rdx )*rdx 
                                        + rRe*( ( vn[i][j+1][k] - vn[i][j][k] )*rdy 
                                              - ( vn[i][j][k] - vn[i][j-1][k] )*rdy )*rdy 
                                        + rRe*( ( vn[i][j][k+1] - vn[i][j][k] )*rdy 
                                              - ( vn[i][j][k] - vn[i][j][k-1] )*rdy )*rdy ; 
      }
      
    }
    
  }

}
void wRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &wTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re)
{
  long double uaver1 = 0.0, uaver2 = 0.0, vaver1 = 0.0, vaver2 = 0.0, waver1 = 0.0, waver2 = 0.0;
  long double rdx = 1.0/dx , rdy = 1.0/dy, rdz = 1.0/dz, rRe = 1.0/Re;
  //wTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  
  for (size_t k = 1; k < nzp2-2; k++)
  {
    for (size_t j = 1; j < nyp2-1; j++)
    {
      for (size_t i = 1; i < nxp2-1; i++)
      {
        waver2 = 0.5 * ( wn[i][j][k+1] + wn[i][j][k] );
        waver1 = 0.5 * ( wn[i][j][k] + wn[i][j][k-1] );
        wTemp[i][j][k] = wTemp[i][j][k] + ( waver1*waver1 - waver2*waver2) * rdz;

        waver2 = 0.5 * ( wn[i+1][j][k] + wn[i][j][k] );
        waver1 = 0.5 * ( wn[i][j][k] + wn[i-1][j][k] );
        uaver2 = 0.5 * ( un[i][j][k+1] + un[i][j][k] );
        uaver1 = 0.5 * ( un[i-1][j][k+1] + un[i-1][j][k] );
        wTemp[i][j][k] = wTemp[i][j][k] + ( uaver1*waver1 - uaver2*waver2) * rdx;
        
        waver2 = 0.5 * ( wn[i][j+1][k] + wn[i][j][k] );
        waver1 = 0.5 * ( wn[i][j][k] + wn[i][j-1][k] );
        vaver2 = 0.5 * ( vn[i][j][k+1] + vn[i][j][k] );
        vaver1 = 0.5 * ( vn[i][j-1][k+1] + vn[i][j-1][k] );
        wTemp[i][j][k] = wTemp[i][j][k] + ( vaver1*waver1 - vaver2*waver2) * rdy;

        wTemp[i][j][k] = wTemp[i][j][k] + rRe*( ( wn[i+1][j][k] - wn[i][j][k] )*rdx
                                              - ( wn[i][j][k] - wn[i-1][j][k] )*rdx )*rdx               
                                        + rRe*( ( wn[i][j+1][k] - wn[i][j][k] )*rdy
                                              - ( wn[i][j][k] - wn[i][j-1][k] )*rdy )*rdy
                                        + rRe*( ( wn[i][j][k+1] - wn[i][j][k] )*rdy 
                                              - ( wn[i][j][k] - wn[i][j][k-1] )*rdy )*rdy ; 
      }
      
    }
    
  }

}


void var_advance(std::vector<std::vector<std::vector<long double>>> &v1, 
                    std::vector<std::vector<std::vector<long double>>> &v2,
                    std::vector<std::vector<std::vector<long double>>> &v3,
                    const long double f1, const long double f2, const long double f3,
                    const long double nxp2, const long double nyp2, const long double nzp2)
                    //v1-->qU-->(EU*delT) v2-->UTemp1(c2*qU) v3-->FTemp(EU)
{
  for (size_t k = 0; k < nzp2; k++)
  {
    for (size_t j = 0; j < nyp2; j++)
    {
      for (size_t i = 0; i < nxp2; i++)
      {
        v1[i][j][k] = f1*v1[i][j][k] + v3[i][j][k]*f3;
        v2[i][j][k] = f2*v1[i][j][k];
      } 
    } 
  }  
}

void psource(std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        std::vector<std::vector<std::vector<long double>>> &v1Temp,
                        std::vector<std::vector<std::vector<long double>>> &v2Temp,
                        std::vector<std::vector<std::vector<long double>>> &v3Temp,
                        std::vector<std::vector<std::vector<long double>>> &src,
                        size_t nx, size_t ny, size_t nz,
                        size_t nxp2, size_t nyp2, size_t nzp2,
                        const long double dx, const long double dy, const long double dz, const long double rkdt)
{
  const long double rdx = 1.0/dx;
  const long double rdy = 1.0/dy;
  const long double rdz = 1.0/dz;
  const long double rrkdt = 1.0/rkdt;

  for (size_t i = 0; i < nxp2; i++)
  {
    for (size_t j = 0; j < nyp2; j++)
    {
      for (size_t k = 0; k < nzp2; k++)
      {
        v1Temp[i][j][k] = v1Temp[i][j][k] + un[i][j][k];
        v2Temp[i][j][k] = v2Temp[i][j][k] + vn[i][j][k];
        v3Temp[i][j][k] = v3Temp[i][j][k] + wn[i][j][k];
      }
      
    }
    
  }
  
  for (size_t k = 1; k < nzp2; k++)
  {
    for (size_t j = 1; j < nyp2; j++)
    {
      for (size_t i = 1; i < nzp2; i++)
      {
        src[i][j][k] = src[i][j][k] +( (v1Temp[i][j][k] - v1Temp[i-1][j][k]) * rdx \
                                                  + (v2Temp[i][j][k] - v2Temp[i][j-1][k]) * rdy \
                                                  + (v3Temp[i][j][k] - v3Temp[i][j][k-1]) * rdz ) * rrkdt;
      }
      
    }
    
  }
}
void pgrad(std::vector<std::vector<std::vector<long double>>> &pn,
                        std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        const long double dx, const long double dy, const long double dz,
                        size_t nxp2, size_t nyp2, size_t nzp2, const long double dtrk)
{
  const long double rdx = 1.0/dx, rdy = 1.0/dy, rdz = 1.0/dz;
  for (size_t k = 0; k < nzp2; k++)
  {
    for (size_t j = 0; j < nyp2; j++)
    {
      for (size_t i = 0; i < nxp2-1; i++)
      {
        un[i][j][k] = un[i][j][k] - ( pn[i+1][j][k]-pn[i][j][k] )*rdx * dtrk;
      } 
    } 
  }  
  for (size_t k = 0; k < nzp2; k++)
  {
    for (size_t j = 0; j < nyp2-1; j++)
    {
      for (size_t i = 0; i < nxp2; i++)
      {
        vn[i][j][k] = vn[i][j][k] - ( pn[i][j+1][k]-pn[i][j][k] )*rdy * dtrk;
      } 
    } 
  }
  for (size_t k = 0; k < nzp2-1; k++)
  {
    for (size_t j = 0; j < nyp2; j++)
    {
      for (size_t i = 0; i < nxp2-1; i++)
      {
        wn[i][j][k] = wn[i][j][k] - ( pn[i][j][k+1]-pn[i][j][k] )*rdz * dtrk;
      } 
    } 
  }

}

void rksubstep(long double a1, long double a2, long double a3, long double dt_rk, long double delT,
                    std::vector<std::vector<std::vector<long double>>> &m_un,
                    std::vector<std::vector<std::vector<long double>>> &m_vn,
                    std::vector<std::vector<std::vector<long double>>> &m_wn,
                    std::vector<std::vector<std::vector<long double>>> &m_pn,
                    std::vector<std::vector<std::vector<long double>>> &m_qu,
                    std::vector<std::vector<std::vector<long double>>> &m_qv,
                    std::vector<std::vector<std::vector<long double>>> &m_qw,
                    size_t nx, size_t ny, size_t nz,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    long double dx, long double dy, long double dz, long double Re)
{
  std::vector<std::vector<std::vector<long double>>> FTemp;
  std::vector<std::vector<std::vector<long double>>> rhs;
  std::vector<std::vector<std::vector<long double>>> uTemp1;
  std::vector<std::vector<std::vector<long double>>> vTemp1;
  std::vector<std::vector<std::vector<long double>>> wTemp1;
  
  FTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  rhs = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  uTemp1 = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  vTemp1 = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  wTemp1 = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

  uRhs(m_un,m_vn,m_wn,FTemp,nxp2,nyp2,nzp2,dx,dy,dz,Re);
  var_advance(m_qu, uTemp1, FTemp, a1, a2, delT, nxp2, nyp2, nzp2);

  FTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

  vRhs(m_un,m_vn,m_wn,FTemp,nxp2,nyp2,nzp2,dx,dy,dz,Re);
  var_advance(m_qv, vTemp1, FTemp, a1, a2, delT, nxp2, nyp2, nzp2);
  FTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

  wRhs(m_un,m_vn,m_wn,FTemp,nxp2,nyp2,nzp2,dx,dy,dz,Re);
  var_advance(m_qw, wTemp1, FTemp, a1, a2, delT, nxp2, nyp2, nzp2);
  FTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

  applyFlowVelocityBC(uTemp1, vTemp1, wTemp1,  m_pn, nx, ny, nz, nxp2, nyp2, nzp2);

  psource(m_un, m_vn, m_wn, uTemp1, vTemp1, wTemp1, rhs, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, dt_rk);

  gs(m_pn,rhs,nx,ny,nz,nxp2,nyp2,nzp2,dx,dy,dz);

  pgrad(m_pn,uTemp1, vTemp1, wTemp1,dx,dy,dz,nxp2,nyp2,nzp2,dt_rk);

  applyFlowVelocityBC(uTemp1, vTemp1, wTemp1,m_pn,nx,ny,nz,nxp2,nyp2,nzp2);

  m_un = uTemp1;
  m_vn = vTemp1;
  m_wn = wTemp1;

  std::vector<std::vector<std::vector<long double>>>().swap(FTemp);
  std::vector<std::vector<std::vector<long double>>>().swap(rhs);
  std::vector<std::vector<std::vector<long double>>>().swap(uTemp1);
  std::vector<std::vector<std::vector<long double>>>().swap(vTemp1);
  std::vector<std::vector<std::vector<long double>>>().swap(wTemp1);

}


void NSSolverRK3_CN(std::vector<std::vector<std::vector<long double>>> &un,
                    std::vector<std::vector<std::vector<long double>>> &vn,
                    std::vector<std::vector<std::vector<long double>>> &wn,
                    std::vector<std::vector<std::vector<long double>>> &pn,
                    size_t nx, size_t ny, size_t nz,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    long double dx, long double dy, long double dz, long double Time,const long double initialTimestep, long double Re)
{
  int count = 0;
  long double totaltime = 0.0;
  long double delT = 0.0;
  while(totaltime<Time)
  {
  applyFlowVelocityBC(un,vn,wn,pn,nx,ny,nz,nxp2,nyp2,nzp2);
  
  if(count == 0)
    delT = initialTimestep;
  else
    delT = calculateDelt(un,vn,wn,nxp2,nyp2,nzp2,dx,dy,dz);
  totaltime = totaltime + delT;
  printlevel(count,totaltime,delT);

  std::vector<std::vector<std::vector<long double>>> qu;
  std::vector<std::vector<std::vector<long double>>> qv;
  std::vector<std::vector<std::vector<long double>>> qw;

  qu = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  qv = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  qw = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

  long double a1 = 0.0, a2 = 0.0, a3 = 0.0, dt_rk = 0.0;

  a1 = 0.0;
  a2 = 1.0;
  a3 = 1.0 / 3.0;
  dt_rk = 1.0 / 3.0 * delT;
  rksubstep(a1, a2, a3, dt_rk, delT, un, vn, wn, pn, qu, qv, qw, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, Re);

  a1 = -5.0/9.0;
  a2 = 15.0/16.0;
  a3 = 5.0/12.0;
  dt_rk = 5.0/12.0 * delT;
  rksubstep(a1, a2, a3, dt_rk, delT, un, vn, wn, pn, qu, qv, qw, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, Re);

  a1 = -153.0/128.0;
  a2 = 8.0/15.0;
  a3 = 1.0/4.0;
  dt_rk = 1.0 / 4.0 * delT;
  rksubstep(a1, a2, a3, dt_rk, delT, un, vn, wn, pn, qu, qv, qw, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, Re);

  std::vector<std::vector<std::vector<long double>>>().swap(qu);
  std::vector<std::vector<std::vector<long double>>>().swap(qv);
  std::vector<std::vector<std::vector<long double>>>().swap(qw);

  if((count%5)==0)
  print(1.0,1.0,1.0,dx,dy,dz,nx,ny,nz,nxp2,nyp2,nzp2,un,vn,wn,pn,count,totaltime,delT);




  count++;
  }

  
    
}

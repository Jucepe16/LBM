/* This code simulates the poiseuille flow in 2D using BB boundary conditions
   at top and bottom
 */
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "LBM.h"

// Pressure conditions
double gradP=8.*nu*u_max/(NY*NY);
double rho_outlet=1.;
double rho_inlet=3.*(NX-1)*gradP+rho_outlet;

// Analytical solution: Poiseuille velocity
double ybottom=0.;
double ytop=NY;
double u_analy=-4.*u_max/(NY*NY)*(y-ybottom)*(y-ytop);

// convergence parameters
double tol=1e-12;      // tolerance to steady state convergence
double teval=100.;      // time step to evaluate convergence
double u_old[NX,NY]=0.;

void init_equilibrium(double *f0, double *f1, double *r, double *u, double *v)
{
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
        r[scalar_index(x,y)]=1.;
        u[scalar_index(x,y)]=0.;
        v[scalar_index(x,y)]=0.;

          double rho = r[scalar_index(x,y)];
            double ux  = u[scalar_index(x,y)];
            double uy  = v[scalar_index(x,y)];

                  // load equilibrium
            // feq_i  = w_i rho [1 + 3(ci . u) + (9/2) (ci . u)^2 - (3/2) (u.u)]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u) + (1/2) (ci . 3u)^2]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u){ 1 + (1/2) (ci . 3u) }]

            // temporary variables
            double w0r = w0;//*rho;
            double wsr = ws;//*rho;
            double wdr = wd;//*rho;
            double omusq = 1.0;// - 1.5*(ux*ux+uy*uy);

            double tux = 3.0*ux;
            double tuy = 3.0*uy;

            f0[field0_index(x,y)]    = w0r*(omusq);

            double cidot3u = tux;
            f1[fieldn_index(x,y,1)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy;
            f1[fieldn_index(x,y,2)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tux;
            f1[fieldn_index(x,y,3)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tuy;
            f1[fieldn_index(x,y,4)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));

            cidot3u = tux+tuy;
            f1[fieldn_index(x,y,5)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy-tux;
            f1[fieldn_index(x,y,6)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -(tux+tuy);
            f1[fieldn_index(x,y,7)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tux-tuy;
            f1[fieldn_index(x,y,8)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
        }
    }
}

void stream_collide_save(double *f0, double *f1, double *f2, double *r, double *u, double *v, bool save)
{
    // useful constants
    //const double tauinv = 2.0/(6.0*nu+1.0); // 1/tau
    const double omtauinv = 1.0-omega;//tauinv;     // 1 - 1/tau

    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            unsigned int xp1 = (x+1)%NX;
            unsigned int yp1 = (y+1)%NY;
            unsigned int xm1 = (NX+x-1)%NX;
            unsigned int ym1 = (NY+y-1)%NY;

            // direction numbering scheme
            // 6 2 5
            // 3 0 1
            // 7 4 8

            double ft0 = f0[field0_index(x,y)];

            // load populations from adjacent nodes
            double ft1 = f1[fieldn_index(xm1,y,  1)];
            double ft2 = f1[fieldn_index(x,  ym1,2)];
            double ft3 = f1[fieldn_index(xp1,y,  3)];
            double ft4 = f1[fieldn_index(x,  yp1,4)];
            double ft5 = f1[fieldn_index(xm1,ym1,5)];
            double ft6 = f1[fieldn_index(xp1,ym1,6)];
            double ft7 = f1[fieldn_index(xp1,yp1,7)];
            double ft8 = f1[fieldn_index(xm1,yp1,8)];

            // compute moments
            double rho = ft0+ft1+ft2+ft3+ft4+ft5+ft6+ft7+ft8;
            double rhoinv = 1.0/rho;

            double ux = rhoinv*(ft1+ft5+ft8-(ft3+ft6+ft7));
            double uy = rhoinv*(ft2+ft5+ft6-(ft4+ft7+ft8));

            // only write to memory when needed
            if(save)
            {
                r[scalar_index(x,y)] = rho;
                u[scalar_index(x,y)] = ux;
                v[scalar_index(x,y)] = uy;
            }

            // now compute and relax to equilibrium
            // note that
            // feq_i  = w_i rho [1 + 3(ci . u) + (9/2) (ci . u)^2 - (3/2) (u.u)]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u) + (1/2) (ci . 3u)^2]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u){ 1 + (1/2) (ci . 3u) }]

            // temporary variables
            double tw0r = omega*w0*rho; //   w[0]*rho/tau
            double twsr = omega*ws*rho; // w[1-4]*rho/tau
            double twdr = omega*wd*rho; // w[5-8]*rho/tau
            double omusq = 1.0 - 1.5*(ux*ux+uy*uy); // 1-(3/2)u.u

            double tux = 3.0*ux;
            double tuy = 3.0*uy;

            // Inlet/Outlet BC: PBBC (w/ i=1 and i=NX outside layers)
            for (int k=0,k<NPOP,k++){
                f(1,k,0) = w0*(rho_inlet+3*(cx(k)*u(NX-1,:)+cy(k)*v(NX-1,:)))+(f(NX-1,:,k)-feq(NX-1,:,k));
                f(NX,k,0)= w0*(rho_outlet+3*(cx(k)*u(2,:)+cy(k)*v(2,:)))+(f(2,:,k)-feq(2,:,k));
            }


            f0[field0_index(x,y)]    = omtauinv*ft0  + tw0r*(omusq);

            double cidot3u = tux;
            f2[fieldn_index(x,y,1)]  = omtauinv*ft1  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy;
            f2[fieldn_index(x,y,2)]  = omtauinv*ft2  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tux;
            f2[fieldn_index(x,y,3)]  = omtauinv*ft3  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tuy;
            f2[fieldn_index(x,y,4)]  = omtauinv*ft4  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));

            cidot3u = tux+tuy;
            f2[fieldn_index(x,y,5)]  = omtauinv*ft5  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy-tux;
            f2[fieldn_index(x,y,6)]  = omtauinv*ft6  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -(tux+tuy);
            f2[fieldn_index(x,y,7)]  = omtauinv*ft7  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tux-tuy;
            f2[fieldn_index(x,y,8)]  = omtauinv*ft8  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
        }
    }
}

void compute_flow_properties(unsigned int t, double *r, double *u, double *v, double *prop)
{
    // prop must point to space for 4 doubles:
    // 0: energy
    // 1: L2 error in rho
    // 2: L2 error in ux
    // 3: L2 error in uy

    double E = 0.0; // kinetic energy

    double sumrhoe2 = 0.0; // sum of error squared in rho
    double sumuxe2 = 0.0;  //                         ux
    double sumuye2 = 0.0;  //                         uy

    double sumrhoa2 = 0.0; // sum of analytical rho squared
    double sumuxa2 = 0.0;  //                   ux
    double sumuya2 = 0.0;  //                   uy

    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            double rho = r[scalar_index(x,y)];
            double ux  = u[scalar_index(x,y)];
            double uy  = v[scalar_index(x,y)];
            E += rho*(ux*ux + uy*uy);

            double rhoa, uxa, uya;
            taylor_green(t,x,y,&rhoa,&uxa,&uya);

            sumrhoe2 += (rho-rhoa)*(rho-rhoa);
            sumuxe2  += (ux-uxa)*(ux-uxa);
            sumuye2  += (uy-uya)*(uy-uya);

            sumrhoa2 += (rhoa-rho0)*(rhoa-rho0);
            sumuxa2  += uxa*uxa;
            sumuya2  += uya*uya;
        }
    }

    prop[0] = E;
    prop[1] = sqrt(sumrhoe2/sumrhoa2);
    prop[2] = sqrt(sumuxe2/sumuxa2);
    prop[3] = sqrt(sumuye2/sumuya2);
}

void report_flow_properties(unsigned int t, double *rho, double *ux, double *uy)
{
    double prop[4];
    compute_flow_properties(t,rho,ux,uy,prop);
    printf("%u,%g,%g,%g,%g\n",t,prop[0],prop[1],prop[2],prop[3]);
}

void save_scalar(const char* name, double *scalar, unsigned int n)
{
    // assume reasonably-sized file names
    char filename[128];
    char format[16];

    // compute maximum number of digits
    int ndigits = floor(log10((double)NSTEPS)+1.0);

    // generate format string
    // file name format is name0000nnn.bin
    sprintf(format,"%%s%%0%dd.bin",ndigits);
    sprintf(filename,format,name,n);

    // open file for writing
    FILE *fout = fopen(filename,"wb+");

    // write data
    fwrite(scalar,1,mem_size_scalar,fout);

    // close file
    fclose(fout);

    if(ferror(fout))
    {
        fprintf(stderr,"Error saving to %s\n",filename);
        perror("");
    }
    else
    {
        if(!quiet)
            printf("Saved to %s\n",filename);
    }
}

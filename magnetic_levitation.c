#include <math.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

// work in SI units

const double mu0=1.25663706212e-6;
const double pi=3.1415926535897;
/* magnetization */
const double Bm_4_pi = 1.2/4/pi;
const double sigma = 3.5e7;
const double v0 = 1;
/* The metal */
const double Lx=0.1, Ly=0.1, Lz=0.01;
const int Nx = 100, Ny = 100, Nz =10, N=Nx*Ny*Nz;
const double dx=Lx/(Nx-1),dy=Ly/(Ny-1),dz=Lz/(Nz-1);

const double Mx0=.045, Mx1=.055, My0=.04, My1=.06, Mz0=0.05, Mz1=0.055;
const double MNx=10, MNy=20, MNz=5;
const double Mdx=(Mx1-Mx0)/MNx, Mdy=(My1-My0)/MNy, Mdz=(Mz1-Mz0)/MNz;

const double Ox0=0, Ox1=.1, Oy0=0, Oy1=.1, Oz=0.03;
const int NOx=100, NOy=100, NO=NOx*NOy;

double *B0y, *B0z, *By, *Bz, *jy, *jz,*OBz;
double fx,fz;

double max_change;

const double C1 = mu0 * dx * dy * dz / (4 * pi);

#define GET_DIST() do{\
    deltax = x0-x1;\
    deltay = y0-y1;\
    deltaz = z0-z1;\
    dist3 = deltax*deltax+deltay*deltay+deltaz*deltaz;\
    dist3 *= sqrt(dist3);\
}while(0)

#define FOR(i) for(int idx##i=0; idx##i<N; idx##i++){\
    double x##i=(idx##i/(Ny*Nz)) * dx;\
    double y##i=((idx##i/Nz)%Ny) * dy;\
    double z##i=(idx##i%Nz) * dz;

#define END }

void make_grid(){
    size_t s = sizeof(double)*N;
    B0y = malloc(s); 
    B0z = malloc(s); 
    By = malloc(s); 
    Bz = malloc(s); 
    jy = malloc(s); 
    jz = malloc(s);
    OBz = malloc(sizeof(double)*NOx*NOy);
}

void make_B0(){
#pragma omp parallel for
    FOR(0)
	double C, deltaz, deltay, deltax, dist3;
	B0y[idx0] = 0;
	B0z[idx0] = 0;
	for(int iz=0; iz<MNz; iz++){
	    double z1 = (iz+.5) * Mdz + Mz0;
	    for(int ix=0; ix<MNx; ix++){
		double x1 = (ix+.5) * Mdx + Mx0;
		double y1;

		y1 = My0;
		GET_DIST();
		C = Bm_4_pi * Mdx * Mdz / dist3;
		B0y[idx0] -= deltaz * C;
		B0z[idx0] += deltay * C;

		y1 = My1;
		GET_DIST();
		C = Bm_4_pi * Mdx * Mdz / dist3;
		B0y[idx0] += deltaz * C;
		B0z[idx0] -= deltay * C;
	    }
	    for(int iy=0; iy<MNy; iy++){
		double y1 = (iy+.5) * Mdy + My0;
		double x1;

		x1 = Mx1;
		GET_DIST();
		C = Bm_4_pi * Mdz * Mdy / dist3;
		B0z[idx0] -= deltax * C;

		x1 = Mx0;
		GET_DIST();
		C = Bm_4_pi * Mdz * Mdy / dist3;
		B0z[idx0] += deltax * C;
	    }	
	}
	By[idx0] = B0y[idx0];
	Bz[idx0] = B0z[idx0];
    END
}

void iterate(){
    double change;
    int done;
    max_change = -1;
#pragma omp parallel for reduction(max:max_change)
    for(int i=0; i<N; i++){
	change = fabs(- sigma * v0 * Bz[i] - jy[i]);
	max_change = change > max_change ? change : max_change;
	jy[i] = - sigma * v0 * Bz[i];

	change = fabs(sigma * v0 * By[i] - jz[i]);
	max_change = change > max_change ? change : max_change;
	jz[i] =   sigma * v0 * By[i];
    }
#pragma omp parallel for
    FOR(0)
	double C, deltaz, deltay, deltax, dist3;
	By[idx0] = B0y[idx0];
	Bz[idx0] = B0z[idx0];
	FOR(1)
	    if (idx0==idx1) continue;
	    GET_DIST();
	    C = C1 / dist3;
	    By[idx0] += C * jz[idx1] * deltax;
	    Bz[idx0] -= C * jy[idx1] * deltax;
	END
    END
}

/*
void check_divergence(){
#pragma omp parallel for
    FOR(0)

    END
}
*/

void make_plane(){
#pragma omp parallel for
    for(int idx0=0; idx0<NO; idx0++){
	double deltax, deltay, deltaz, dist3;
	double x0 = Ox0 + (Ox1-Ox0) * (idx0/NOy) / NOx;
	double y0 = Oy0 + (Oy1-Oy0) * (idx0%NOy) / NOy;
	double z0 = Oz;
	OBz[idx0] = 0;
	FOR(1)
	    GET_DIST();
	    OBz[idx0] -= dx * dy * dz * mu0/4/pi * jy[idx1] * deltax / dist3;
	END
    }
}

void get_force(){
    fx = fz = 0;
#pragma omp parallel for reduction(+:fx,fz)
    FOR(1)
	double deltaz,deltay,deltax,dist3;
	for(int ix=0; ix<MNx; ix++){
	    double x0 = (ix+.5)*Mdx + Mx0;
	    for(int iy=0; iy<MNy; iy++){
		double y0 = (iy+.5)*Mdy + My0;
		double z0;
		z0 = Mz0;
		GET_DIST();
		fx -= dx * dy * dz * Bm_4_pi * (jy[idx1] * deltaz - jz[idx1] * deltay) / dist3;
		fz += dx * dy * dz * Bm_4_pi * jy[idx1] * deltax / dist3;
		z0 = Mz1;
		GET_DIST();
		fx += dx * dy * dz * Bm_4_pi * (jy[idx1] * deltaz - jz[idx1] * deltay) / dist3;
		fz -= dx * dy * dz * Bm_4_pi * jy[idx1] * deltax / dist3;
	    }
	}
    END
}

int main(){
    make_grid();
    make_B0();
    for(int step=0; step<10; step++){
	iterate();
	printf("step %d max change %lf\n", step, max_change);
	fflush(stdout);
    }
    make_plane();
    get_force();
    printf("fz=%lf fx=%lf\n", fz,fx);
    return 0;
}

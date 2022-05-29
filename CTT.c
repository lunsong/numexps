#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int     n_theta;
double  dtheta;
int     n_r;
double  dr;
int     n_legendre;
double* legendre_matrix;
double  M,J,k;
double  alpha;
double  max_diff;

#define get_theta(i) (i+.5)*dtheta
#define get_r(i) (i+1)*dr

void make_params(){
    dtheta = 1./n_theta;
    dr = 1./(n_r+1);
    k = 9*J*J/(16*M_PI);
}

double _pow(double x, int n){
    int i=1; for(;i<n; i<<=1);
    double y = 1;
    for(;i>1;i>>=1){
	if(n&i) y *= x;
	y *= y;
    }
    if(n&i) y*=x;
    return y;
}

void calculate_legendre_matrix(){
    legendre_matrix = malloc(sizeof(double)*n_legendre*n_theta);
    double* a = malloc(sizeof(double)*(n_legendre+1));
    for(int l=0; l<n_legendre; l++){
        double L = 2*l;
        double c = sqrt(2*L+1);
        if(l==0)
            a[0] = 1;
        else{
            for(int m=0; m<l; m++)
                a[m] = a[m] * (L-m) * (2*(L-m)-1) * (2*(L-m)-2) * (2*(L-m)-3) / ( 2 * (L-m-1) * (L-m) * (L-2*m) * (L-2*m-1) );
            a[l] = - a[l-1] * (l+1) * 2 / (l * (L+2) * (L+1));
        }
        for(int i=0; i<n_theta; i++){
            double x = get_theta(i);
            double y = a[0];
            for(int m=1; m<=l; m++) y = x*x*y + a[m];
            legendre_matrix[l*n_theta+i] = y*c;
        }
    }
    free(a);
}

void legendre_trans(double*x, double*y){
#pragma omp parallel for
    for(int i=0; i<n_r*n_legendre; i++){
        int l = i%n_legendre;
        int ir = i/n_legendre;
        double s = .5 * (x[ir*n_theta]*legendre_matrix[l*n_theta]
                +x[ir*(n_theta+1)-1]*legendre_matrix[l*(n_theta+1)-1]);
        for(int itheta=1; itheta<n_theta-1; itheta++)
            s += x[ir*n_theta + itheta] *
                legendre_matrix[l*n_theta+itheta];
        y[i] = s * dtheta;
    }
}

void legendre_expand(double*x, double*y){
#pragma omp parallel for
    for(int i=0; i<n_r*n_theta; i++){
        int itheta = i%n_theta;
        int ir     = i/n_theta;
        double s = 0;
	for(int l=0; l<n_legendre; l++)
	    s += x[ir*n_legendre+l] * legendre_matrix[l*n_theta+itheta];
        y[i] = s;
    }
}

void nonlinear_f(double *u){
#pragma omp parallel for
    for(int i=0; i<n_r*n_theta; i++){
        int ir = i/n_theta;
        double c = (i%n_theta)*dtheta;
        double r = alpha * (ir+1)*dr / (1-(ir+1)*dr);
        double x = (u[i]+1)*r+M;
        double x3 = x*x*x;
        double x7 = x3*x3*x;
        u[i] = k * r * (1-c*c) / x7;
    }
}

void integrate_r(double *x, double *y){
    max_diff = 0;
#pragma omp parallel for reduction(max:max_diff)
    for(int i=0; i<n_r*n_legendre; i++){
        int l = i%n_legendre;
        int ir = i/n_legendre;
        double xi = get_r(ir);
        double ri = xi/(1-xi);
	double A = _pow(1/ri,l+1);
        double B = _pow(ri,l);
        double s = 0;
        for(int j=0; j<ir; j++){
            double xj = (j+1)*dr;
            s += x[j*n_legendre+l] * _pow(xj/(1-xj), l+2) * A / ((1-xj)*(1-xj));
        }
        for(int j=ir; j<n_r; j++){
            double xj = (j+1)*dr;
            s += x[j*n_legendre+l] * _pow((1-xj)/xj, l-3) * B / (xj*xj);
        }
        s = s * alpha*alpha * dr / (2*l+1);
	double diff = fabs(y[i] - s);
	max_diff = max_diff > diff ? max_diff : diff ;
	y[i] = s;
    }
}

int main(){
    alpha = 5;
    M = 1;
    J = 1;
    n_theta=32;
    n_legendre = 16;
    n_r=128;
    make_params();
    calculate_legendre_matrix();
    double *u = malloc(sizeof(double)*n_r*n_theta);
    for(int i=0; i<n_r*n_theta; i++) u[i] = 0;
    double *ul = malloc(sizeof(double)*n_r*n_legendre);
    double *vl = malloc(sizeof(double)*n_r*n_legendre);
    for(int i=0; i<200; i++){
#if 0
        printf("u\n");
        for(int I=0; I<n_r; I++){
            for(int J=0; J<n_theta; J++)
                printf("%lf ", u[I*n_theta+J]);
            printf("\n");
        }
#endif
        nonlinear_f(u);
#if 0
        printf("f\n");
        for(int I=0; I<n_r; I++){
            for(int J=0; J<n_theta; J++)
                printf("%lf ", u[I*n_theta+J]);
            printf("\n");
        }
#endif
        legendre_trans(u,ul);
#if 0
        printf("ul\n");
        for(int I=0; I<n_r; I++){
            for(int J=0; J<n_theta; J++)
                printf("%lf ", ul[I*n_theta+J]);
            printf("\n");
        }
#endif
        integrate_r(ul,vl);
#if 0
        printf("vl\n");
        for(int I=0; I<n_r; I++){
            for(int J=0; J<n_theta; J++)
                printf("%lf ", vl[I*n_theta+J]);
            printf("\n");
        }
#endif
        legendre_expand(vl,u);
        fprintf(stderr, "%d %e\n", i, max_diff);
    }
    for(int l=0; l<n_legendre; l++){
	for(int ir=0; ir<n_r; ir++) printf("%e ",vl[ir*n_legendre+l]);
	printf("\n");
    }
    return 0;
}

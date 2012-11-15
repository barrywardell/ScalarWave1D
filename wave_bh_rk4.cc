/* wave_bh.cc
 *
 * Solve the second order in space, first order in time spherical harmonic
 * mode decomposed scalar wave equation in Schwarzschild spacetime using
 * RK4 time integration and spatial finite differencing.
 *
 * As an initial condition, choose a gaussian centered about r_* = 10.
 *
 * The spatial resolution is controlled by the parameter N below.
 *
 * Compiling: g++ -lgsl -o wave_bh wave_bh.cc
 *
 * Running: ./wave_bh
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "precision.h"

extern "C" {
#include <gsl/gsl_sf_lambert.h>
}

/* User specified parameters */
const int  N        = 40001;   /* Number of points */
const REAL rstarMax = 4000.0L; /* Size of grid */
const REAL T        = 100L;   /* Final time */
const int  l        = 1;       /* Spherical harmonic mode */
const REAL M        = 1.0L;    /* Mass of the black hole */

/* Derived parameters */
const REAL h = 2.0L*rstarMax/(N-1); /* Grid spacing */
const REAL k = 0.25L*h;             /* Minimum time step size */

/* Calculate the right hand sides */
void rhs(const REAL t, const REAL V[], const REAL phi[], const REAL rho[], REAL phidot[], REAL rhodot[])
{
  /* Calculate interior */
#pragma omp parallel for
  for(int i=2; i<N-2; i++)
  {
    REAL phi_xx = (-phi[i+2] + 16L*phi[i+1] - 30L*phi[i] + 16L*phi[i-1] - phi[i-2])/(12.0L*h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - phi[i] * V[i];
  }

  /* Second last points */
  {
    int i=1;
    REAL phi_xx = (phi[i+1] - 2L*phi[i] + phi[i-1])/(h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - phi[i] * V[i];
  }

  {
    int i=N-2;
    REAL phi_xx = (phi[i+1] - 2L*phi[i] + phi[i-1])/(h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - phi[i] * V[i];
  }
  
  /* Horizon boundary */
  {
	phidot[0] = - (-1.0L)*(-3L*phi[0]+4L*phi[1]-phi[2])/(2.0L*h);
	rhodot[0]  = - (-1.0L)*(-3L*rho[0]+4L*rho[1]-rho[2])/(2.0L*h);
  }
	
  /* Spatial infinity boundary */
  {
	phidot[N-1] = - 1.0L*(phi[N-3]-4L*phi[N-2]+3L*phi[N-1])/(2.0L*h);
	rhodot[N-1]  = - 1.0L*(rho[N-3]-4L*rho[N-2]+3L*rho[N-1])/(2.0L*h);
  }
}

/* Output values of fields */
void output(REAL t, REAL rstar[], REAL phi[], REAL rho[]) {
//   for (int j = 0; j < N; j++) {
//     printf("%.19f\t", phi[j]);
//   }
  int i = (N-1)/2;

#ifdef FP_PRECISION_QUAD
  char* y = new char[1000];
  quadmath_snprintf(y, 1000, "%.19Qg", t);
  printf("%s\t", y);
  quadmath_snprintf(y, 1000, "%.19Qg", rstar[i]);
  printf("%s\t", y);
  quadmath_snprintf(y, 1000, "%.19Qg", phi[i]);
  printf("%s\t", y);
  quadmath_snprintf(y, 1000, "%.19Qg", rho[i]);
  printf("%s", y);
  printf("\n");
  delete y;
#elif defined FP_PRECISION_EXTENDED_DOUBLE
  printf("%.19Lg\t%.19Lg\t%.19Lg\t%.19Lg\n", t, rstar[i], phi[i], rho[i]);
#else
  printf("%.16g\t%.16g\t%.16g\t%.16g\n", t, rstar[i], phi[i], rho[i]);
#endif
}

int main()
{
  /* Setup grid */
  REAL rstar0 = 12.7725887222397812376689284858327062723020005374410210164827L;
  REAL *V     = new REAL[N];
  REAL *rstar = new REAL[N];

//   FILE * potential;
//   potential = fopen ("potential.dat","r");
  for (int j = 0; j < N; j++) {
    rstar[j] = rstar0-rstarMax + j * h;
    
    /* Pre-compute the potential */
    REAL rm2M = 2.0L*M*gsl_sf_lambert_W0(exp(-1.0L+rstar[j]/(2.0L*M)));
    REAL r = rm2M + 2.0L*M;
	V[j] = rm2M*(l*(l+1)*r + 2.0L*M)/(r*r*r*r);
// #ifdef FP_PRECISION_QUAD
//     char tmp[1000];
// 	fscanf(potential, "%s\n", tmp);
// 	V[j] = strtoflt128 (tmp, NULL);
// #elif defined FP_PRECISION_EXTENDED_DOUBLE
// 	fscanf(potential, "%Lg\n", &V[j]);
// #else
// 	fscanf(potential, "%lg\n", &V[j]);
// #endif
  }
//   fclose (potential);

  REAL *phi = new REAL[N];
  REAL *rho = new REAL[N];

  /* Initial Conditions: a gaussian about r=10 */
  for (int j = 0; j < N; j++) {
    phi[j] = EXP(-(rstar[j]-rstar0)*(rstar[j]-rstar0)/2.0L);
    rho[j] = -(rstar[j]-rstar0)*phi[j];
  }

  /* Allocate scratch space for integration */
  REAL *phidot = new REAL[N];
  REAL *rhodot = new REAL[N];
  REAL *phi_tmp = new REAL[N];
  REAL *rho_tmp = new REAL[N];
  
  REAL *k1_phi = new REAL[N];
  REAL *k1_rho = new REAL[N];
  REAL *k2_phi = new REAL[N];
  REAL *k2_rho = new REAL[N];
  REAL *k3_phi = new REAL[N];
  REAL *k3_rho = new REAL[N];
  REAL *k4_phi = new REAL[N];
  REAL *k4_rho = new REAL[N];

  /* Fourth order Runge-Kutta time integration */
  for (REAL t = 0; t <= T; t+=k)
  {
    /* Output the results */
    output(t, rstar, phi, rho);

  	/* Calculate k1 term in RK4 formula */
    rhs(t, V, phi, rho, phidot, rhodot);

#pragma omp parallel for
    for (int j = 0; j < N; j++)
    {
      k1_phi[j] = k * phidot[j];
      k1_rho[j] = k * rhodot[j];
    }

    /* Calculate k2 term in RK4 formula */
#pragma omp parallel for
    for (int j = 0; j < N; j++)
    {
      phi_tmp[j] = phi[j] + 0.5L*k1_phi[j];
      rho_tmp[j] = rho[j] + 0.5L*k1_rho[j];
    }

    rhs(t + 0.5L*k, V, phi_tmp, rho_tmp, phidot, rhodot);

#pragma omp parallel for
    for (int j = 0; j < N; j++)
    {
      k2_phi[j] = k * phidot[j];
      k2_rho[j] = k * rhodot[j];
    }

    /* Calculate k3 term in RK4 formula */
#pragma omp parallel for
    for (int j = 0; j < N; j++)
    {
      phi_tmp[j] = phi[j] + 0.5L*k2_phi[j];
      rho_tmp[j] = rho[j] + 0.5L*k2_rho[j];
    }

    rhs(t + 0.5L*k, V, phi_tmp, rho_tmp, phidot, rhodot);

#pragma omp parallel for
    for (int j = 0; j < N; j++)
    {
      k3_phi[j] = k * phidot[j];
      k3_rho[j] = k * rhodot[j];
    }

    /* Calculate k4 term in RK4 formula */
#pragma omp parallel for
    for (int j = 0; j < N; j++)
    {
      phi_tmp[j] = phi[j] + k3_phi[j];
      rho_tmp[j] = rho[j] + k3_rho[j];
    }

    rhs(t + k, V, phi_tmp, rho_tmp, phidot, rhodot);

#pragma omp parallel for
    for (int j = 0; j < N; j++)
    {
      k4_phi[j] = k * phidot[j];
      k4_rho[j] = k * rhodot[j];
    }

	/* Apply RK4 formula to update variables */
#pragma omp parallel for
    for (int j = 0; j < N; j++)
    {
      phi[j] = phi[j] + (k1_phi[j]+2.0L*k2_phi[j]+2.0L*k3_phi[j]+k4_phi[j])/6.0L;
      rho[j] = rho[j] + (k1_rho[j]+2.0L*k2_rho[j]+2.0L*k3_rho[j]+k4_rho[j])/6.0L;
    }
  }
  
  delete phi;
  delete rho;
  delete rstar;
  delete V;

  delete phidot;
  delete rhodot;
  delete phi_tmp;
  delete rho_tmp;

  delete k1_phi;
  delete k1_rho;
  delete k2_phi;
  delete k2_rho;
  delete k3_phi;
  delete k3_rho;
  delete k4_phi;
  delete k4_rho;

  return 0;
}
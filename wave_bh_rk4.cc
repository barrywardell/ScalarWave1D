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
extern "C" {
#include <quadmath.h>
}

/* User specified parameters */
const int    N = 40001;				/* Number of points */
const __float128 rstarMax = 4000.0L;		/* Size of grid */
const __float128 T = 8000L;				/* Final time */
const int    l = 1;					/* Spherical harmonic mode */
const __float128 M = 1.0L;				/* Mass of the black hole */
const int scaling = 0;

/* Derived parameters */
const __float128 h = 2.0L*rstarMax/(N-1); /* Grid spacing */
const __float128 k = 0.25L*h;			 /* Minimum time step size */

/* Calculate the right hand sides */
void rhs(const __float128 t, const __float128 V[], const __float128 phi[], const __float128 rho[], __float128 phidot[], __float128 rhodot[])
{
  __float128 rho_coeff, phi_coeff;
  __float128 t_scaling_2, t_scaling_1, t_scaling;

  __float128 t3 = t*t*t;
  __float128 t4 = t3*t;
  __float128 t5 = t4*t;

  switch(scaling) {
  case 0:
    phi_coeff = 0.0L;
    rho_coeff = 0.0L;
    break;
  case 1:
    phi_coeff = -2.0L / ((1.0L + t) * (1.0L + t));
    rho_coeff = 2.0L / (1.0L + t);
    break;
  case 2:
    phi_coeff = (2.0L - 6.0L*t*t) / ((1.0L + t) * (1.0L + t));
    rho_coeff = 4.0L * t / (1.0L + t*t);
    break;
  case 5:
    phi_coeff = -(10.0L*t3)*(3.0L*t5-2.0L) / ((1.0L + t5) * (1.0L + t5));
    rho_coeff = 10.0L * t4 / (1.0L + t5);
    break;
  default:
    t_scaling_2 = pow(t, scaling-2.0L);
    t_scaling_1 = t*t_scaling_2;
    t_scaling   = t*t_scaling_1;
    phi_coeff = - scaling * t_scaling_2 * (1.0L + t_scaling + scaling * (t_scaling - 1.0L)) /
        ((1.0L + t_scaling) * (1.0L + t_scaling));
    rho_coeff = 2.0L * scaling * t_scaling_1 / (1.0L + t_scaling);
    break;
  }

  /* Calculate interior */
#pragma omp parallel for
  for(int i=2; i<N-2; i++)
  {
    __float128 phi_xx = (-phi[i+2] + 16L*phi[i+1] - 30L*phi[i] + 16L*phi[i-1] - phi[i-2])/(12.0L*h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - phi[i] * V[i]
      + rho[i] * rho_coeff
      + phi[i] * phi_coeff;
  }

  /* Second last points */
  {
    int i=1;
    __float128 phi_xx = (phi[i+1] - 2L*phi[i] + phi[i-1])/(h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - phi[i] * V[i]
      + rho[i] * rho_coeff
      + phi[i] * phi_coeff;
  }

  {
    int i=N-2;
    __float128 phi_xx = (phi[i+1] - 2L*phi[i] + phi[i-1])/(h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - phi[i] * V[i]
      + rho[i] * rho_coeff
      + phi[i] * phi_coeff;
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
void output(__float128 t, __float128 rstar[], __float128 phi[], __float128 rho[]) {
//   for (int j = 0; j < N; j++) {
//     printf("%.19f\t", phi[j]);
//   }
  int i = (N-1)/2;
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
}

int main()
{
  /* Setup grid */
  __float128 rstar0 = 12.7725887222397812376689284858327062723020005374410210164827L;
  __float128 *V     = new __float128[N];
  __float128 *rstar = new __float128[N];

  FILE * potential;
  potential = fopen ("potential.dat","r");
  for (int j = 0; j < N; j++) {
    rstar[j] = rstar0-rstarMax + j * h;
    char tmp[1000];
    
    /* Pre-compute the potential */
    /*__float128 rm2M = 2.0L*M*gsl_sf_lambert_W0(exp(-1.0L+rstar[j]/(2.0L*M)));
    __float128 r = rm2M + 2.0L*M;
	V[j] = rm2M*(l*(l+1)*r + 2.0L*M)/(r*r*r*r);*/
	fscanf(potential, "%s\n", tmp);
	V[j] = strtoflt128 (tmp, NULL);
  }
  fclose (potential);

  __float128 *phi = new __float128[N];
  __float128 *rho = new __float128[N];

  /* Initial Conditions: a gaussian about r=10 */
  for (int j = 0; j < N; j++) {
    phi[j] = expq(-(rstar[j]-rstar0)*(rstar[j]-rstar0)/2.0L);
    rho[j] = -(rstar[j]-rstar0)*phi[j];
  }

  /* Allocate scratch space for integration */
  __float128 *phidot = new __float128[N];
  __float128 *rhodot = new __float128[N];
  __float128 *phi_tmp = new __float128[N];
  __float128 *rho_tmp = new __float128[N];
  
  __float128 *k1_phi = new __float128[N];
  __float128 *k1_rho = new __float128[N];
  __float128 *k2_phi = new __float128[N];
  __float128 *k2_rho = new __float128[N];
  __float128 *k3_phi = new __float128[N];
  __float128 *k3_rho = new __float128[N];
  __float128 *k4_phi = new __float128[N];
  __float128 *k4_rho = new __float128[N];

  /* Fourth order Runge-Kutta time integration */
  for (__float128 t = 0; t <= T; t+=k)
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
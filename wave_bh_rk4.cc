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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_lambert.h>

/* User specified parameters */
const int    N = 20001;				/* Number of points */
const long double rstarMax = 1000.0L;		/* Size of grid */
const long double T = 2000L;				/* Final time */
const int    l = 1;					/* Spherical harmonic mode */
const long double M = 1.0L;				/* Mass of the black hole */

/* Derived parameters */
const long double h = 2.0L*rstarMax/(N-1); /* Grid spacing */
const long double k = 0.25L*h;			 /* Minimum time step size */

/* Calculate the right hand sides */
void rhs(const long double V[], const long double phi[], const long double rho[], long double phidot[], long double rhodot[])
{
  /* Calculate interior */
  for(int i=2; i<N-2; i++)
  {
    long double phi_xx = (-phi[i+2] + 16L*phi[i+1] - 30L*phi[i] + 16L*phi[i-1] - phi[i-2])/(12.0L*h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - V[i]*phi[i];
  }

  /* Second last points */
  {
    int i=1;
    long double phi_xx = (phi[i+1] - 2L*phi[i] + phi[i-1])/(h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - V[i]*phi[i];
  }

  {
    int i=N-2;
    long double phi_xx = (phi[i+1] - 2L*phi[i] + phi[i-1])/(h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - V[i]*phi[i];
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
void output(long double t, long double rstar[], long double phi[], long double rho[]) {
//   for (int j = 0; j < N; j++) {
//     printf("%.19f\t", phi[j]);
//   }
  int i = (N-1)/2;
  printf("%.19Lg\t%.19Lg\t%.19Lg\t%.19Lg\t", t, rstar[i], phi[i], rho[i]);
  printf("\n");
}

int main()
{
  /* Setup grid */
  long double rstar0 = 12.7725887222397812376689284858327062723020005374410210164827L;
  long double *V     = new long double[N];
  long double *rstar = new long double[N];
  for (int j = 0; j < N; j++) {
    rstar[j] = rstar0-rstarMax + j * h;
    
    /* Pre-compute the potential */
    long double rm2M = 2.0L*M*gsl_sf_lambert_W0(exp(-1.0L+rstar[j]/(2.0L*M)));
    long double r = rm2M + 2.0L*M;
	V[j] = rm2M*(l*(l+1)*r + 2.0L*M)/(r*r*r*r);
  }

  long double *phi = new long double[N];
  long double *rho = new long double[N];

  /* Initial Conditions: a gaussian about r=10 */
  long double t = 0;
  for (int j = 0; j < N; j++) {
    phi[j] = exp(-(rstar[j]-rstar0)*(rstar[j]-rstar0)/2.0L);
    rho[j] = -(rstar[j]-rstar0)*phi[j];
  }

  /* Output the initial data */
  output(t, rstar, phi, rho);

  /* Allocate scratch space for integration */
  long double *phidot = new long double[N];
  long double *rhodot = new long double[N];
  long double *phi_tmp = new long double[N];
  long double *rho_tmp = new long double[N];
  
  long double *k1_phi = new long double[N];
  long double *k1_rho = new long double[N];
  long double *k2_phi = new long double[N];
  long double *k2_rho = new long double[N];
  long double *k3_phi = new long double[N];
  long double *k3_rho = new long double[N];
  long double *k4_phi = new long double[N];
  long double *k4_rho = new long double[N];

  /* Fourth order Runge-Kutta time integration */
  for (long double t = 0; t <= T; t+=k)
  {
  	/* Calculate k1 term in RK4 formula */
  	rhs(V, phi, rho, phidot, rhodot);

    for (int j = 0; j < N; j++)
    {
      k1_phi[j] = k * phidot[j];
      k1_rho[j] = k * rhodot[j];
    }

    /* Calculate k2 term in RK4 formula */
    for (int j = 0; j < N; j++)
    {
      phi_tmp[j] = phi[j] + 0.5L*k1_phi[j];
      rho_tmp[j] = rho[j] + 0.5L*k1_rho[j];
    }

  	rhs(V, phi_tmp, rho_tmp, phidot, rhodot);

    for (int j = 0; j < N; j++)
    {
      k2_phi[j] = k * phidot[j];
      k2_rho[j] = k * rhodot[j];
    }

    /* Calculate k3 term in RK4 formula */
    for (int j = 0; j < N; j++)
    {
      phi_tmp[j] = phi[j] + 0.5L*k2_phi[j];
      rho_tmp[j] = rho[j] + 0.5L*k2_rho[j];
    }

  	rhs(V, phi_tmp, rho_tmp, phidot, rhodot);

    for (int j = 0; j < N; j++)
    {
      k3_phi[j] = k * phidot[j];
      k3_rho[j] = k * rhodot[j];
    }

    /* Calculate k4 term in RK4 formula */
    for (int j = 0; j < N; j++)
    {
      phi_tmp[j] = phi[j] + k3_phi[j];
      rho_tmp[j] = rho[j] + k3_rho[j];
    }

  	rhs(V, phi_tmp, rho_tmp, phidot, rhodot);

    for (int j = 0; j < N; j++)
    {
      k4_phi[j] = k * phidot[j];
      k4_rho[j] = k * rhodot[j];
    }

	/* Apply RK4 formula to update variables */
    for (int j = 0; j < N; j++)
    {
      phi[j] = phi[j] + (k1_phi[j]+2.0L*k2_phi[j]+2.0L*k3_phi[j]+k4_phi[j])/6.0L;
      rho[j] = rho[j] + (k1_rho[j]+2.0L*k2_rho[j]+2.0L*k3_rho[j]+k4_rho[j])/6.0L;
    }
    
    /* Output the results */
    output(t, rstar, phi, rho);
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
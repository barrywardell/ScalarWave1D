/* wave_bh_rk2.cc
 *
 * Solve the spherical harmonic mode-decomposed scalar wave equation in
 * Schwarzschild spacetime using RK2 (midpoint method) and second order
 * centered finite differencing.
 *
 * As an initial condition, choose a gaussian centered about r_* = 10.
 *
 * For simplicity, place boundaries far away so that we can ignore them.
 *
 * The resolution is controlled by the parameter N below.
 *
 * Compiling: g++ -lgsl -o wave_bh_rk2 wave_bh_rk2.cc
 *
 * Running: ./wave_bh_rk2
 *
 */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_lambert.h>

const int    N = 6002;				/* Number of points */
const double rstarMax = 1000.0;		/* Size of grid */
const double h = 2.0*rstarMax/(N-2);/* Grid spacing */
const double k = 0.065*h;			/* Time step size */
const double T = 2000;				/* Final time */
const int    l = 1;					/* Spherical harmonic mode */
const double M = 1.0;				/* Mass of the black hole */

/* Output values of phi and pi */
void output(double rstar[], double phi[], double pi[]) {
  for (int j = 0; j < N-1; j++) {
    printf("%.19f\t%.19f\t%.19f\n", rstar[j], phi[j], pi[j]);
  }
  printf("\n");
}

/* Calculate the right hand sides */
void f(double *phi, double *pi, double *r, double *phidot, double *pidot, int j) {
  phidot[j] = pi[j];
  pidot[j]  = (phi[j+1] - 2*phi[j] + phi[j-1])/(h*h)
    - phi[j]*(r[j]-2.0*M)/r[j]*(l*(l+1)*r[j] + 2.*M)/pow(r[j],3);
}

int main()
{
  double phi[N], pi[N], r[N], rstar[N];
  double rstar0 = 12.772588722239782;

  /* Initial data - a gaussian about rstar=10 */
  for (int j = 0; j < N; j++) {
    rstar[j]   = rstar0-rstarMax + j * h;
    r[j] = 2.0*M*(1.0 + gsl_sf_lambert_W0(exp(-1.0+rstar[j]/(2.0*M))));
    phi[j] = exp(-(rstar[j]-rstar0)*(rstar[j]-rstar0)/2.0);
    pi[j]  = -(rstar[j]-rstar0)*phi[j];
  }

  /* Output the initial data */
  // output(rstar, phi, pi);
  printf("0\t%.19f\t%.19f\n", phi[3000], pi[3000]);
  double k1phi[N], phipk1phi[N], k2phi[N]; 
  double k1pi[N],  pipk1pi[N],   k2pi[N]; 
  double phidot[N], pidot[N];

  /* Second order Runge-Kutta time integration */
  for (double t = 0; t <= T; t+=k) {

	/* Calculate k1 term in RK2 formula */
    for (int j = 1; j < N-1; j++)
    {
      f(phi, pi, r, phidot, pidot, j);
      k1phi[j] = k * phidot[j];
      k1pi[j]  = k * pidot[j];
    }

	/* Flat boundary conditions */
    k1phi[0]   = k1phi[1];
    k1phi[N-1] = k1phi[N-2];

    k1pi[0]   = k1pi[1];
    k1pi[N-1] = k1pi[N-2];

    /* Calculate k2 term in RK2 formula */
    for (int j = 0; j < N; j++)
    {
      phipk1phi[j] = phi[j] + 0.5*k1phi[j];
      pipk1pi[j] = pi[j] + 0.5*k1pi[j];
    }

    for (int j = 1; j < N-1; j++)
    {
      f(phipk1phi, pipk1pi, r, phidot, pidot, j);
      k2phi[j] = k * phidot[j];
      k2pi[j]  = k * pidot[j];
    }

	/* Flat boundary conditions */
    k2phi[0]   = k2phi[1];
    k2phi[N-1] = k2phi[N-2];

    k2pi[0]   = k2pi[1];
    k2pi[N-1] = k2pi[N-2];

	/* Apply RK2 formula to update variables */
    for (int j = 0; j < N; j++)
    {
      phi[j] = phi[j] + k2phi[j];
      pi[j] = pi[j] + k2pi[j];
    }
    
    // output(rstar, phi, pi);
	printf("%.19f\t%.19f\n", t, phi[3000], pi[3000]);
  }
}

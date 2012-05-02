/* wave_bh.cc
 *
 * Solve the second order in space, first order in time spherical harmonic
 * mode decomposed scalar wave equation in Schwarzschild spacetime using GSL
 * time integration and spatial finite differencing.
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_lambert.h>

const int    N = 301;				/* Number of points */
const double rstarMax = 50.0;		/* Size of grid */
const double h = 2.0*rstarMax/N;    /* Grid spacing */
const double k0 = 0.5*h;			/* Minimum time step size */
const double T = 100;				/* Final time */
const int    l = 1;					/* Spherical harmonic mode */
const double M = 1.0;				/* Mass of the black hole */

/* Parameters which will be passed to the RHS calculation by the GSL */
struct wave_params {
  const double *r;
  const double *rstar;
  double l; /* Spherical harmonic mode */
  double M; /* Black hole mass */
};

/* Define the potential */
static double V(int l, double r, double M)
{
  return (r-2.0*M)*(l*(l+1)*r + 2.*M)/gsl_pow_4(r);
}

/* Calculate the right hand sides */
int rhs(double tau, const double y[], double f[], void *params) {
  /* Unpack the data */
  const double * phi = y;
  const double * pi  = y + N;

  double *phidot = f;
  double *pidot  = f + N;

  struct wave_params p = *((struct wave_params*) params);
  const double *r = p.r;
  const double *rstar = p.rstar;
  const int l = p.l;
  const double M = p.M;

  /* Calculate interior */
  for(int i=1; i<N-1; i++)
  {
    double phi_xx = (phi[i+1] - 2*phi[i] + phi[i-1])/(h*h);
    phidot[i] = pi[i];
    pidot[i]  = phi_xx - V(l, r[i], M)*phi[i];
  }
  
  /* Horizon boundary */
  {
	phidot[0]   = 0;
	pidot[0]    = 0;
  }
	
  /* Spatial infinity boundary */
  {
	phidot[N-1]   = 0;
	pidot[N-1]    = 0;
  }

  return GSL_SUCCESS;
}

/* Output values of fields */
void output(double rstar[], double phi[], double pi[]) {
  for (int j = 0; j < N; j++) {
    //printf("%.19f\t%.19f\t%.19f\n", rstar[j], phi[j], pi[j]);
    printf("%.19f\t", phi[j]);
  }
  printf("\n");
}

int main()
{

  /* Use a Runge-Kutta integrator with adaptive step-size */
  const gsl_odeiv_step_type * type = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s 			= gsl_odeiv_step_alloc (type, 2*N);
  gsl_odeiv_control * c 		= gsl_odeiv_control_standard_new (1e-6, 1e-6, 1.0, 1.0);
  gsl_odeiv_evolve * e 			= gsl_odeiv_evolve_alloc (2*N);

  /* Setup grid */
  double rstar0 = 12.772588722239782;
  double *r     = new double[N];
  double *rstar = new double[N];
  for (int j = 0; j < N; j++) {
    rstar[j] = rstar0-rstarMax + j * h;
    r[j] = 2.0*M*(1.0 + gsl_sf_lambert_W0(exp(-1.0+rstar[j]/(2.0*M))));
  }
  
  double *y   = new double[2*N];
  double *phi = y;
  double *pi  = y + N;

  /* Parameters for the evolution */
  struct wave_params params = {r, rstar, l, M};
  gsl_odeiv_system sys = {rhs, NULL, 2*N, &params};
  
  double k = k0;

  /* Initial Conditions: a gaussian about r=10 */
  double t = 0;
  for (int j = 0; j < N; j++) {
    phi[j] = exp(-(rstar[j]-rstar0)*(rstar[j]-rstar0)/2.0);
    pi[j]  = -(rstar[j]-rstar0)*phi[j];
  }

  /* Output the initial data */
  output(rstar, phi, pi);

  /* GSL time integration */
  while (t < T)
  {
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, T, &k, y);

    if (status != GSL_SUCCESS)
      break;

    /* Output the results */
    output(rstar, phi, pi);

    /* Enforce CFL condition */
    if (k > k0)
    {
      k=k0;
    }
    
    /* Keep a fixed time step */
    if (k < k0)
    {
      k=k0;
    }
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  
  delete y;
  delete r;
  delete rstar;

  return 0;
}
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

/* User specified parameters */
const int    N = 20001;				/* Number of points */
const double rstarMax = 1000.0;		/* Size of grid */
const double T = 2000;				/* Final time */
const int    l = 1;					/* Spherical harmonic mode */
const double M = 1.0;				/* Mass of the black hole */

/* Derived parameters */
const double h  = 2.0*rstarMax/(N-1); /* Grid spacing */
const double k0 = 0.25*h;			  /* Minimum time step size */

/* Parameters which will be passed to the RHS calculation by the GSL */
struct wave_params {
  const double *rm2M;
  double l; /* Spherical harmonic mode */
  double M; /* Black hole mass */
};

/* Define the potential */
static double V(int l, double rm2M, double M)
{
  double r = rm2M + 2.0*M;
  return rm2M*(l*(l+1)*r + 2.*M)/gsl_pow_4(r);
}

/* Calculate the right hand sides */
int rhs(double tau, const double y[], double f[], void *params) {
  /* Unpack the data */
  const double * phi = y;
  const double * pi  = y + N;

  double *phidot = f;
  double *pidot  = f + N;

  struct wave_params p = *((struct wave_params*) params);
  const double *rm2M = p.rm2M;
  const int l = p.l;
  const double M = p.M;

  /* Calculate interior */
  for(int i=2; i<N-2; i++)
  {
    double phi_xx = (-phi[i+2] + 16*phi[i+1] - 30*phi[i] + 16*phi[i-1] - phi[i-2])/(12.0*h*h);
    phidot[i] = pi[i];
    pidot[i]  = phi_xx - V(l, rm2M[i], M)*phi[i];
  }

  /* Second last points */
  {
    int i=1;
    double phi_xx = (phi[i+1] - 2*phi[i] + phi[i-1])/(h*h);
    phidot[i] = pi[i];
    pidot[i]  = phi_xx - V(l, rm2M[i], M)*phi[i];
  }

  {
    int i=N-2;
    double phi_xx = (phi[i+1] - 2*phi[i] + phi[i-1])/(h*h);
    phidot[i] = pi[i];
    pidot[i]  = phi_xx - V(l, rm2M[i], M)*phi[i];
  }
  
  /* Horizon boundary */
  {
	phidot[0] = - (-1.0)*(-3*phi[0]+4*phi[1]-phi[2])/(2.0*h);
	pidot[0]  = - (-1.0)*(-3*pi[0]+4*pi[1]-pi[2])/(2.0*h);
  }
	
  /* Spatial infinity boundary */
  {
	phidot[N-1] = - 1.0*(phi[N-3]-4*phi[N-2]+3*phi[N-1])/(2.0*h);
	pidot[N-1]  = - 1.0*(pi[N-3]-4*pi[N-2]+3*pi[N-1])/(2.0*h);
  }

  return GSL_SUCCESS;
}

/* Output values of fields */
void output(double t, double rstar[], double phi[], double pi[]) {
//   for (int j = 0; j < N; j++) {
//     printf("%.19f\t", phi[j]);
//   }
  int i = (N-1)/2;
  printf("%.19g\t%.19g\t%.19g\t%.19g\t", t, rstar[i], phi[i], pi[i]);
  printf("\n");
}

int main()
{

  /* Use a Runge-Kutta integrator with adaptive step-size */
  const gsl_odeiv_step_type * type = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s 			= gsl_odeiv_step_alloc (type, 2*N);
  gsl_odeiv_control * c 		= gsl_odeiv_control_standard_new (1e-10, 1e-10, 1.0, 1.0);
  gsl_odeiv_evolve * e 			= gsl_odeiv_evolve_alloc (2*N);

  /* Setup grid */
  double rstar0 = 12.772588722239782;
  double *rm2M  = new double[N];
  double *rstar = new double[N];
  for (int j = 0; j < N; j++) {
    rstar[j] = rstar0-rstarMax + j * h;
    rm2M[j]  = 2.0*M*gsl_sf_lambert_W0(exp(-1.0+rstar[j]/(2.0*M)));
  }
  
  double *y   = new double[2*N];
  double *phi = y;
  double *pi  = y + N;

  /* Parameters for the evolution */
  struct wave_params params = {rm2M, l, M};
  gsl_odeiv_system sys = {rhs, NULL, 2*N, &params};
  
  double k = k0;

  /* Initial Conditions: a gaussian about r=10 */
  double t = 0;
  for (int j = 0; j < N; j++) {
    phi[j] = exp(-(rstar[j]-rstar0)*(rstar[j]-rstar0)/2.0);
    pi[j]  = -(rstar[j]-rstar0)*phi[j];
  }

  /* Output the initial data */
  output(t, rstar, phi, pi);

  /* GSL time integration */
  while (t < T)
  {
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, T, &k, y);

    if (status != GSL_SUCCESS)
      break;

    /* Output the results */
    output(t, rstar, phi, pi);

    /* Enforce CFL condition */
    if (k > k0)
    {
      k=k0;
    }
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  
  delete y;
  delete rm2M;
  delete rstar;

  return 0;
}
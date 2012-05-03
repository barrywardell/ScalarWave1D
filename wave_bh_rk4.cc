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
const int    N = 1001;				/* Number of points */
const double rstarMax = 100.0;		/* Size of grid */
const double T = 200;				/* Final time */
const int    l = 1;					/* Spherical harmonic mode */
const double M = 1.0;				/* Mass of the black hole */

/* Derived parameters */
const double h = 2.0*rstarMax/(N-1); /* Grid spacing */
const double k = 0.25*h;			 /* Minimum time step size */

/* Calculate the right hand sides */
void rhs(const double V[], const double phi[], const double rho[], double phidot[], double rhodot[])
{
  /* Calculate interior */
  for(int i=2; i<N-2; i++)
  {
    double phi_xx = (-phi[i+2] + 16*phi[i+1] - 30*phi[i] + 16*phi[i-1] - phi[i-2])/(12.0*h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - V[i]*phi[i];
  }

  /* Second last points */
  {
    int i=1;
    double phi_xx = (phi[i+1] - 2*phi[i] + phi[i-1])/(h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - V[i]*phi[i];
  }

  {
    int i=N-2;
    double phi_xx = (phi[i+1] - 2*phi[i] + phi[i-1])/(h*h);
    phidot[i] = rho[i];
    rhodot[i]  = phi_xx - V[i]*phi[i];
  }
  
  /* Horizon boundary */
  {
	phidot[0] = - (-1.0)*(-3*phi[0]+4*phi[1]-phi[2])/(2.0*h);
	rhodot[0]  = - (-1.0)*(-3*rho[0]+4*rho[1]-rho[2])/(2.0*h);
  }
	
  /* Spatial infinity boundary */
  {
	phidot[N-1] = - 1.0*(phi[N-3]-4*phi[N-2]+3*phi[N-1])/(2.0*h);
	rhodot[N-1]  = - 1.0*(rho[N-3]-4*rho[N-2]+3*rho[N-1])/(2.0*h);
  }
}

/* Output values of fields */
void output(double t, double rstar[], double phi[], double rho[]) {
//   for (int j = 0; j < N; j++) {
//     printf("%.19f\t", phi[j]);
//   }
  int i = (N-1)/2;
  printf("%.19g\t%.19g\t%.19g\t%.19g\t", t, rstar[i], phi[i], rho[i]);
  printf("\n");
}

int main()
{
  /* Setup grid */
  double rstar0 = 12.772588722239782;
  double *V     = new double[N];
  double *rstar = new double[N];
  for (int j = 0; j < N; j++) {
    rstar[j] = rstar0-rstarMax + j * h;
    
    /* Pre-compute the potential */
    double rm2M = 2.0*M*gsl_sf_lambert_W0(exp(-1.0+rstar[j]/(2.0*M)));
    double r = rm2M + 2.0*M;
	V[j] = rm2M*(l*(l+1)*r + 2.*M)/gsl_pow_4(r);
  }

  double *phi = new double[N];
  double *rho = new double[N];

  /* Initial Conditions: a gaussian about r=10 */
  double t = 0;
  for (int j = 0; j < N; j++) {
    phi[j] = exp(-(rstar[j]-rstar0)*(rstar[j]-rstar0)/2.0);
    rho[j] = -(rstar[j]-rstar0)*phi[j];
  }

  /* Output the initial data */
  output(t, rstar, phi, rho);

  /* Allocate scratch space for integration */
  double *phidot = new double[N];
  double *rhodot = new double[N];
  double *phi_tmp = new double[N];
  double *rho_tmp = new double[N];
  
  double *k1_phi = new double[N]; 
  double *k1_rho = new double[N];
  double *k2_phi = new double[N];
  double *k2_rho = new double[N];  
  double *k3_phi = new double[N]; 
  double *k3_rho = new double[N];
  double *k4_phi = new double[N];
  double *k4_rho = new double[N];

  /* Fourth order Runge-Kutta time integration */
  for (double t = 0; t <= T; t+=k)
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
      phi_tmp[j] = phi[j] + 0.5*k1_phi[j];
      rho_tmp[j] = rho[j] + 0.5*k1_rho[j];
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
      phi_tmp[j] = phi[j] + 0.5*k2_phi[j];
      rho_tmp[j] = rho[j] + 0.5*k2_rho[j];
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
      phi[j] = phi[j] + (k1_phi[j]+2.0*k2_phi[j]+2*k3_phi[j]+k4_phi[j])/6.0;
      rho[j] = rho[j] + (k1_rho[j]+2.0*k2_rho[j]+2*k3_rho[j]+k4_rho[j])/6.0;
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
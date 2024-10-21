#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
//#include "view.h"
#include "two-phase.h"
#include "fractions.h"
#include "curvature.h"
#include "tension.h"
#include "heights.h"
#include "reduced.h"
#include "tracer.h"
//#include "diffusion.h"
#include <math.h>
// Defined parameters (Default values)
#define EPS 1e-14        // Small number
#define PI 3.14159       // Numerical value of Pi
#define RESOLUTION 100   // Grid resolution
#define H0 0.            // Initial interface height, half filled
#define DEPTH H0+0.5     // Depth of fluid
#define rho_ratio   (rho_gas/rho_liquid) // density ratio
#define mu_ratio    (mu_gas/mu_liquid)   // viscosity ratio 

// Dimensional parametrs to be given as cmd line args
double sig;         // Surface tension coefficient 
double mu_gas;   // viscosity of gas
double rho_gas;  // density of gas 
double omega;        // Shaking frequency
double width;       // Size of the box
double rho_liquid;       // Density of water
double mu_liquid;       // Viscosity of water
double grav;       // Acceleration due to gravity
int level;      // Resolution level
double t_end;   // Simulation run time in dimensionless units

// Weber number
double We; 

// Reynolds number
double Re;

// Forcing parameter
double F;

// Tracer field
scalar T[];
scalar * tracers = {T};

// Initial/Max sigma_sq
double sigma_sq_init=0;

// Other variables needed for mixing coefficient calculation
double ka = 0;
double sum_sq = 0;
double sum_mean = 0;
double sigma_sq = 0;

// Code to find interfacial height taken from http://basilisk.fr/sandbox/oaholroyd/controlled_film.c
#define NP 30
vector hei[]; // heights
double interfacial_height(double xp) {
  double dh[NP];
  double yh[NP];
  double yp;
  double y0 = -0.4, y1 = 0.4;

  /* try and find range of possible heights */
  for (int i = 0; i < NP; i++) {
    yp = y0 + i*(y1-y0)/(NP-1);
    Point point = locate(xp,yp);

    if (hei.y[] != nodata) {
      yh[i] = y + height(hei.y[])*Delta;
      dh[i] = fabs(y-yh[i]);
      return yh[i];
    } else {
      yh[i] = -1000;
      dh[i] = 1000;
    }
  } // i end

  /* find the closest one */
  int j = 0;
  for (int i = 1; i < NP; i++) {
    if (dh[i] < dh[j]) { j = i; }
  } // i end

  return yh[j];
}

int main(int argc, char * argv[])
{
  size (1 [0]);
  DT = 1. [0];
  origin (-0.5,-0.5);
  TOLERANCE=1e-5;

  // Import command line parameters
  rho_liquid = atof(argv[1]);
  rho_gas = atof(argv[2]);
  mu_liquid = atof(argv[3]);
  mu_gas = atof(argv[4]);
  sig = atof(argv[5]);
  omega = atof(argv[6]);
  width = atof(argv[7]);
  grav = atof(argv[8]);
  level = atof(argv[9]);
  t_end = atof(argv[10]);

  // Calculate dimensionless numbers
  Re = rho_liquid * omega * width * width / mu_liquid;
  We = rho_liquid * width * omega * width * omega / sig;
  F = width * omega * omega / grav;

  N=pow(2,level);
  // Calculate the viscosities from Reynolds number
  // Has to be defined in main
  mu1 = 1./Re;
  mu2 = mu_ratio*mu1;
  f.sigma = 1./We;

  // Set the densities of the two different fluids
  rho1 = 1.;
  rho2 = rho1*rho_ratio;

  // Set initial gravitational acceleration
  G.y = 1./F;
  run();
}

scalar un[];

// Just a square box for now
event init (t=0) {

  /* Interfacial shape */
  fraction (f, -y + H0);

  /* // Place tracer in top half of bottom fluid
  fraction (T, intersection(0.25+y,0-y)); */

  /* // Place ball of tracer in middle of bottom fluid
  fraction (T, (pow(x,2)+pow((y+0.25),2)) < pow(0.1,2)); */

  // Place two balls of tracer on LHS and RHS of bottom fluid (actually just places second one, so LHS)
  fraction (T, (pow(x-0.19,2)+pow((y+0.25),2)) < pow(0.1,2));
  fraction (T, (pow(x+0.19,2)+pow((y+0.25),2)) < pow(0.1,2));

  // Set boundary conditions, no slip and no pen on all surfaces
  u.t[top] = dirichlet(0.);
  u.n[top] = dirichlet(0.);
  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);
  foreach()
    un[]=u.x[];

}

/* // Adapt the mesh to increase speed
event adapt (i++) {
  adapt_wavelet ({u.x}, (double []){4e-3}, maxlevel = 8);
} */

// Update "gravity" to produce our forcing
event update (i++; t<=t_end) {

  double theta_max = (1./3.)*PI/4.0;

  double theta = 0;

  double t_prime = 0;
  t_prime = fmod(t,1.0);
  // Time for one full oscillation (next period always starts rocking left but that means g points right)
  double t_max_prime = 1.;

  // Time reactor starts rocking back right
  double t_1_prime = t_max_prime/4.;

  // Time reactor starts rocking back left
  double t_2_prime = 3.*t_max_prime/4.;

  // Use piecewise function to define theta in terms of t_prime
  if(t_prime<=t_1_prime){
    theta = t_prime * theta_max / t_1_prime;
    } else if(t_prime<=t_2_prime){
      theta = theta_max - (4.*(t_prime - t_1_prime)/t_max_prime)*theta_max;
      } else {
        theta = -theta_max + (4.*(t_prime - t_2_prime)/t_max_prime)*theta_max;
      }
  G.x = (1./F) * sin(theta);
  G.y = -(1./F) * cos(theta);
  /* G.x = (1./F) * (1./sqrt(2.)) * sin(OMEGA * t);
  G.y = -1.0*(1./F) * ((1. - 1./sqrt(2)) * cos(2 * OMEGA * t) + (1/sqrt(2))); */
}

/* // Diffusion given peclet number of 1, also will diffuse into air?
event tracer_diffusion(i++) {
  diffusion (T, dt);
} */

// Output data during the calculation
event intermediate (t+=0.25) {
  heights(f,hei);
  int j, a, b;
  double u_x = 0;
  double u_y = 0;
  double y_step = 0.375;
  double x_step = 0;
  double x_pos = 0;
  double T_out = 0;
  double f_out = 0;
  char str[80];


  // Output data required to 
  // evaluate the mixing coefficient

  sum_sq = 0;
  sum_mean = 0;

  // Create the filename
  sprintf(str, "mix%.3f.dat", t);

  FILE * fp_3 = fopen(str,"w");

  for (b = 0; b<101; b++){
    for (a=0; a<101; a++) {
      T_out = interpolate(T,-0.5+a*0.01,-0.5+b*0.01);
      f_out = interpolate(f,-0.5+a*0.01,-0.5+b*0.01);
      u_x = interpolate(u.x,-0.5+a*0.01,-0.5+b*0.01);
      u_y = interpolate(u.y,-0.5+a*0.01,-0.5+b*0.01);
      fprintf(fp_3, "%lf %lf %lf %lf %lf %lf\n", -0.5+a*0.01, -0.5+b*0.01, T_out, f_out, u_x, u_y);
    }
  }
  fclose(fp_3);
  // Output the interface height

  // Create the filename
  /* sprintf(str, "interface%.3f.dat", t);

  FILE * fp_2 = fopen (str, "w");
  x_step = (double)1/(RESOLUTION);
  for (int j = 0; j < RESOLUTION; j++) {
    x_pos = -0.5 + j*x_step;
    fprintf(fp_2, "%lf %lf \n", x_pos, interfacial_height(x_pos));
  } // i end

  fclose(fp_2); */
}

/* // Output animation of tracer
event movies (t+=0.01) {
  char timestring[100];
	squares("T", map = cool_warm, min = 0, max = 1);
  draw_vof("f", lw=2);
  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Animations/Reactor_with_interface_new.mp4");  
} */

/* event profile (t=end)
{
  FILE * fp = fopen ("out.dat", "w");
  foreach() {
      fprintf (fp, "%g %g %g %g\n",
	       x, y, u.x[], u.y[]);
  }
  draw_vof ("cs","fs");
  squares ("u.x",linear = true, spread =-0.2);
  save ("u.x_tp.png");
} */

    

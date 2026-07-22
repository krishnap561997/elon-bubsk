#include "axi.h"
#include "navier-stokes/centered.h"

#include "two-phase.h"
#include "tension.h"

#include "navier-stokes/conserving.h"
//#include "reduced.h"
#include "view.h"
#include "embed.h"

vector h[];

//#include "hdf5_headers/output_xdmf_params.h"

double t_out = 0.01;
double t_dump = 0.1;
double t_end = 2.8;

double Rel = 92.8;        // 2*rho_l*U_l*R/mu_l
double Cal = 0.00464;     // mu_l*U_l/gamma
double Bo = 1.25;         // rho_l*g*R^2/gamma

double rho_ratio = 1./1000.,
       mu_ratio = 1./100.;

double bubL = 20.,
       LX = 50.;
double R = 1.,
       Ul = 1.;

double G;


int MAXlevel = 14,
    MINlwcwl = 7;

scalar f0[], profile[];


u.n[left]  = dirichlet(Ul*profile[]);
u.t[left]  = dirichlet(0.);
f[left] = dirichlet(0.);

/*u.n[right] = f[] > 1e-6 ? neumann(0.):dirichlet(0.) ;
u.t[right] = f[] > 1e-6 ? neumann(0.):dirichlet(0.) ;
p[right] = dirichlet(0);
pf[right] = dirichlet(0); */

p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);
u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
f[right] = neumann(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
f[embed] = dirichlet(0.);

/*u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
p[top] = dirichlet(0.);*/

void read_params(const char * fname);

int main (int argc, char * argv[])
{
  const char * fname = "params.in";
  if (argc > 1)
    fname = argv[1];

  TOLERANCE = 1e-4;
  NITERMIN = 2;
  NITERMAX = 100;
  CFL = 0.25;
//  DT = 5e-5;

  read_params(fname);

  size(LX);
  //dimensions(nx = AR, ny = 1);
  
  init_grid(1<<MINlevel);
  X0 = 0;
  Y0 = 0;

  
  rho1 = 1., rho2 =rho_ratio;
  mu1 = 2./Rel, mu2 = mu_ratio*mu1;

  f.sigma = 2./(Rel*Cal);
  f.height = h;

  //G.x = grav*sin(angle);
  //G.y = -grav*cos(angle);
  //Z.y = H0;
  G = -2.*B0/(Rel*Cal);

  char comm[80];
  sprintf(comm, "mkdir -p images");
  system(comm);
  
  sprintf(comm, "mkdir -p output");
  system(comm);

  sprintf(comm, "mkdir -p infc");
  system(comm);

  fprintf(stderr, "Re: %.8f\n", Rel);
  fprintf(stderr, "Ca: %.8f\n", Cal);
  fprintf(stderr, "Bo: %.8f\n", Bo);

  run();
}

void read_params(const char * fname)
{
  FILE * fp;
  if ((fp = fopen(fname, "rt"))) {
    char line[100];
    char key[80], val[80];

    while(fgets(line,100,fp)) {
      sscanf(line, "%15s = %15s", key, val);
      if (strcmp(key,"LX") == 0)              { LX         =  atof(val);        }
      else if (strcmp(key, "MAXLEVEL") == 0)  { MAXlevel   = atoi(val);         }
    }
    fclose(fp);
  } else {
    fprintf(stdout, "file %s not found\n", fname);
    //exit(0);
  }
}



event init (t = 0) {
  if (!restore (file = "dump")) { 
    solid(cs, fs, y-R);

    fraction (f0, H0 - y);
    restriction ({f0}); // for boundary conditions on levels

    foreach(){
      profile[] = 2.0*(1. - sq(y/R));
    }
    restriction ({profile}); // for boundary conditions on levels
   
   
    foreach() {
      //f[] = f0[];
      u.x[] = profile[]; // + 1-f0[]);
      u.y[] = 0;
    }
    boundary({d, u});
    
  }
}

event init_noise (i = 0)
{
  srand(noise_seed);

  for (int k = 0; k < MNOISE; k++)
    noise_phase[k] = 2.*pi*(rand()/((double) RAND_MAX));

  inlet_F = chang_noise_signal(t);
  inlet_U = US*(1. + inlet_F);
}

event check_grid(i=1)
{
  double xmax=0., ymax = 0., maxDelta = 0., minDelta = 10.;
  foreach(reduction(max:xmax) reduction(max:ymax) reduction(max:maxDelta) reduction(min:minDelta)){
    if(x > xmax) xmax = x;
    if(y > ymax) ymax = y;
    if(maxDelta < Delta) maxDelta = Delta;
    if(minDelta > Delta) minDelta = Delta;
  }

  fprintf(stderr, "N: %ld\n", grid->tn);
  fprintf(stderr, "Delta: %g , %g\n", maxDelta, minDelta);
  fprintf(stderr, "Domain: \nx : %g -> %g. \ny : %g -> %g\n", X0, xmax, Y0, ymax);
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(x){
    av.x[] += G;
  }
}

void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stdout, " \t - \t %d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
}


event logfile (i++) {
  if (i == 0)
    fprintf (stderr,
	     "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  fprintf (stderr, "%g %g %d %d %d %ld %g %g\n", 
	   t, dt, mgp.i, mgpf.i, mgu.i,
	   grid->tn, perf.t, perf.speed);
  fprintf (stdout, "\nPressure Residuals ");
  mg_print (mgp);
  fprintf (stdout, "\nVelocity Residuals ");
  mg_print (mgu);
  fprintf (stdout, "\n");
  fflush (stdout);
}


event interfacevel (t += t_out)
{
  char name[80];

  if (i==0)
  {
	clear();
        view (tx = -0.5, ty = -0.5, sx = zoomy, sy = 2*zoomy);
	draw_vof ("f", lw = 6);
	cells ();
	sprintf (name, "images/dimcheck-%5.4f.png", t);
	save (name);      
  }
  clear();
  view (tx = -0.5, ty = -0.5, sy = zoomy);
  draw_vof ("f", lw = 2);
  squares ("u.x", min = 0, max = 1.5*US, linear = true);
  colorbar(min = 0, max = 1.5*US);
  //isoline ("u.x", 1., lc = {1,1,1}, lw = 2);
  sprintf (name, "images/ux-%5.4f.png", t);
  save (name);

  clear();
  view (tx = -0.5, ty = -0.5, sy = 2*zoomy);
  draw_vof ("f", lw = 2);
  squares ("p", linear = true, spread=10);
  sprintf (name, "images/pfp-%5.4f.png", t);
  save (name);
}

/*event interface (t += t_out) {

   char names[80];
   sprintf(names, "infc/interface%d", pid());
   FILE * fp = fopen (names, "w");
   output_facets (f,fp);
   fclose(fp);
   char command[80];
   sprintf(command, "LC_ALL=C  cat infc/interfa* > infc/infc%05.4f.dat",t);
   system(command);
}*/

/* event velocityprofile (t += 0.0001) {

   char names[80];
   sprintf(names, "velocity%d", pid());
   FILE * fp = fopen (names, "w");
   foreach() {
     fprintf(fp, "%.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",
	     x, y, u.x[], u.y[], f[], d[], Delta);
   }

  fclose(fp);

   fclose(fp);
   char command[80];
   sprintf(command, "LC_ALL=C  cat velocity* > ux-vel%07.4f.dat",t);
   system(command);
} */

event savedump(t += t_dump; t <= t_end)
{
  char name[80];
  //scalar kappa[];
  //curvature(f, kappa);
  sprintf (name, "dump-%06.4f", t);
  p.nodump = false;
  dump (file = name); // so that we can restart
}

/*event output_h5(t += t_out; t<=t_end)
{
  char fname[256];
  scalar kappa[], kappa_minus[];

  sprintf(fname, "output/snapshot_%06.4f", t);

  output_xmf((scalar *){f,d,p,f_minus, d_minus, p_minus}, (vector *){u, u_minus}, fname);

  const char *parameter_names[] = {
    "Re",
    "Ca",
    "h0",
    "u0",
    "t0",
    "angle",
    "time",
    "time_minus",
    "dt"
  };

  double parameter_values[] = {
    RE,
    CA,
    H0,
    U0,
    T0,
    angle,
    t,
    t_minus,
    dt
  };

  int number_of_parameters =
    sizeof(parameter_values)/sizeof(parameter_values[0]);

  if (append_xmf_parameters (fname,
                             parameter_names,
                             parameter_values,
                             number_of_parameters) < 0) {
    if (pid() == 0)
      fprintf (stderr,
               "Could not append parameters to %s.h5\n",
               fname);
  }

  sprintf(fname, "infc/interface_%06.4f", t);
  foreach()
        kappa[] = distance_curvature (point, d);

  output_facets_xmf(f, kappa, fname);
  
  sprintf(fname, "infc/interface_minus_%06.4f", t);
  foreach()
        kappa_minus[] = distance_curvature (point, d_minus);

  output_facets_xmf(f_minus, kappa_minus, fname);
}*/


/*event save_previous (i++)
{
  t_minus = t;

  foreach() {
    f_minus[] = f[];
    p_minus[] = p[];
    d_minus[] = d[];
    foreach_dimension()
      u_minus.x[] = u.x[];
  }
}*/

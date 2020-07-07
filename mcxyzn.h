#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <omp.h>
#include <stdbool.h>

#define STRLEN 		32          /* String length. */
#define Ntiss		19          /* Number of tissue types. */

typedef struct mcxyz_config {
  int   mcflag, launchflag, boundaryflag, gradientflag;
  float	xfocus, yfocus, zfocus;
  float	ux0, uy0, uz0;
  float	radius;
  float	waist;
  float Nphotons;     /* number of photons in simulation */
  float	dx, dy, dz;     /* bin size [cm] */
  int  Nx, Ny, Nz, Nt; /* # of bins */
  long NN;
  float	xs, ys, zs;		/* launch position */
  char  myname[STRLEN];		// Holds the user's choice of myname, used in input and output files.
  char	filename[STRLEN];     // temporary filename for writing output.
  char  buf[32];                // buffer for reading header.dat
  /* time */
  float	time_min;               // Requested time duration of computation.
  time_t now;
  double start_time, finish_time, temp_time; /* for clock() */
  /* tissue parameters */
  char	tissuename[50][32];
  float muav[Ntiss];            // muav[0:Ntiss-1], absorption coefficient of ith tissue type
  float musv[Ntiss];            // scattering coeff.
  float gv[Ntiss];              // anisotropy of scattering
  float nv[Ntiss];
  long	i_photon;       /* current photon */
  FILE* fid; // file ID pointer
} mcconfig;

typedef struct mcxyz_float3 {
  float x,y,z;
} float3;

typedef struct mcxyz_float4 {
  float x,y,z,w;
} float4;

/* DECLARE FUNCTIONS */
void mcxyz_init(mcconfig *cfg,int argc, const char * argv[]);
void mcxyz_launchsimulation(mcconfig *cfg);
void mcxyz_kernel(mcconfig *cfg,unsigned char *v,float *F,int *dseed,float4 *g,const int Nphotons);

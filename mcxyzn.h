#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <omp.h>
#include <stdbool.h>
#include <omp.h>

#define STRLEN 		32          /* String length. */
#define Ntiss		19          /* Number of tissue types. */
#define ls          5E-6      /* Moving photon a little bit off the voxel face */
#define	PI          3.1415926
#define	LIGHTSPEED	2.997925E10 /* in vacuo speed of light [cm/s] */
#define ALIVE       1   		/* if photon not yet terminated */
#define DEAD        0    		/* if photon is to be terminated */
#define THRESHOLD   0.001		/* used in roulette */
#define CHANCE      10  		/* used in roulette */
#define Boolean     char
#define SQR(x)		pow(x,2)
#define SIGN(x)     ((x)>=0 ? 1:-1)
#define MAX_TIR     500

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
void mcxyzn_init(mcconfig *cfg,int argc, const char * argv[]);
void mcxyzn_launchsimulation(mcconfig *cfg);
void mcxyzn_kernel(mcconfig *cfg,unsigned char *v,float *F,int *dseed,float4 *g,const int Nphotons);
static inline unsigned long rotl(const unsigned long x, int k);
float RandomGen(unsigned long* s);
void LaunchPhoton(mcconfig* cfg, float4* pos, float4* u, float* rnd, unsigned long* seed);
int SameVoxel(mcconfig* cfg, float x1, float y1, float z1, float x2, float y2, float z2);
float FindVoxelFace2(mcconfig* cfg, float x1, float y1, float z1, float x2, float y2, float z2, float ux, float uy, float uz);
float RFresnel(float4* u, float4* g, float4* gb, float n1, float n2, unsigned long* seed, int* TIR_events, float* status);
int getindex(mcconfig* cfg, int x, int y, int z);
void InterpGradient(mcconfig* cfg, float4* g, unsigned char* v, float4* pos, float4* n, unsigned char tissue);
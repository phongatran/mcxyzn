/********************************************
 *  mcxyzn,	in ANSI Standard C programing language
 *      Usage:  mcxyz myname and myname_T.bin
 *      which loads myname_H.mci, and saves myname_F.bin.
 * 
 *  Initial version is an extension of mcxyz.c and the methodology is described in:
 *  A.P. Tran and S.L. Jacques, 2020. 
 *  Modeling voxel-based Monte Carlo light transport with curved and oblique boundary surfaces. 
 *  Journal of Biomedical Optics, 25(2), p.025001.
 *  
 *  USAGE   mcxyzn myname
 *              where myname is the user's choice.
 *          The program reads two files prepared by user:
 *                  myname_H.mci    = header input file for mcxyz
 *                  myname_T.bin    = tissue structure file
 *          The output will be written to 3 files:
 *                  myname_OP.m     = optical properties  (mua, mus, g for each tissue type)
 *                  myname_F.bin    = fluence rate output F[i] [W/cm^2 per W delivered]
 *
 *  The MATLAB program maketissue.m can create the two input files (myname_H.mci, myname_T.bin).
 *
 *  The MATLAB program lookmcxyzn.m can read the output files and display
 *          1. Fluence rate F [W/cm^2 per W delivered]
 *          2. Deposition rate A [W/cm^3 per W delivered].
 *
 *  Log:
 *  Original mcxyz.c was created by Steven L. Jacques and Ting Li (Oregon Health & Science University), 2010/2012.
 *  Written by Ting based on Steve's mcsub.c., 2010.
 *      Use Ting's FindVoxelFace().
 *	Use Steve's FindVoxelFace(), Dec. 30, 2010.
 *  Reorganized by Steve. May 8, 2012:
 *      Reads input files, outputs binary files.
 **********/

#include "mcxyzn.h"

int main(int argc, const char * argv[]) {
    if (argc==0) {
        printf("assuming you've compiled mcxyz.c as gomcxyz ...\n");
        printf("USAGE: gomcxyz name\n");
        printf("which will load the files name_H.mci and name_T.bin\n");
        printf("and run the Monte Carlo program.\n");
        printf("Yields  name_F.bin, which holds the fluence rate distribution.\n");
        return 0;
    }
    mcconfig cfg;

    mcxyzn_init(&cfg,argc,argv);
    mcxyzn_launchsimulation(&cfg);
    return 0;
} /* end of main */

void mcxyzn_init(mcconfig *cfg,int argc, const char * argv[])
{
    long int i;
    /* Input/Output */
    strcpy(cfg->myname, argv[1]);    // acquire name from argument of function call by user.
    printf("name = %s\n",cfg->myname);

    /**** INPUT FILES *****/
    /* IMPORT myname_H.mci */
    strcpy(cfg->filename,cfg->myname);
    strcat(cfg->filename, "_H.mci");
    cfg->fid = fopen(cfg->filename,"r");
    fgets(cfg->buf, 32, cfg->fid);
    // run parameters
    sscanf(cfg->buf, "%f", &cfg->time_min); // desired time duration of run [min]
    fgets(cfg->buf, 32, cfg->fid);
    sscanf(cfg->buf, "%d", &cfg->Nx);  // # of bins
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%d", &cfg->Ny);  // # of bins
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%d", &cfg->Nz);  // # of bins

    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->dx);	 // size of bins [cm]
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->dy);	 // size of bins [cm]
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->dz);	 // size of bins [cm]

    // launch parameters
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%d", &cfg->mcflag);  // mcflag, 0 = uniform, 1 = Gaussian, 2 = iso-pt
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%d", &cfg->launchflag);  // launchflag, 0 = ignore, 1 = manually set
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%d", &cfg->boundaryflag);  // 0 = no boundaries, 1 = escape at all boundaries, 2 = escape at surface only
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%d", &cfg->gradientflag);
  
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->xs);  // initial launch point
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->ys);  // initial launch point
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->zs);  // initial launch point

    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->xfocus);  // xfocus
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->yfocus);  // yfocus
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->zfocus);  // zfocus

    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->ux0);  // ux trajectory
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->uy0);  // uy trajectory
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->uz0);  // uz trajectory

    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->radius);  // radius
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%f", &cfg->waist);  // waist


    // tissue optical properties
    fgets(cfg->buf, 32,cfg->fid);
    sscanf(cfg->buf, "%d", &cfg->Nt);				// # of tissue types in tissue list
    for (i=1; i<=cfg->Nt; i++) {
        fgets(cfg->buf, 32, cfg->fid);
        sscanf(cfg->buf, "%f", &cfg->muav[i]);	// absorption coeff [cm^-1]
        fgets(cfg->buf, 32, cfg->fid);
        sscanf(cfg->buf, "%f", &cfg->musv[i]);	// scattering coeff [cm^-1]
        fgets(cfg->buf, 32, cfg->fid);
        sscanf(cfg->buf, "%f", &cfg->gv[i]);		// anisotropy of scatter [dimensionless]
        fgets(cfg->buf, 32, cfg->fid);
        sscanf(cfg->buf, "%f", &cfg->nv[i]);
    }
    fclose(cfg->fid);

    printf("time_min = %0.2f min\n",cfg->time_min);
    printf("Nx = %d, dx = %0.4f [cm]\n",cfg->Nx,cfg->dx);
    printf("Ny = %d, dy = %0.4f [cm]\n",cfg->Ny,cfg->dy);
    printf("Nz = %d, dz = %0.4f [cm]\n",cfg->Nz,cfg->dz);

    printf("xs = %0.4f [cm]\n",cfg->xs);
    printf("ys = %0.4f [cm]\n",cfg->ys);
    printf("zs = %0.4f [cm]\n",cfg->zs);
    printf("mcflag = %d\n",cfg->mcflag);
    if (cfg->mcflag==0) printf("launching uniform flat-field beam\n");
    if (cfg->mcflag==1) printf("launching Gaussian beam\n");
    if (cfg->mcflag==2) printf("launching isotropic point source\n");
    if (cfg->mcflag==3) printf("launching square source\n");
    printf("xfocus = %0.4f [cm]\n",cfg->xfocus);
    printf("yfocus = %0.4f [cm]\n",cfg->yfocus);
    printf("zfocus = %0.2e [cm]\n",cfg->zfocus);
    if (cfg->launchflag==1) {
        printf("Launchflag ON, so launch the following:\n");
        printf("ux0 = %0.4f [cm]\n",cfg->ux0);
        printf("uy0 = %0.4f [cm]\n",cfg->uy0);
        printf("uz0 = %0.4f [cm]\n",cfg->uz0);
    }
    else {
        printf("Launchflag OFF, so program calculates launch angles.\n");
        printf("radius = %0.4f [cm]\n",cfg->radius);
        printf("waist  = %0.4f [cm]\n",cfg->waist);
    }
    if (cfg->boundaryflag==0)
        printf("boundaryflag = 0, so no boundaries.\n");
    else if (cfg->boundaryflag==1)
        printf("boundaryflag = 1, so escape at all boundaries.\n");
    else if (cfg->boundaryflag==2)
        printf("boundaryflag = 2, so escape at surface only.\n");
    else{
        printf("improper boundaryflag. quit.\n");
        //return 0;
    }
    printf("# of tissues available, Nt = %d\n",cfg->Nt);
    for (i=1; i<=cfg->Nt; i++) {
        printf("muav[%ld] = %0.4f [cm^-1]\n",i,cfg->muav[i]);
        printf("musv[%ld] = %0.4f [cm^-1]\n",i,cfg->musv[i]);
        printf("  gv[%ld] = %0.4f [--]\n",i,cfg->gv[i]);
        printf("  nv[%ld] = %0.4f [--]\n\n",i,cfg->nv[i]);
    }

    // SAVE optical properties, for later use by MATLAB.
    strcpy(cfg->filename,cfg->myname);
    strcat(cfg->filename,"_props.m");
    cfg->fid = fopen(cfg->filename,"w");
    for (i=1; i<=cfg->Nt; i++) {
        fprintf(cfg->fid,"muav(%ld) = %0.4f;\n",i,cfg->muav[i]);
        fprintf(cfg->fid,"musv(%ld) = %0.4f;\n",i,cfg->musv[i]);
        fprintf(cfg->fid,"gv(%ld) = %0.4f;\n",i,cfg->gv[i]);
        fprintf(cfg->fid,"nv(%ld) = %0.4f;\n\n",i,cfg->nv[i]);
    }
    fclose(cfg->fid);

    /* IMPORT BINARY TISSUE FILE */
    cfg->NN = cfg->Nx*cfg->Ny*cfg->Nz;

    return;
}

void mcxyzn_launchsimulation(mcconfig *cfg)
{
    unsigned char *v;
    float *F;
    float4 *g;
    //float *R_host;
    
    /* Initializing storing elements */
    v  = (unsigned char *)malloc(cfg->NN*sizeof(unsigned char));  /* tissue structure */
    F  = (float *)malloc(cfg->NN*sizeof(float));	/* relative fluence rate [W/cm^2/W.delivered] */
    g = (float4 *)malloc(cfg->NN*sizeof(float4));
    //R  = (float *)malloc(cfg->Nx*cfg->Ny*sizeof(cl_mem));	/* escaping flux [W/cm^2/W.delivered] */
    //for (i=0; i<Ny*Nx; i++) R[i] = 0;

    //prop_host = (float4 *)malloc(cfg->NN*sizeof(float4));

    //cl_float *R  = (cl_float *)(R_host);
    /* read binary file */
    strcpy(cfg->filename,cfg->myname);
    strcat(cfg->filename, "_T.bin");
    cfg->fid = fopen(cfg->filename, "rb");
    fread(v, sizeof(unsigned char), cfg->NN, cfg->fid);
    fclose(cfg->fid);

    if (cfg->gradientflag > 0){
        strcpy(cfg->filename,cfg->myname);
        strcat(cfg->filename, "_Gx.bin");
        cfg->fid = fopen(cfg->filename, "rb");
        fread(F, sizeof(float), cfg->NN, cfg->fid);
        fclose(cfg->fid);
        for(int j=0; j<cfg->NN;j++)
        {
            g[j].x = F[j]; // ensure F[] starts empty.
        }
        strcpy(cfg->filename,cfg->myname);
        strcat(cfg->filename, "_Gy.bin");
        cfg->fid = fopen(cfg->filename, "rb");
        fread(F, sizeof(float), cfg->NN, cfg->fid);
        fclose(cfg->fid);
        for(int j=0; j<cfg->NN;j++)
        {
            g[j].y = F[j]; // ensure F[] starts empty.
        }
        strcpy(cfg->filename,cfg->myname);
        strcat(cfg->filename, "_Gz.bin");
        cfg->fid = fopen(cfg->filename, "rb");
        fread(F, sizeof(float), cfg->NN, cfg->fid);
        fclose(cfg->fid);
        for(int j=0; j<cfg->NN;j++)
        {
            g[j].z = F[j]; // ensure F[] starts empty.
        }
    }
    
    for(int j=0; j<cfg->NN;j++)
    {
        F[j] = 0.f; // ensure F[] starts empty.
    }
    /* Show tissue on screen, along central z-axis, by listing tissue type #'s.*/
    printf("central axial profile of tissue types:\n");
    for (int iz=0; iz<cfg->Nz; iz++) {
        int i = (long)(iz*cfg->Ny*cfg->Nx + (cfg->Ny/2)*cfg->Nx + cfg->Nx/2);
        printf("%d",v[i]);
    }
    printf("\n\n");

    /*************************************
     * == Setting up OpenCL structure == *
     *************************************/
    /* Number of photons launched */
    cfg->Nphotons = 50000;
    int nb_threads = omp_get_max_threads();
    /* Create seeds for threads */
    int *seed;
    seed = (int *)malloc(nb_threads*sizeof(int)*2);
    for (int i=0;i<nb_threads;i++){
      seed[i*2]= rand();
      seed[i*2+1] = rand();
    }
    printf("Creating random seed of length (2 seeds per thread): %i \n",(int)nb_threads*2);
    

    /*******************************************************************
     * ============================ MAJOR CYCLE ========================
     *******************************************************************/
    cfg->start_time = clock();
    cfg->now = time(NULL);
    printf("\n%s\n", ctime(&cfg->now));
    /* Launch main kernel */

    printf("[====== Main kernel ======]\n");
    printf("Launching %i photons with %i threads. \n",(int)cfg->Nphotons,(int)nb_threads);
    
    double start,start2;
    double end,end2;
    start = omp_get_wtime();
    mcxyz_kernel(cfg,v,F,seed,g,cfg->Nphotons);
    end = omp_get_wtime();
    printf("Test kernel of %i photons took %f sec.\n",(int) cfg->Nphotons, end-start);
    cfg->Nphotons = ceil(cfg->Nphotons/(end-start)*cfg->time_min*60)- cfg->Nphotons;
    printf("Launching remaining %i photons.\n", (int) cfg->Nphotons);
    start2 = omp_get_wtime();
    mcxyz_kernel(cfg,v,F,seed,g,cfg->Nphotons);
    end2 = omp_get_wtime();
    printf("Main kernel took %f sec.\n",end2-start2);
    printf("Total running time of %f sec for %i photons. \n",end+end2-start-start2,(int)cfg->Nphotons+50000);

    /* printf("------------------------------------------------------\n"); */
    /* printf("Elapsed Time for %i photons = %f sec\n",(int)cfg->Nphotons,(float)(end-start)/(1e9)); */
    /* printf("%i photons per minute\n", (int) (cfg->Nphotons/(end-start)*(1e9)*60)); */
    /* printf("------------------------------------------------------\n"); */

    /**************
     * == Save == *
     **************/

    // Normalize deposition (A) to yield fluence rate (F).
    float temp = cfg->dx*cfg->dy*cfg->dz*(cfg->Nphotons);
    for (int i=0; i<cfg->NN;i++){
      F[i] = (F[i]/(temp*cfg->muav[v[i]]));
    }
    // Save the binary file
    strcpy(cfg->filename,cfg->myname);
    strcat(cfg->filename,"_F.bin");
    printf("saving %s\n",cfg->filename);
    cfg->fid = fopen(cfg->filename, "wb");   /* 3D voxel output */
    fwrite(F, sizeof(float), cfg->NN, cfg->fid);
    fclose(cfg->fid);
    
    /* save reflectance */
    /*float temp = cfg->dx*cfg->dy*(cfg->Nphotons+50000);
    for (int i=0; i<cfg->Nx*cfg->Ny;i++){
      R[i] = (F[i]/(temp));
    }
    strcpy(filename,myname);
    strcat(filename,"_Ryx.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 2D voxel output */
    /*fwrite(R, sizeof(float), cfg->Ny*cfg->Nx, fid);
    fclose(fid);
    printf("%s is done.\n",myname);*/

    printf("------------------------------------------------------\n");
    cfg->now = time(NULL);
    printf("%s\n", ctime(&cfg->now));

    free(F);
    free(v);
    free(g);
    //free(R_host);

    return;
}


/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */
/* SUBROUTINES */

static inline unsigned long rotl(const unsigned long x, int k) {
    return (x << k) | (x >> (64 - k));
}

/*********************
 *	RandomGen    *
 *********************/
float RandomGen(unsigned long* s) {
    union {
        unsigned long i;
        unsigned int  u[2];
        float f[2];
    } result;

    result.i = s[0] + s[1];
    s[1] ^= s[0];
    s[0] = rotl(s[0], 24) ^ s[1] ^ (s[1] << 16); // a, b
    s[1] = rotl(s[1], 37); // c
    result.u[0] = 0x3F800000U | (result.u[0] >> 9);

    return result.f[0] - 1.f;
}


/************* SET SOURCE***************
 * Launch collimated beam at x,y center.
 ***************************************/
void LaunchPhoton(mcconfig* cfg, float4* pos, float4* u, float* rnd, unsigned long* seed)
{
    float r, phi, temp;
    /****************************/
    /* Initial position. */
    /* trajectory */
    if (cfg->launchflag == 1) { // manually set launch
        pos->x = cfg->xs;
        pos->y = cfg->ys;
        pos->z = cfg->zs;
        u->x = cfg->ux0;
        u->y = cfg->uy0;
        u->z = cfg->uz0;
    }
    else { // use mcflag
        if (cfg->mcflag == 0) { // uniform beam
            //set launch point and width of beam
            while ((*rnd = RandomGen(seed)) <= 0.0); // avoids rnd = 0
            r = cfg->radius * sqrt(*rnd); // radius of beam at launch point
            while ((*rnd = RandomGen(seed)) <= 0.0); // avoids rnd = 0
            phi = (*rnd) * 2.0 * PI;
            pos->x = cfg->xs + r * cos(phi);
            pos->y = cfg->ys + r * sin(phi);
            pos->z = cfg->zs;
            // set trajectory toward focus
            while ((*rnd = RandomGen(seed)) <= 0.0); // avoids rnd = 0
            r = cfg->waist * sqrt(*rnd); // radius of beam at focus
            while ((*rnd = RandomGen(seed)) <= 0.0); // avoids rnd = 0
            phi = (*rnd) * 2.0 * PI;
            float xfocus = cfg->xs + r * cos(phi); //SLJ add cfg->xs
            float yfocus = cfg->ys + r * sin(phi); //SLJ add cfg->ys
            temp = 1 / sqrt((pos->x - xfocus) * (pos->x - xfocus) + (pos->y - yfocus) * (pos->y - yfocus) + (pos->z - cfg->zfocus) * (pos->z - cfg->zfocus));
            u->x = -(pos->x - xfocus) * temp;
            u->y = -(pos->y - yfocus) * temp;
            u->z = sqrt(1 - u->x * u->x - u->y * u->y);
        }
        else if (cfg->mcflag == 2) { // isotropic pt source
            float ctheta = 1.0 - 2.0 * RandomGen(seed);
            float stheta = sqrt(1.0 - ctheta * ctheta);
            float psi = 2.0 * PI * RandomGen(seed);
            float cpsi = cos(psi);
            float spsi;
            if (psi < PI)
                spsi = sqrt(1.0 - cpsi * cpsi);
            else
                spsi = -sqrt(1.0 - cpsi * cpsi);
            pos->x = cfg->xs;
            pos->y = cfg->ys;
            pos->z = cfg->zs;
            u->x = stheta * cpsi;
            u->y = stheta * spsi;
            u->z = ctheta;
        }
        else if (cfg->mcflag == 3) { // rectangular source collimated
            while ((*rnd = RandomGen(seed)) <= 0.0); // avoids rnd = 0
            pos->x = cfg->radius * ((*rnd) * 2 - 1); // use radius to specify x-halfwidth of rectangle
            while ((*rnd = RandomGen(seed)) <= 0.0); // avoids rnd = 0
            pos->y = cfg->radius * ((*rnd) * 2 - 1); // use radius to specify y-halfwidth of rectangle
            pos->z = cfg->zs;
            u->x = 0.0;
            u->y = 0.0;
            u->z = 1.0; // collimated beam
        }
    } // end  use mcflag
    pos->x = cfg->Nx / 2 + pos->x / cfg->dx;
    pos->y = cfg->Ny / 2 + pos->y / cfg->dy;
    pos->z = pos->z / cfg->dz;
    /****************************/
}


/***********************************************************
 *  Determine if the two position are located in the same voxel
 *	Returns 1 if same voxel, 0 if not same voxel.
 ****/
int SameVoxel(mcconfig* cfg, float x1, float y1, float z1, float x2, float y2, float z2)
{
    float xmin = fmin((floor)(x1), (floor)(x2));
    float ymin = fmin((floor)(y1), (floor)(y2));
    float zmin = fmin((floor)(z1), (floor)(z2));
    float xmax = xmin + 1;
    float ymax = ymin + 1;
    float zmax = zmin + 1;
    return ((x1 <= xmax && x2 <= xmax && y1 <= ymax && y2 <= ymax && z1 < zmax && z2 <= zmax));
}

/********************
 * my version of FindVoxelFace for no scattering.
 * s = ls + FindVoxelFace2(x,y,z, tempx, tempy, tempz, dx, dy, dz, ux, uy, uz);
 ****/
float FindVoxelFace2(mcconfig* cfg, float x1, float y1, float z1, float x2, float y2, float z2, float ux, float uy, float uz)
{
    int ix1 = floor(x1);
    int iy1 = floor(y1);
    int iz1 = floor(z1);

    int ix2, iy2, iz2;
    if (ux >= 0)
        ix2 = ix1 + 1;
    else
        ix2 = ix1;

    if (uy >= 0)
        iy2 = iy1 + 1;
    else
        iy2 = iy1;

    if (uz >= 0)
        iz2 = iz1 + 1;
    else
        iz2 = iz1;

    float xs = fabs((ix2 - x1) / ux);
    float ys = fabs((iy2 - y1) / uy);
    float zs = fabs((iz2 - z1) / uz);

    float s = fmin(xs, fmin(ys, zs));

    return (s * cfg->dx);
}

/***********************************************************
 *	FRESNEL REFLECTANCE
 * Computes reflectance as photon passes from medium 1 to
 * medium 2 with refractive indices n1,n2. Incident
 * angle a1 is specified by cosine value ca1 = cos(a1).
 * Program returns value of transmitted angle a1 as
 * value in *ca2_Ptr = cos(a2).
 ****/
float RFresnel(float4* u, float4* g, float4* gb, float n1, float n2, unsigned long* seed, int* TIR_events, float* status)
{
    if (n1 == n2) {
        return 1.0;
    }
    else {
        if ((g->x * g->x + g->y * g->y + g->z * g->z) == 0) {
            *g = *gb;
        }
        float rand = RandomGen(seed);
        float cos_i = -(u->x) * g->x - (u->y) * g->y - (u->z) * g->z;
        if (cos_i > 0.99999) {
            float r = (n2 - n1) / (n2 + n1);
            r *= r;
            if (rand > r) {
                //u->x = -g->x, u->y = -g->y, u->z = -g->z;
                return 1.0;
            }
            else {
                u->x = -u->x, u->y = -u->y, u->z = -u->z;
                return 0.0;
            }
        }
        else if (cos_i < 1e-5) {
            u->x = u->x + 2 * cos_i * g->x;
            u->y = u->y + 2 * cos_i * g->y;
            u->z = u->z + 2 * cos_i * g->z;
            return 0.0;
        }
        else {
            float sin_t2 = pow(n1 / n2, 2) * (1 - cos_i * cos_i);
            if (sin_t2 >= 1.0) {
                if (*TIR_events < MAX_TIR) {
                    u->x = u->x + 2 * cos_i * g->x;
                    u->y = u->y + 2 * cos_i * g->y;
                    u->z = u->z + 2 * cos_i * g->z;
                    (*TIR_events)++;
                    return 0.0;
                }
                else {
                    *status = DEAD;
                    u->w = 0.0;
                    return 1.0;
                }
            }
            else {
                float cos_t = sqrt(1.0 - sin_t2);
                float temp1 = n1 * cos_i;
                float temp2 = n2 * cos_t;
                temp1 = (temp1 - temp2) / (temp1 + temp2);
                float r = 0.5 * temp1 * temp1;
                temp1 = n2 * cos_i;
                temp2 = n1 * cos_t;
                temp1 = (temp1 - temp2) / (temp1 + temp2);
                r += 0.5 * temp1 * temp1;
                if (rand > r) {
                    temp1 = n1 / n2;
                    temp2 = temp1 * cos_i - cos_t;
                    u->x = temp1 * (u->x) + temp2 * g->x;
                    u->y = temp1 * (u->y) + temp2 * g->y;
                    u->z = temp1 * (u->z) + temp2 * g->z;
                    return 1.0;
                }
                else {
                    u->x = u->x + 2 * cos_i * g->x;
                    u->y = u->y + 2 * cos_i * g->y;
                    u->z = u->z + 2 * cos_i * g->z;
                    return 0.0;
                }
            }
        }
    }
} /******** END SUBROUTINE **********/

int getindex(mcconfig* cfg, int x, int y, int z)
{
    return z * cfg->Ny * cfg->Nx + x * cfg->Ny + y;
}

void InterpGradient(mcconfig* cfg, float4* g, unsigned char* v, float4* pos, float4* n, unsigned char tissue)
{
    if (pos->x >= cfg->Nx - 0.5) { pos->x = cfg->Nx - 0.51; }
    if (pos->y >= cfg->Ny - 0.5) { pos->y = cfg->Ny - 0.51; }
    if (pos->z >= cfg->Nz - 0.5) { pos->z = cfg->Nz - 0.51; }
    if (pos->x < 0.5) { pos->x = 0.51; }
    if (pos->y < 0.5) { pos->y = 0.51; }
    if (pos->z < 0.5) { pos->z = 0.51; }

    float x = round(pos->x);
    float y = round(pos->y);
    float z = round(pos->z);

    float xd = pos->x - x + 0.5;
    float yd = pos->y - y + 0.5;
    float zd = pos->z - z + 0.5;

    float v000, v001, v010, v011, v100, v101, v110, v111;
    v000 = (v[getindex(cfg, x - 1, y - 1, z - 1)] == tissue);
    v001 = (v[getindex(cfg, x - 1, y - 1, z)] == tissue);
    v010 = (v[getindex(cfg, x - 1, y, z - 1)] == tissue);
    v011 = (v[getindex(cfg, x - 1, y, z)] == tissue);
    v100 = (v[getindex(cfg, x, y - 1, z - 1)] == tissue);
    v101 = (v[getindex(cfg, x, y - 1, z)] == tissue);
    v110 = (v[getindex(cfg, x, y, z - 1)] == tissue);
    v111 = (v[getindex(cfg, x, y, z)] == tissue);

    float c00 = (1 - xd) * g[getindex(cfg, x - 1, y - 1, z - 1)].x * v000 + xd * g[getindex(cfg, x, y - 1, z - 1)].x * v100;
    float c01 = (1 - xd) * g[getindex(cfg, x - 1, y - 1, z)].x * v001 + xd * g[getindex(cfg, x, y - 1, z)].x * v101;
    float c10 = (1 - xd) * g[getindex(cfg, x - 1, y, z - 1)].x * v010 + xd * g[getindex(cfg, x, y, z - 1)].x * v110;
    float c11 = (1 - xd) * g[getindex(cfg, x - 1, y, z)].x * v011 + xd * g[getindex(cfg, x, y, z)].x * v111;
    float c0 = (1 - yd) * c00 + yd * c10;
    float c1 = (1 - yd) * c01 + yd * c11;
    n->x = c0 * (1 - zd) + c1 * zd;

    c00 = (1 - xd) * g[getindex(cfg, x - 1, y - 1, z - 1)].y * v000 + xd * g[getindex(cfg, x, y - 1, z - 1)].y * v100;
    c01 = (1 - xd) * g[getindex(cfg, x - 1, y - 1, z)].y * v001 + xd * g[getindex(cfg, x, y - 1, z)].y * v101;
    c10 = (1 - xd) * g[getindex(cfg, x - 1, y, z - 1)].y * v010 + xd * g[getindex(cfg, x, y, z - 1)].y * v110;
    c11 = (1 - xd) * g[getindex(cfg, x - 1, y, z)].y * v011 + xd * g[getindex(cfg, x, y, z)].y * v111;
    c0 = (1 - yd) * c00 + yd * c10;
    c1 = (1 - yd) * c01 + yd * c11;
    n->y = c0 * (1 - zd) + c1 * zd;

    c00 = (1 - xd) * g[getindex(cfg, x - 1, y - 1, z - 1)].z * v000 + xd * g[getindex(cfg, x, y - 1, z - 1)].z * v100;
    c01 = (1 - xd) * g[getindex(cfg, x - 1, y - 1, z)].z * v001 + xd * g[getindex(cfg, x, y - 1, z)].z * v101;
    c10 = (1 - xd) * g[getindex(cfg, x - 1, y, z - 1)].z * v010 + xd * g[getindex(cfg, x, y, z - 1)].z * v110;
    c11 = (1 - xd) * g[getindex(cfg, x - 1, y, z)].z * v011 + xd * g[getindex(cfg, x, y, z)].z * v111;
    c0 = (1 - yd) * c00 + yd * c10;
    c1 = (1 - yd) * c01 + yd * c11;
    n->z = c0 * (1 - zd) + c1 * zd;

    float magn = sqrt(n->x * n->x + n->y * n->y + n->z * n->z);
    n->x = n->x / magn;
    n->y = n->y / magn;
    n->z = n->z / magn;
    return;
}

void mcxyzn_kernel(mcconfig* cfg, unsigned char* v, float* F, int* dseed, float4* g, const int Nphotons)
{
    #pragma omp parallel
    {
        int idx = omp_get_thread_num();
        unsigned long seed[2];
        seed[0] = dseed[idx * 2];
        seed[1] = dseed[idx * 2 + 1];

        #pragma omp for
        for (int k = 0; k < Nphotons; k++)
        {
            //if (idx == 0)
            //printf("Thread %i photon %i/%i Rand %f %i %i\n",idx,k+1,Nphotons,RandomGen(seed),seed[0],seed[1]);

            /**** LAUNCH
              Initialize photon position and trajectory.
            *****/
            //if (fmod(i_photon,10)==0) printf("photon %ld took %d steps\n",i_photon,CNT);
            float  rnd;            /* assigned random value 0-1 */			/* dummy values */
            float  temp;           /* dummy variable */
            int  ix, iy, iz;     /* Added. Used to track photons */
            float4 pos;            /* photon position .w = weight */
            float4 u;              /* photon trajectory .w = sleft */
            unsigned char  type;              /* absorption coef [cm^-1] scattering coef [cm^-1] anisotropy [-] refractivity index [-] */
            pos.w = 1.0;                    /* set photon weight to one */
            float status = ALIVE;      /* Launch an ALIVE photon */
            int TIR_events = 0;
            LaunchPhoton(cfg, &pos, &u, &rnd, seed);

            /* Get tissue voxel properties of launchpoint.
             * If photon beyond outer edge of defined voxels,
             * the tissue equals properties of outermost voxels.
             * Therefore, set outermost voxels to infinite background value.
             */
            ix = (int)(pos.x);
            iy = (int)(pos.y);
            iz = (int)(pos.z);
            if (ix >= cfg->Nx) ix = cfg->Nx - 1;
            if (iy >= cfg->Ny) iy = cfg->Ny - 1;
            if (iz >= cfg->Nz) iz = cfg->Nz - 1;
            if (ix < 0)   ix = 0;
            if (iy < 0)   iy = 0;
            if (iz < 0)   iz = 0;
            /* Get the tissue type of located voxel */
            int i = (int)(iz * cfg->Ny * cfg->Nx + ix * cfg->Ny + iy);//(iz*cfg->Ny*cfg->Nx + ix*cfg->Ny + iy);
            type = v[i];
            int bflag = 1;

            /* HOP_DROP_SPIN_CHECK
               Propagate one photon until it dies as determined by ROULETTE.
            *******/
            do {
                /**** HOP
                  Take step to new position
                  s = dimensionless stepsize
                  x, uy, uz are cosines of current photon trajectory
                *****/
                while ((rnd = RandomGen(seed)) <= ls);   /* yields 0 < rnd <= 1 */
                u.w = -log(rnd);				/* dimensionless step */

                do {  // while sleft>0
                    float s = u.w / cfg->musv[type];			/* Step size [cm].*/
                    float tempx = pos.x + s * u.x / cfg->dx;				/* Update positions. [cm] */
                    float tempy = pos.y + s * u.y / cfg->dx;
                    float tempz = pos.z + s * u.z / cfg->dx;
                    if (SameVoxel(cfg, pos.x, pos.y, pos.z, tempx, tempy, tempz)) /* photon in same voxel */
                    {
                        pos.x = tempx;					/* Update positions. */
                        pos.y = tempy;
                        pos.z = tempz;

                        /**** DROP
                          Drop photon weight (W) into local bin.
                        *****/
                        float absorb = pos.w * (1 - exp(-cfg->muav[type] * s));	/* photon weight absorbed at this step */
                        if (absorb != absorb) {
                            status = DEAD;
                            u.w = 0;
                        }
                        else {
                            //atomicadd(&(F[i]),absorb);
                            if (bflag) {
                                #pragma omp atomic
                                F[i] += absorb;
                            }
                            pos.w -= absorb;					/* decrement WEIGHT by amount absorbed */
                        }
                        // If photon within volume of heterogeneity, deposit energy in F[].
                        // Normalize F[] later, when save output.


                        /* Update sleft */
                        u.w = 0;		/* dimensionless step remaining */
                    }
                    else /* photon has crossed voxel boundary */
                    {
                        /* step to voxel face + "littlest step" so just inside new voxel. */
                        s = ls + FindVoxelFace2(cfg, pos.x, pos.y, pos.z, tempx, tempy, tempz, u.x, u.y, u.z);

                        float temp_px, temp_py, temp_pz;
                        /* Update positions. */
                        temp_px = pos.x + s * u.x / cfg->dx;
                        temp_py = pos.y + s * u.y / cfg->dx;
                        temp_pz = pos.z + s * u.z / cfg->dx;

                        /**** DROP
                          Drop photon weight (W) into local bin.
                        *****/
                        float absorb = pos.w * (1 - exp(-cfg->muav[type] * s));   /* photon weight absorbed at this step */
                        if (absorb != absorb) {
                            status = DEAD;
                            u.w = 0;
                        }
                        else {
                            //atomicadd(&(F[i]),absorb);
                            if (bflag) {
                                #pragma omp atomic
                                F[i] += absorb;
                            }
                            pos.w -= absorb;					/* decrement WEIGHT by amount absorbed */
                        }

                        /* Update sleft */
                        u.w -= s * cfg->musv[type];  /* dimensionless step remaining */
                        if (u.w <= ls) u.w = 0;

                        int temp_ix = (int)floor(temp_px);
                        int temp_iy = (int)floor(temp_py);
                        int temp_iz = (int)floor(temp_pz);
                        bflag = 1; //boundary flag. Initalize as 1 = inside volume, then check;
                        if (cfg->boundaryflag == 0) {
                            if (temp_iz >= cfg->Nz) { temp_iz = cfg->Nz - 1; bflag = 0; }
                            if (temp_ix >= cfg->Nx) { temp_ix = cfg->Nx - 1; bflag = 0; }
                            if (temp_iy >= cfg->Ny) { temp_iy = cfg->Ny - 1; bflag = 0; }
                            if (temp_iz < 0) { temp_iz = 0; bflag = 0; }
                            if (temp_ix < 0) { temp_ix = 0; bflag = 0; }
                            if (temp_iy < 0) { temp_iy = 0; bflag = 0; }
                        }
                        else if (cfg->boundaryflag == 1) {
                            if (temp_iz >= cfg->Nz) { temp_iz = cfg->Nz - 1; status = DEAD; u.w = 0; }
                            if (temp_ix >= cfg->Nx) { temp_ix = cfg->Nx - 1; status = DEAD; u.w = 0; }
                            if (temp_iy >= cfg->Ny) { temp_iy = cfg->Ny - 1; status = DEAD; u.w = 0; }
                            if (temp_iz < 0) { temp_iz = 0; status = DEAD; u.w = 0; }
                            if (temp_ix < 0) { temp_ix = 0; status = DEAD; u.w = 0; }
                            if (temp_iy < 0) { temp_iy = 0; status = DEAD; u.w = 0; }
                        }
                        else if (cfg->boundaryflag == 2) {
                            if (temp_iz >= cfg->Nz) { temp_iz = cfg->Nz - 1; bflag = 0; }
                            if (temp_ix >= cfg->Nx) { temp_ix = cfg->Nx - 1; bflag = 0; }
                            if (temp_iy >= cfg->Ny) { temp_iy = cfg->Ny - 1; bflag = 0; }
                            if (temp_iz < 0) { temp_iz = 0; status = DEAD; u.w = 0; }
                            if (temp_ix < 0) { temp_ix = 0; bflag = 0; }
                            if (temp_iy < 0) { temp_iy = 0; bflag = 0; }
                        }

                        int p = (int)(temp_iz * cfg->Ny * cfg->Nx + temp_ix * cfg->Ny + temp_iy);

                        int fstatus = 1;
                        float4 ga = g[i];
                        float4 gb = { ix - temp_ix,iy - temp_iy,iz - temp_iz,0.0 };
                        if (cfg->gradientflag == 2) {
                            InterpGradient(cfg, g, v, &pos, &ga, v[i]);
                            fstatus = RFresnel(&u, &ga, &gb, cfg->nv[type], cfg->nv[(int)v[p]], seed, &TIR_events, &status);
                        }
                        else if (cfg->gradientflag == 1) {
                            fstatus = RFresnel(&u, &ga, &gb, cfg->nv[type], cfg->nv[(int)v[p]], seed, &TIR_events, &status);
                        }
                        else if (cfg->gradientflag == 0) {
                            fstatus = RFresnel(&u, &gb, &gb, cfg->nv[type], cfg->nv[(int)v[p]], seed, &TIR_events, &status);
                        }
                        if (fstatus == 0) {
                            pos.x = temp_px + (pos.x - temp_px) * ls * 2;
                            pos.y = temp_py + (pos.y - temp_py) * ls * 2;
                            pos.z = temp_pz + (pos.z - temp_pz) * ls * 2;
                        }
                        else {
                            ix = temp_ix; iy = temp_iy; iz = temp_iz;
                            // update pointer to tissue type
                            type = v[p];
                            pos.x = temp_px;
                            pos.y = temp_py;
                            pos.z = temp_pz;
                            i = p;
                        }
                    } //(sv) /* same voxel */

                } while (u.w > 0.f);

                /**** SPIN
                  Scatter photon into new trajectory defined by theta and psi.
                  Theta is specified by cos(theta), which is determined
                  based on the Henyey-Greenstein scattering function.
                  Convert theta and psi into cosines ux, uy, uz.
                *****/
                /* Sample for costheta */
                rnd = RandomGen(seed);
                float ctheta, stheta, psi, cpsi, spsi;
                if (cfg->gv[type] == 0.0)
                    ctheta = 2.0 * rnd - 1.0;
                else {
                    temp = (1.0 - cfg->gv[type] * cfg->gv[type]) / (1.0 - cfg->gv[type] + 2 * cfg->gv[type] * rnd);
                    ctheta = (1.0 + cfg->gv[type] * cfg->gv[type] - temp * temp) / (2.0 * cfg->gv[type]);
                }
                stheta = sqrt(1.0 - ctheta * ctheta); /* sqrtf() is faster than sin(). */

                /* Sample psi. */
                psi = 2.0 * PI * RandomGen(seed);
                cpsi = cos(psi);
                if (psi < PI)
                    spsi = sqrt(1.0 - cpsi * cpsi);     /* sqrtf() is faster than sin(). */
                else
                    spsi = -sqrt(1.0 - cpsi * cpsi);

                /* New trajectory. */
                if (1 - fabs(u.z) <= ls) {      /* close to perpendicular. */
                    u.x = stheta * cpsi;
                    u.y = stheta * spsi;
                    u.z = ctheta * SIGN(u.z);   /* SIGN() is faster than division. */
                }
                else {					/* usually use this option */
                    temp = sqrt(1.0 - u.z * u.z);
                    float ux, uy, uz;
                    ux = stheta * (u.x * u.z * cpsi - u.y * spsi) / temp + u.x * ctheta;
                    uy = stheta * (u.y * u.z * cpsi + u.x * spsi) / temp + u.y * ctheta;
                    uz = -stheta * cpsi * temp + u.z * ctheta;
                    u.x = ux;
                    u.y = uy;
                    u.z = uz;
                }

                /**** CHECK ROULETTE
                  If photon weight below THRESHOLD, then terminate photon using Roulette technique.
                  Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
                  and 1-CHANCE probability of terminating.
                *****/
                if (pos.w < THRESHOLD) {
                    if (RandomGen(seed) * CHANCE <= 1.f)
                        pos.w *= CHANCE;
                    else
                        status = DEAD;
                }

            } while (status == ALIVE);  /* end STEP_CHECK_HOP_SPIN */
            /* if ALIVE, continue propagating */
            /* If photon DEAD, then launch new photon. */

        } //#pragma omp for
    } //#pragma omp parallel
    return;
}

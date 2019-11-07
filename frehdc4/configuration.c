#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

#include "configuration.h"
#include "utilities.h"

// -----------------------------------------------------------------------------
// This file contains all user-defined parameters for the FREHD model. This
// should be the only file user needs to modify. A complete description of
// all the fields in setting could be found in configuration.h
// - Zhi Li 2017-05-05 -
// -----------------------------------------------------------------------------

// ========= Assign setting values to its fields =========
void ReadUserSettings(Config **setting)
{
  int ii;
  *setting = malloc(sizeof(Config));
  // -------------------- GRID SETTINGS ----------------------------------------



  // UCBV8V4 1x1
  //(*setting)->dx = 1.0;    (*setting)->dy = 1.0;
  //(*setting)->NX = 480;   (*setting)->NY = 1984;
  // UCBV8V4 16x16
  (*setting)->dx = 16.0;    (*setting)->dy = 16.0;
  (*setting)->NX = 30;    (*setting)->NY = 124;
  // DC 1x1
  //(*setting)->dx = 1.0;    (*setting)->dy = 1.0;
  //(*setting)->NX = 200;  (*setting)->NY = 1200;
  //(*setting)->NX = 732;  (*setting)->NY = 720;
  // DC 12x12
  //(*setting)->dx = 12.0;    (*setting)->dy = 12.0;
  //(*setting)->NX = 41;     (*setting)->NY = 60;
  // -------------------- OPERATION SETTINGS -----------------------------------
  // nuecesUpChannelBay
  (*setting)->dt = 2.0;    (*setting)->Nt = 1.0*86400/(*setting)->dt;
  (*setting)->OutItvl = 0.5*3600/(*setting)->dt;
  
  strcpy((*setting)->tStart, "2013-06-10");
  strcpy((*setting)->saveFolder, "outputV8V5MM2CD/output_16x16subAmVmRecOnly3/");
  strcpy((*setting)->inputFolder, "input/");
  (*setting)->useCellEdge = 0;
  (*setting)->savesurface = 1;
  (*setting)->saveuu = 1;
  (*setting)->savevv = 1;
  (*setting)->savedepth = 1;
  (*setting)->savescalar = 1;
  (*setting)->saveCD = 0;
  (*setting)->savesub = 1;
  (*setting)->checkConservation = 0;
  // -------------------- PHYSICAL PROPERTIES ----------------------------------
  (*setting)->g = 9.8066;
  (*setting)->NUx = 0.0001;
  (*setting)->NUy = 0.0001;
  (*setting)->CDnotN = 0;
  (*setting)->CDx = 0.01;
  (*setting)->CDy = 0.01;
  (*setting)->manningN = 0.025;
  (*setting)->z0 = 0.1;
  (*setting)->rhoa = 1.25;
  (*setting)->Cw = 0.0013;
  // -------------------- BOUNDARY CONDITIONS ----------------------------------
  (*setting)->bcType = 1;
  
  // --------- UCBV8V4 1x1 ----------
  /*(*setting)->tideLocLengthM = 10;
  int tideLocM[10];
  for (ii = 0; ii < 10; ii++) {tideLocM[ii] = 160 + ii;}
  (*setting)->tideLocLengthP = 200;
  int tideLocP[200];
  for (ii = 0; ii < 200; ii++) {tideLocP[ii] = 60 + ii;}
  (*setting)->inflowLocLength = 1;
  int jj, kk;
  int inflowLoc[1] = {0};
  jj = 1;*/
  // ---------- UCBV8V4 16x16 ----------
  (*setting)->tideLocLengthM = 1;
  int tideLocM[1];
  for (ii = 0; ii < 1; ii++) {tideLocM[ii] = 10 + ii;}
  (*setting)->tideLocLengthP = 10;
  int tideLocP[10];
  for (ii = 0; ii < 10; ii++) {tideLocP[ii] = ii + 5;}
  (*setting)->inflowLocLength = 1;
  int inflowLoc[1] = {0};
  int jj;
  jj = 1;
    // ------------------------------------
  (*setting)->useWind = 0;
  (*setting)->northAngle = 0;
  (*setting)->windspdN = 119279;
  (*setting)->winddirN = 119279;
  (*setting)->tideNP = 721;//119280;//119280; //721;
  (*setting)->tideNM = 721;
  (*setting)->inflowN = 487;//730;//730; //487;
  (*setting)->useEvap = 0;
  (*setting)->evapN = 100;
  (*setting)->useRain = 0;
  (*setting)->rainN = 100;
  // -------------------- INITIAL CONDITIONS -----------------------------------
  (*setting)->initU = 0;
  (*setting)->initV = 0;
  (*setting)->initSurf = 0.4;//0.445;//0.654;//0.485; //0.557;//0.234;//0.570;//1.173;//1.000;//0.564;//0.600; //1.173;
  (*setting)->useConstSurf0 = 1;
  (*setting)->useConstU0 = 1;
  (*setting)->useConstV0 = 1;
  // -------------------- SCALAR SETTINGS --------------------------------------
  (*setting)->useScalar = 1;
  (*setting)->useConstInitS = 1;
  (*setting)->useConstInflowS = 1;
  (*setting)->useConstTidePS = 1;
  (*setting)->useConstTideMS = 1;
  (*setting)->tidalPSN = 23856;
  (*setting)->tidalMSN = 1;
  (*setting)->initS = 0;
  (*setting)->tidePS = 20;
  (*setting)->tideMS = 0;
  (*setting)->inflowS = 0;
  (*setting)->scalarAdv = 1;
  (*setting)->Kx = 0.0001;    (*setting)->Ky = 0.0001;
  // -------------------- SUBGRID SETTINGS -------------------------------------
  (*setting)->useSubgrid = 1;
  (*setting)->useSubDrag = 0;
  (*setting)->lambda1 = 0.0;
  (*setting)->lambda2 = 0.0;
  (*setting)->useCorrector = 0;
  strcpy((*setting)->subgridFolder, "subdata_UCBV8V5MM216x16/");
  (*setting)->dxf = 1.0;   (*setting)->dyf = 1.0;
  (*setting)->dA = (*setting)->dxf * (*setting)->dyf;
  (*setting)->subR = (*setting)->dx / (*setting)->dxf;
  (*setting)->surfmax = 0.75;//0.4;//0.8;//1.0; //0.8;
  (*setting)->surfmin = 0.25;//0.0;//-0.2;//-0.2; //0.2;
  (*setting)->dsurf = 0.01;//0.005;//0.002;
  (*setting)->staggeredV = 1;
  (*setting)->phiSurface = 0;
  (*setting)->advChannel = 0;
  (*setting)->useminA = 1;
  (*setting)->beta = 1.0;
  (*setting)->useminV = 3;
  (*setting)->Cv = 0.0;
  (*setting)->b = 0.0;
  (*setting)->recircZone = 0;
  (*setting)->NindN = 4;
  (*setting)->OindN = 3;
  int minNind[4] = {3273,3338,3403,3470};
  int minOind[3] = {3404,3405,3471};
  (*setting)->divider = 1.0;
  // -------------------- MPI SETTINGS -----------------------------------------
  (*setting)->useMPI = 1;
  (*setting)->np = 2;
  // -------------------- PERFORMANCE SETTINGS ---------------------------------
  (*setting)->minDepth = 0.00001;
  (*setting)->eps = 0.00000001;
  (*setting)->maxIter = 10000000;
  (*setting)->useThinLayer = 0;
  (*setting)->CDmax = 1.0;
  (*setting)->CwT = 10.0;
  (*setting)->hD = (*setting)->z0;
  (*setting)->wtfh = 0.00001;
  (*setting)->CFLl = 0.5;
  (*setting)->CFLh = 0.7;
  // ---------------------------------------------------------------------------
  // -------------------- NO USER CHANGE BELOW ---------------------------------
  // ---------------------------------------------------------------------------
  // ---------- Derived Values ----------
  // sizes of data arrays
  (*setting)->N2CI = (*setting)->NX*(*setting)->NY;
  (*setting)->N2CT = ((*setting)->NX+2)*((*setting)->NY+2);
  if ((*setting)->NY % (*setting)->np != 0)
  {printf("WARNING: The computation domain cannot be evenly divided by nranks!\n");}
  else
  {
    (*setting)->nx = (*setting)->NX;    (*setting)->ny = round((*setting)->NY/(*setting)->np);
    (*setting)->N2ci = (*setting)->nx*(*setting)->ny;
    (*setting)->N2ct = ((*setting)->nx+2)*((*setting)->ny+2);
  }
  // start and end time represented by date number
  (*setting)->tNStart = dateNum((*setting)->tStart);
  (*setting)->tNEnd = (*setting)->tNStart + (*setting)->dt*(*setting)->Nt;
  // locations of boundary conditions
  for (ii = 0; ii < (*setting)->tideLocLengthP; ii++)
  {(*setting)->tideLocP[ii] = tideLocP[ii] + (*setting)->nx*((*setting)->ny-1);}
  for (ii = 0; ii < (*setting)->tideLocLengthM; ii++)
  {(*setting)->tideLocM[ii] = tideLocM[ii];}
  (*setting)->inflow_rank = floor(jj / ((*setting)->NY / (*setting)->np));
  //printf("inflow rank = %d\n",(*setting)->inflow_rank);
  for (ii = 0; ii < (*setting)->inflowLocLength; ii++)
    {(*setting)->inflowLoc[ii] = inflowLoc[ii] - (*setting)->inflow_rank * (*setting)->nx * (*setting)->ny;}// + (*setting)->nx;}
  for (ii = 0; ii < (*setting)->NindN; ii++)
    {(*setting)->minNind[ii] = minNind[ii];}
  for (ii = 0; ii < (*setting)->OindN; ii++)
    {(*setting)->minOind[ii] = minOind[ii];}
  // check if subgrid bathymetry can be created
  if ((*setting)->useSubgrid == 1)
  {
    if (fmod((*setting)->dx,(*setting)->dxf) != 0 | fmod((*setting)->dy,(*setting)->dyf) != 0)
    {printf("WARNING: The coarse grid cannot be divided by the fine grid!\n");}
  }
}

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
  (*setting)->dx = 10.0;    (*setting)->dy = 10.0;
  (*setting)->NX = 60;  (*setting)->NY = 80;
  // -------------------- OPERATION SETTINGS -----------------------------------
  (*setting)->dt = 1.0;    (*setting)->Nt = (1.0/3.0)*86400/(*setting)->dt;
  (*setting)->OutItvl = 0.1*3600/(*setting)->dt;
  strcpy((*setting)->tStart, "2013-06-10");
  strcpy((*setting)->saveFolder, "WeirV5noinflow/output_10x10sub/");
  strcpy((*setting)->inputFolder, "WeirV5noinflow/input_10x10/");
  (*setting)->isRestart = 0;
  (*setting)->ttRestart = 43200;
  strcpy((*setting)->restartFile, "StairV4/output_spinup/restart_43200.dat");
  (*setting)->useCellEdge = 0;
  (*setting)->savesurface = 1;
  (*setting)->saveuu = 1;
  (*setting)->savevv = 1;
  (*setting)->savedepth = 1;
  (*setting)->savescalar = 1;
  (*setting)->saveCD = 1;
  (*setting)->savesub = 1;
  // -------------------- PHYSICAL PROPERTIES ----------------------------------
  (*setting)->g = 9.8066;
  (*setting)->NUx = 0.0001;
  (*setting)->NUy = 0.0001;
  (*setting)->CDnotN = 0;
  (*setting)->CDx = 0.01;
  (*setting)->CDy = 0.01;
  (*setting)->manningN = 0.03;
  (*setting)->z0 = 0.1;
  (*setting)->rhoa = 1.225;
  (*setting)->Cw = 0.0013;
  // -------------------- BOUNDARY CONDITIONS ----------------------------------
  (*setting)->bcType = 1;
  (*setting)->tideLocLengthP = 40;
  int tideLocP[40];
  for (ii = 0; ii < 40; ii++) {tideLocP[ii] = 10 + ii;}
  (*setting)->tideLocLengthM = 7;
  int tideLocM[7];
  for (ii = 0; ii < 7; ii++) {tideLocM[ii] = 18 + ii;}
  (*setting)->inflowLocLength = 7;
  int inflowLoc[7] = {18,19,20,21,22,23,24};
  
  (*setting)->useWind = 0;
  (*setting)->northAngle = 0;
  (*setting)->windspdN = 119279;
  (*setting)->winddirN = 119279;
  (*setting)->tideNP = 721;//119280; //721;
  (*setting)->tideNM = 1;
  (*setting)->inflowN = 487;//730; //487;
  (*setting)->useEvap = 0;
  (*setting)->evapN = 100;
  (*setting)->useRain = 0;
  (*setting)->rainN = 100;
  // -------------------- INITIAL CONDITIONS -----------------------------------
  (*setting)->initU = 0;
  (*setting)->initV = 0;
  (*setting)->initSurf = 1.087;//0.500;//1.000;//0.564;//0.600; //1.173;
  // -------------------- SCALAR SETTINGS --------------------------------------
  (*setting)->useScalar = 1;
  (*setting)->useConstInitS = 0;
  (*setting)->useConstInflowS = 1;
  (*setting)->useConstTidePS = 1;
  (*setting)->useConstTideMS = 1;
  (*setting)->tidalPSN = 23856;
  (*setting)->tidalMSN = 1;
  (*setting)->initS = 10;
  (*setting)->tidePS = 25;
  (*setting)->tideMS = 0;
  (*setting)->inflowS = 0;
  (*setting)->Kx = 0.0001;    (*setting)->Ky = 0.0001;
  // -------------------- SUBGRID SETTINGS -------------------------------------
  (*setting)->useSubgrid = 1;
  (*setting)->useSubDrag = 0;
  (*setting)->useCorrector = 0;
  strcpy((*setting)->subgridFolder, "subdata/");
  (*setting)->dxf = 1.0;   (*setting)->dyf = 1.0;
  (*setting)->dA = (*setting)->dxf * (*setting)->dyf;
  (*setting)->subR = (*setting)->dx / (*setting)->dxf;
  (*setting)->surfmax = 1.3;//0.9;//1.3;//1.0; //0.8;
  (*setting)->surfmin = 0.7;//0.1;//0.7;//-0.2; //0.2;
  (*setting)->dsurf = 0.005;//0.005;//0.002;
  (*setting)->subAlpha = 0.0;
  // -------------------- MPI SETTINGS -----------------------------------------
  (*setting)->useMPI = 0;
  (*setting)->np = 1;
  // -------------------- PERFORMANCE SETTINGS ---------------------------------
  (*setting)->minDepth = 0.002;
  (*setting)->eps = 0.000001;
  (*setting)->maxIter = 100000;
  (*setting)->useThinLayer = 1;
  (*setting)->CDmax = 1;
  (*setting)->CwT = 5;
  (*setting)->hD = (*setting)->z0;
  (*setting)->wtfh = 0.002;
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
  for (ii = 0; ii < (*setting)->inflowLocLength; ii++)
  {(*setting)->inflowLoc[ii] = inflowLoc[ii];}// + (*setting)->nx;}
  // check if subgrid bathymetry can be created
  if ((*setting)->useSubgrid == 1)
  {
    if (fmod((*setting)->dx,(*setting)->dxf) != 0 | fmod((*setting)->dy,(*setting)->dyf) != 0)
    {printf("WARNING: The coarse grid cannot be divided by the fine grid!\n");}
  }
}

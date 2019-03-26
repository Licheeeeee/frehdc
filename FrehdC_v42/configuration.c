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
void ReadUserSettings(Config **setting, double *configArr, char *date, char *inputFolder, char *outputFolder, char *subgridFolder)
{
    int ii, jj, kk = 0;
    *setting = malloc(sizeof(Config));
    // -------------------- GRID SETTINGS ----------------------------------------
    (*setting)->dx = configArr[4];    (*setting)->dy = configArr[5];
    (*setting)->NX = configArr[6];  (*setting)->NY = configArr[7];
    // -------------------- OPERATION SETTINGS -----------------------------------
    (*setting)->dt = configArr[8];    (*setting)->Nt = configArr[9]*86400/(*setting)->dt;
    (*setting)->OutItvl = configArr[11]*3600/(*setting)->dt;
    strcpy((*setting)->tStart, date);
    strcpy((*setting)->saveFolder, outputFolder);
    strcpy((*setting)->inputFolder, inputFolder);
    (*setting)->useCellEdge = configArr[10];
    (*setting)->savesurface = configArr[12];
    (*setting)->saveuu = configArr[13];
    (*setting)->savevv = configArr[14];
    (*setting)->savedepth = configArr[15];
    (*setting)->savescalar = configArr[16];
    (*setting)->saveCD = configArr[17];
    (*setting)->savesub = configArr[18];
    // -------------------- BOUNDARY CONDITIONS ----------------------------------
    (*setting)->bcType = configArr[19];
    // --- Tide P ---
    (*setting)->tideLocLengthP = configArr[20];;
    int tideLocP[(*setting)->tideLocLengthP];
    for (ii = 0; ii < (*setting)->tideLocLengthP; ii++)
    {tideLocP[ii] = configArr[21] + ii;}
    (*setting)->tideNP = configArr[22];
    // --- Tide M ---
    (*setting)->tideLocLengthM = configArr[23];;
    int tideLocM[(*setting)->tideLocLengthM];
    for (ii = 0; ii < (*setting)->tideLocLengthM; ii++)
    {tideLocM[ii] = configArr[24] + ii;}
    (*setting)->tideNM = configArr[25];
    // --- Inflow ---
    (*setting)->inflowLocLength = configArr[26];
    int inflowLoc[(*setting)->inflowLocLength];
    for (ii = configArr[27]; ii < configArr[28]; ii++)
    {
        for (jj = configArr[29]; jj < configArr[30]; jj++)
        {inflowLoc[kk] = jj * (*setting)->NX + ii - 1;      kk += 1;}
    }
    (*setting)->inflowN = configArr[31];
    // --- Wind ---
    (*setting)->useWind = configArr[32];
    (*setting)->northAngle = configArr[33];
    (*setting)->windspdN = configArr[34];
    (*setting)->winddirN = configArr[35];
    // --- Scalar ---
    (*setting)->useScalar = configArr[36];
    (*setting)->scalarAdv = configArr[37];
    (*setting)->useConstInflowS = configArr[38];
    (*setting)->inflowS = configArr[39];
    (*setting)->useConstTidePS = configArr[40];
    (*setting)->tidePS = configArr[41];
    (*setting)->tidalPSN = configArr[42];
    (*setting)->useConstTideMS = configArr[43];
    (*setting)->tideMS = configArr[44];
    (*setting)->tidalMSN = configArr[45];
    // -------------------- INITIAL CONDITIONS -----------------------------------
    (*setting)->initU = configArr[46];
    (*setting)->initV = configArr[47];
    (*setting)->initSurf = configArr[48];
    (*setting)->initS = configArr[49];
    (*setting)->useConstSurf0 = configArr[50];
    (*setting)->useConstU0 = configArr[51];
    (*setting)->useConstV0 = configArr[52];
    (*setting)->useConstInitS = configArr[53];
    // -------------------- SUBGRID SETTINGS -------------------------------------
    (*setting)->useSubgrid = configArr[54];
    (*setting)->useSubDrag = configArr[55];
    strcpy((*setting)->subgridFolder, subgridFolder);
    (*setting)->dxf = configArr[56];
    (*setting)->dyf = configArr[57];
    (*setting)->dA = (*setting)->dxf * (*setting)->dyf;
    (*setting)->subR = (*setting)->dx / (*setting)->dxf;
    (*setting)->surfmax = configArr[58];
    (*setting)->surfmin = configArr[59];
    (*setting)->dsurf = configArr[60];
    (*setting)->staggeredV = configArr[61];
    (*setting)->useminA = configArr[62];
    (*setting)->beta = configArr[63];
    (*setting)->useminV = configArr[64];
    // -------------------- MPI SETTINGS -----------------------------------------
    (*setting)->useMPI = configArr[65];
    (*setting)->np = configArr[66];
    // -------------------- PHYSICAL PROPERTIES ----------------------------------
    (*setting)->g = configArr[67];
    (*setting)->NUx = configArr[68];
    (*setting)->NUy = configArr[69];
    (*setting)->CDnotN = configArr[70];
    (*setting)->CDx = configArr[71];
    (*setting)->CDy = configArr[72];
    (*setting)->manningN = configArr[73];
    (*setting)->useThinLayer = configArr[74];
    (*setting)->z0 = configArr[75];
    (*setting)->hD = (*setting)->z0;
    (*setting)->CDmax = configArr[76];
    (*setting)->CwT = configArr[77];
    (*setting)->rhoa = configArr[78];
    (*setting)->Cw = configArr[79];
    (*setting)->Kx = configArr[80];
    (*setting)->Ky = configArr[81];
    (*setting)->minDepth = configArr[82];
    (*setting)->wtfh = configArr[83];
    (*setting)->CFLl = configArr[84];
    (*setting)->CFLh = configArr[85];
    (*setting)->eps = configArr[86];
    (*setting)->maxIter = configArr[87];
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

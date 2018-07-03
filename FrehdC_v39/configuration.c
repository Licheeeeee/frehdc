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
  // CCWD 1x1
  //(*setting)->dx = 1.0;    (*setting)->dy = 1.0;
  //(*setting)->NX = 600;  (*setting)->NY = 800;
  // CCWD 4x4
  //(*setting)->dx = 4.0;    (*setting)->dy = 4.0;
  //(*setting)->NX = 150;  (*setting)->NY = 200;
  // CCWD 10x10
  //(*setting)->dx = 10.0;    (*setting)->dy = 10.0;
  //(*setting)->NX = 60;  (*setting)->NY = 80;
  // --------------- UCBV3 ---------------
  // nuecesUCBV3 8x8
  // (*setting)->dx = 8.0;    (*setting)->dy = 8.0;
  // (*setting)->NX = 76;  (*setting)->NY = 348;
  // nuecesUCBV3 16x16
   (*setting)->dx = 16.0;    (*setting)->dy = 16.0;
   (*setting)->NX = 38;  (*setting)->NY = 174;
  // --------------- UCBV4 ---------------
  // nuecesUCBV4 8x8
//   (*setting)->dx = 8.0;    (*setting)->dy = 8.0;
//   (*setting)->NX = 104;  (*setting)->NY = 162;
  // nuecesUCBV4 16x16
//  (*setting)->dx = 16.0;    (*setting)->dy = 16.0;
//  (*setting)->NX = 52;  (*setting)->NY = 81;
  // -------------------- OPERATION SETTINGS -----------------------------------
  // Stair and Weir
  //(*setting)->dt = 1.0;    (*setting)->Nt = (1.0/3.0)*86400/(*setting)->dt;
  //(*setting)->OutItvl = 0.1*3600/(*setting)->dt;
  // nuecesUpChannelBay
  (*setting)->dt = 0.5;    (*setting)->Nt = 1.5*86400/(*setting)->dt;
  (*setting)->OutItvl = 0.5*3600/(*setting)->dt;
  // ChannelBlockage
  //(*setting)->dt = 1.0;    (*setting)->Nt = 0.5*86400/(*setting)->dt;
  //(*setting)->OutItvl = 0.1*3600/(*setting)->dt;
  // NDHM
  //(*setting)->dt = 10.0;    (*setting)->Nt = 61*86400/(*setting)->dt;
  //(*setting)->OutItvl = 2*3600/(*setting)->dt;
  strcpy((*setting)->tStart, "2013-06-10");
  strcpy((*setting)->saveFolder, "nuecesUCBV3/output_16x16mpidebug/");
  strcpy((*setting)->inputFolder, "nuecesUCBV3/input_16x16/");
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
  (*setting)->checkConservation = 1;
  // -------------------- PHYSICAL PROPERTIES ----------------------------------
  (*setting)->g = 9.8066;
  (*setting)->NUx = 0.0001;
  (*setting)->NUy = 0.0001;
  (*setting)->CDnotN = 0;
  (*setting)->CDx = 0.01;
  (*setting)->CDy = 0.01;
  (*setting)->manningN = 0.05;
  (*setting)->z0 = 0.1;
  (*setting)->rhoa = 1.225;
  (*setting)->Cw = 0.0013;
  // -------------------- BOUNDARY CONDITIONS ----------------------------------
  (*setting)->bcType = 1;
  // ---------- for CCWD 1x1 ----------
  /*(*setting)->tideLocLengthP = 400;
  int tideLocP[400], inflowLoc[70];
  for (ii = 0; ii < 400; ii++) {tideLocP[ii] = 100 + ii;}
  (*setting)->tideLocLengthM = 7;
  int tideLocM[7];
  for (ii = 0; ii < 7; ii++) {tideLocM[ii] = 18 + ii;}
  (*setting)->inflowLocLength = 70;
  for (ii = 0; ii < 70; ii++) {inflowLoc[ii] = 181 + ii;}*/
  // ---------- for CCWD 4x4 ----------
  /*(*setting)->tideLocLengthP = 100;
  int tideLocP[100];
  for (ii = 0; ii < 100; ii++) {tideLocP[ii] = 25 + ii;}
  (*setting)->tideLocLengthM = 7;
  int tideLocM[7];
  for (ii = 0; ii < 7; ii++) {tideLocM[ii] = 18 + ii;}
  (*setting)->inflowLocLength = 17;
  int inflowLoc[17] = {45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61};*/
  // ---------- for CCWD 10x10 ----------
  /*(*setting)->tideLocLengthP = 40;
  int tideLocP[40];
  for (ii = 0; ii < 40; ii++) {tideLocP[ii] = 10 + ii;}
  (*setting)->tideLocLengthM = 7;
  int tideLocM[7];
  for (ii = 0; ii < 7; ii++) {tideLocM[ii] = 18 + ii;}
  (*setting)->inflowLocLength = 7;
  int inflowLoc[7] = {18,19,20,21,22,23,24};*/
  //int inflowLoc[7] = {78,79,80,81,82,83,84};*/
  // ---------- nuecesUCBV3 8x8 ----------
  // (*setting)->tideLocLengthP = 50;
  // int tideLocP[50];
  // for (ii = 0; ii < 50; ii++) {tideLocP[ii] = 10 + ii;}
  // ---------- nuecesUCBV3 16x16 ----------
   (*setting)->tideLocLengthP = 25;
   int tideLocP[25];
   for (ii = 0; ii < 25; ii++) {tideLocP[ii] = 5 + ii;}
  // ------------------------------------
  (*setting)->tideLocLengthM = 1;
  int tideLocM[1];
  for (ii = 0; ii < 1; ii++) {tideLocM[ii] = 1 + ii;}
   (*setting)->inflowLocLength = 1;
   int inflowLoc[1] = {1};
  // ---------- nuecesUCBV4 8x8 ----------
//   (*setting)->tideLocLengthP = 5;
//   int tideLocP[5];
//   for (ii = 0; ii < 5; ii++) {tideLocP[ii] = 29 + ii;}
//   (*setting)->inflowLocLength = 8;
//   int inflowLoc[8] = {1464,1465,1568,1569,1550,1551,1654,1655};
  // ---------- nuecesUCBV4 16x16 ----------
//  (*setting)->tideLocLengthP = 2;
//  int tideLocP[2];
//  for (ii = 0; ii < 2; ii++) {tideLocP[ii] = 15 + ii;}
//  (*setting)->inflowLocLength = 2;
//  int inflowLoc[2] = {368,411}; // [5,8], [48,8]
  // ------------------------------------
  (*setting)->useWind = 0;
  (*setting)->northAngle = 0;
  (*setting)->windspdN = 721;//119279
  (*setting)->winddirN = 721;
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
  (*setting)->initSurf = 0.500;//1.087;//0.500;//1.000;//0.600; //1.173;
  // -------------------- SCALAR SETTINGS --------------------------------------
  (*setting)->useScalar = 1;
  (*setting)->useConstInitS = 0;
  (*setting)->useConstInflowS = 1;
  (*setting)->useConstTidePS = 1;
  (*setting)->useConstTideMS = 1;
  (*setting)->tidalPSN = 23856;
  (*setting)->tidalMSN = 1;
  (*setting)->initS = 45;
  (*setting)->tidePS = 45; //45 for nuecesUCBV2, 25 for Weir
  (*setting)->tideMS = 0;
  (*setting)->inflowS = 0;
  (*setting)->Kx = 0.0001;    (*setting)->Ky = 0.0001;
  (*setting)->scalarAdv = 1;
  // -------------------- SUBGRID SETTINGS -------------------------------------
  (*setting)->useSubgrid = 1; // 1 for subMass+Mome, 2 for subMass only
  (*setting)->useSubDrag = 0; // 1 for curvature-based Cd, 2 for Volp2013, 3 for ZhiLi20180329, 4 for Casas2010
  (*setting)->lambda1 = 0.0;
  (*setting)->lambda2 = 0.0;
  (*setting)->useCorrector = 0;
  strcpy((*setting)->subgridFolder, "subdata/");
  (*setting)->dxf = 1.0;   (*setting)->dyf = 1.0;
  (*setting)->dA = (*setting)->dxf * (*setting)->dyf;
  (*setting)->subR = (*setting)->dx / (*setting)->dxf;
    (*setting)->surfmax = 0.9;//0.6;//0.9;//1.3;//1.0; //0.8;
    (*setting)->surfmin = 0.1;//0.4;//0.1;//0.7;//-0.2; //0.2;
  (*setting)->dsurf = 0.005;//0.005;//0.002;
  (*setting)->staggeredV = 1;
  (*setting)->useAVcutoff = 0;
  (*setting)->effHmin = 0.99*(1.0/(*setting)->dx);
  // -------------------- MPI SETTINGS -----------------------------------------
  (*setting)->useMPI = 1;
  (*setting)->np = 2;
  // -------------------- PERFORMANCE SETTINGS ---------------------------------
  (*setting)->minDepth = 0.00001;
  (*setting)->eps = 0.00000001;
  (*setting)->maxIter = 10000000;
  (*setting)->useThinLayer = 1;
  (*setting)->CDmax = 1.0;
  (*setting)->CwT = 2.0;
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
  for (ii = 0; ii < (*setting)->inflowLocLength; ii++)
  {(*setting)->inflowLoc[ii] = inflowLoc[ii];}// + (*setting)->nx;}
  // check if subgrid bathymetry can be created
  if ((*setting)->useSubgrid == 1)
  {
    if (fmod((*setting)->dx,(*setting)->dxf) != 0 | fmod((*setting)->dy,(*setting)->dyf) != 0)
    {printf("WARNING: The coarse grid cannot be divided by the fine grid!\n");}
  }
}

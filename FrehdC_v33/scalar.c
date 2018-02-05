#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

// -----------------------------------------------------------------------------
// This file contains functions used for scalar computations
// - ZhiLi 2017-05-13 -
// -----------------------------------------------------------------------------

#include "configuration.h"
#include "bathymetry.h"
#include "map.h"
#include "fileio.h"
#include "mpifunctions.h"
#include "scalar.h"
#include "subgrid.h"
#include "utilities.h"

void scalarTransport(Data **data, Maps *map, Bath *bath, BC *bc, Sub *sub, Config *setting, int tt, int irank, int nrank);
void scalarMass(Data **data, Config *setting);
void scalarAdvection(Data **data, Maps *map, Config *setting);
void scalarDiffusion(Data **data, Maps *map, Config *setting);
void updateScalar(Data **data, BC *bc, Bath *bath, Maps *map, Config *setting, int tt, int irank, int nrank);
void enforceScalarBC(Data **data, BC *bc, Maps *map, Config *setting, int tt, int irank, int nrank);

// ==================== Scalar top level script ====================
void scalarTransport(Data **data, Maps *map, Bath *bath, BC *bc, Sub *sub, Config *setting, int tt, int irank, int nrank)
{
  /*int ii = 21, jj = ii + 1, kk = map->icjM[ii];
  if (irank == 1)
  {
    printf("scalar = %lf,%lf,%lf\n",(*data)->S[ii],(*data)->S[jj],(*data)->S[kk]);
    printf("scalarmass = %lf,%lf,%lf\n",(*data)->Sm[ii],(*data)->Sm[jj],(*data)->Sm[kk]);
  }*/
  scalarAdvection(data, map, setting);
  if (setting->useSubgrid != 0)
  {
    scalarDiffusionSub(data, map, sub, setting);
    updateScalarSub(data, bc, bath, map, sub, setting, tt, irank, nrank);
  }
  else
  {
    scalarDiffusion(data, map, setting);
    updateScalar(data, bc, bath, map, setting, tt, irank, nrank);
  }
  enforceScalarBC(data, bc, map, setting, tt, irank, nrank);
  //scalarMass(data, setting);
  if (setting->useMPI == 1)
  {
    mpiexchange((*data)->S, map, setting, irank, nrank);
    mpiexchange((*data)->Sm, map, setting, irank, nrank);
  }
  int ii = 624, jj = ii+60, kk = jj+60, pp = ii+1, qq = jj+1;
  //printf("scalar = %lf,%lf,%lf\n",(*data)->S[ii],(*data)->S[jj],(*data)->S[pp]);
}

// =============== compute mass of scalar in each cell ===============
void scalarMass(Data **data, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {(*data)->Sm[ii] = (*data)->S[ii] * (*data)->cellV[ii];}
}

// ==================== advective scalar flux ====================
void scalarAdvection(Data **data, Maps *map, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    // advective flux in x direction
    (*data)->Sm[ii] = (*data)->Sm[ii] + setting->dt * \
      (-0.5 * ((*data)->Fuu[ii] + fabs((*data)->Fuu[ii])) * (*data)->S[ii] - \
      0.5 * ((*data)->Fuu[ii] - fabs((*data)->Fuu[ii])) * (*data)->S[map->iPjc[ii]] +\
      0.5 * ((*data)->Fuu[map->iMjc[ii]] + fabs((*data)->Fuu[map->iMjc[ii]])) * (*data)->S[map->iMjc[ii]] +\
      0.5 * ((*data)->Fuu[map->iMjc[ii]] - fabs((*data)->Fuu[map->iMjc[ii]])) * (*data)->S[ii]);
    // advective flux in y direction
    (*data)->Sm[ii] = (*data)->Sm[ii] + setting->dt * \
      (-0.5 * ((*data)->Fvv[ii] + fabs((*data)->Fvv[ii])) * (*data)->S[ii] - \
      0.5 * ((*data)->Fvv[ii] - fabs((*data)->Fvv[ii])) * (*data)->S[map->icjP[ii]] +\
      0.5 * ((*data)->Fvv[map->icjM[ii]] + fabs((*data)->Fvv[map->icjM[ii]])) * (*data)->S[map->icjM[ii]] +\
      0.5 * ((*data)->Fvv[map->icjM[ii]] - fabs((*data)->Fvv[map->icjM[ii]])) * (*data)->S[ii]);
  }
}

// ==================== diffusive scalar flux ====================
void scalarDiffusion(Data **data, Maps *map, Config *setting)
{
  int ii;
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    (*data)->Sm[ii] = (*data)->Sm[ii] + setting->dt * \
      ((setting->Kx * setting->dy * (*data)->depthXP[ii] / setting->dx) * \
      ((*data)->S[map->iPjc[ii]] - (*data)->S[ii]) - \
      (setting->Kx * setting->dy * (*data)->depthXP[map->iMjc[ii]] / setting->dx) * \
      ((*data)->S[ii] - (*data)->S[map->iMjc[ii]]) + \
      (setting->Ky * setting->dx * (*data)->depthYP[ii] / setting->dy) * \
      ((*data)->S[map->icjP[ii]] - (*data)->S[ii]) - \
      (setting->Ky * setting->dx * (*data)->depthYP[map->icjM[ii]] / setting->dy) * \
      ((*data)->S[ii] - (*data)->S[map->icjM[ii]]));
  }
}

// ================ update scalar concentration ===============
void updateScalar(Data **data, BC *bc, Bath *bath, Maps *map, Config *setting, int tt, int irank, int nrank)
{
  int ii, jj, *N, flag = 0;
  int iP = 0, iM = 0, jP = 0, jM = 0;
  double *Sarr, *Smax, *Smin;
  double dXPold, dXMold, dYPold, dYMold;
  // local scalar limiter
  N = malloc(setting->N2ci * sizeof(int));
  Smax = malloc(setting->N2ci * sizeof(double));
  Smin = malloc(setting->N2ci * sizeof(double));
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    N[ii] = 0, Smax[ii] = 0, Smin[ii] = 0, jj = 0;
    // check if the target cell is connected to its neighbors
    // comment out this section to disable the limiter, ZhiLi 20170607
    dXPold = (*data)->surfOld[map->iPjc[ii]] - bath->bottomZ[map->iPjc[ii]];
    //if (dXPold > 0) {N[ii] += 1; iP = 1; dXPold = 0;}
    if ((*data)->depthXP[ii] > 0 & dXPold > 0) {N[ii] += 1; iP = 1; dXPold = 0;}

    dXMold = (*data)->surfOld[map->iMjc[ii]] - bath->bottomZ[map->iMjc[ii]];
    //if (dXMold > 0) {N[ii] += 1; iM = 1; dXMold = 0;}
    if ((*data)->depthXP[map->iMjc[ii]] > 0 & dXMold > 0) {N[ii] += 1; iM = 1; dXMold = 0;}

    if (irank != nrank - 1)
    {
      dYPold = (*data)->surfOld[map->icjP[ii]] - bath->bottomZ[map->icjP[ii]];
      //if (dYPold > 0) {N[ii] += 1; jP = 1; dYPold = 0;}
      if ((*data)->depthYP[ii] > 0 & dYPold > 0) {N[ii] += 1; jP = 1; dYPold = 0;}
    }
    else
    {
      if (map->icjP[ii] <= setting->N2ci)
      {
        dYPold = (*data)->surfOld[map->icjP[ii]] - bath->bottomZ[map->icjP[ii]];
        //if (dYPold > 0) {N[ii] += 1; jP = 1; dYPold = 0;}
        if ((*data)->depthYP[ii] > 0 & dYPold > 0) {N[ii] += 1; jP = 1; dYPold = 0;}
      }
    }

    if (irank != 0)
    {
      dYMold = (*data)->surfOld[map->icjM[ii]] - bath->bottomZ[map->icjM[ii]];
      //if (dYMold > 0) {N[ii] += 1; jM = 1; dYMold = 0;}
      if ((*data)->depthYP[map->icjM[ii]] > 0 & dYMold > 0) {N[ii] += 1; jM = 1; dYMold = 0;}
    }
    else
    {
      if (map->icjM[ii] <= setting->N2ci)
      {
        dYMold = (*data)->surfOld[map->icjM[ii]] - bath->bottomZ[map->icjM[ii]];
        //if (dYMold > 0) {N[ii] += 1; jM = 1; dYMold = 0;}
        if ((*data)->depthYP[map->icjM[ii]] > 0 & dYMold > 0) {N[ii] += 1; jM = 1; dYMold = 0;}
      }
    }


    /*if ((*data)->depthXP[ii] > 0) {N[ii] += 1; iP = 1;}
    if ((*data)->depthXP[map->iMjc[ii]] > 0) {N[ii] += 1; iM = 1;}
    if ((*data)->depthYP[ii] > 0) {N[ii] += 1; jP = 1;}
    if ((*data)->depthYP[map->icjM[ii]] > 0) {N[ii] += 1; jM = 1;}*/
    // if yes, get the max and min values of its neighbors

    if (N[ii] != 0)
    {
      Sarr = malloc(N[ii] * sizeof(double));
      if (iP == 1) {Sarr[jj] = (*data)->S[map->iPjc[ii]]; jj++;}
      if (iM == 1) {Sarr[jj] = (*data)->S[map->iMjc[ii]]; jj++;}
      if (jP == 1) {Sarr[jj] = (*data)->S[map->icjP[ii]]; jj++;}
      if (jM == 1) {Sarr[jj] = (*data)->S[map->icjM[ii]]; jj++;}
      Smax[ii] = getMax(Sarr, N[ii]);
      Smin[ii] = getMin(Sarr, N[ii]);
      iP = 0;
      iM = 0;
      jP = 0;
      jM = 0;
      free(Sarr);
    }
  }
  // scalar from inflow
  if (irank == 0 & setting->bcType != 3)
  {
    for (ii = 0; ii < setting->inflowLocLength; ii++)
    {(*data)->Sm[setting->inflowLoc[ii]] = (*data)->Sm[setting->inflowLoc[ii]] + \
      setting->inflowS * bc->inflow[tt] * setting->dt / setting->inflowLocLength;}
  }
  // update scalar concentration
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->cellV[ii] > 0)
    {(*data)->S[ii] = (*data)->Sm[ii] / (*data)->cellV[ii];}
    else
    {(*data)->S[ii] = 0;}
    // apply the limiter
    if (irank == 0 & setting->bcType != 3)
    {
      for (jj = 0; jj < setting->inflowLocLength; jj++)
      {if (ii == setting->inflowLoc[jj]) {flag = 1; break;}}
    }
    // Note: The scalar limiter was not compatible with MPI, so it was disabled
    // by ZhiLi 20170607
    // debugged to allow salinity transferred across jP and jM boundary
    // by ZhiLi 20170815

    if (N[ii] != 0 & (*data)->S[ii] > Smax[ii] & flag == 0)
    {
      /*if (map->icjP[ii] < setting->N2ci & map->icjM[ii] < setting->N2ci)
      {(*data)->S[ii] = Smax[ii];}
      else if (irank == nrank-1 & map->icjP[ii] >= setting->N2ci)
      {(*data)->S[ii] = Smax[ii];}
      else if (irank == 0 & map->icjM[ii] >= setting->N2ci)
      {(*data)->S[ii] = Smax[ii];}*/
      (*data)->S[ii] = Smax[ii];
    }
    else if (N[ii] != 0 & (*data)->S[ii] < Smin[ii] & flag == 0)
    {
      /*if (map->icjP[ii] < setting->N2ci & map->icjM[ii] < setting->N2ci)
      {(*data)->S[ii] = Smin[ii];}
      else if (irank == nrank-1 & map->icjP[ii] >= setting->N2ci)
      {(*data)->S[ii] = Smin[ii];}
      else if (irank == 0 & map->icjM[ii] >= setting->N2ci)
      {(*data)->S[ii] = Smin[ii];}*/
      (*data)->S[ii] = Smin[ii];
    }
    flag = 0;
  }

  free(N);
  free(Smax);
  free(Smin);
  // remove high scalar singularities
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->S[ii] > 100 | (*data)->S[ii] < 0)
    {
      if (!(irank == 0 & map->icjM[ii] > setting->N2ci) & \
       !(irank == nrank - 1 & map->icjP[ii] > setting->N2ci))
      {
        (*data)->S[ii] = 0;
        (*data)->surf[ii] = (*data)->surf[ii] - (*data)->depth[ii];
        (*data)->depth[ii] = 0;
      }
    }
  }
}

// =============== enforce boundary condition ===============
void enforceScalarBC(Data **data, BC *bc, Maps *map, Config *setting, int tt, \
  int irank, int nrank)
{
  int ii;
  // Neumann bc at all boundaries
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*data)->S[map->jPgt[ii]] = (*data)->S[map->jPbd[ii]];
    (*data)->S[map->jMgt[ii]] = (*data)->S[map->jMbd[ii]];
  }
  for (ii = 0; ii < setting->ny; ii++)
  {
    (*data)->S[map->iPgt[ii]] = (*data)->S[map->iPbd[ii]];
    (*data)->S[map->iMgt[ii]] = (*data)->S[map->iMbd[ii]];
  }
  // Direchlet bc at YP tidal boundaries
  if (irank == nrank - 1 & setting->bcType != 2)
  {
    if (setting->useConstTidePS == 1)
    {
      for (ii = 0; ii < setting->tideLocLengthP; ii++)
      {(*data)->S[setting->tideLocP[ii]] = setting->tidePS;}
    }
    else
    {
      for (ii = 0; ii < setting->tideLocLengthP; ii++)
      {(*data)->S[setting->tideLocP[ii]] = bc->tidalPS[tt];}
    }
  }
  // Direchlet bc at YM tidal boundaries
  if (irank == 0 & setting->bcType != 1)
  {
    if (setting->useConstTideMS == 1)
    {
      for (ii = 0; ii < setting->tideLocLengthM; ii++)
      {(*data)->S[setting->tideLocM[ii]] = setting->tideMS;}
    }
    else
    {
      for (ii = 0; ii < setting->tideLocLengthM; ii++)
      {(*data)->S[setting->tideLocM[ii]] = bc->tidalMS[tt];}
    }
  }
  // Direchlet bc at inflow boundaries
  //for (ii = 0; ii < setting->inflowLocLength; ii++)
  //{(*data)->S[setting->inflowLoc[ii]] = setting->inflowS;}
  // zero scalar over dry lands
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    if ((*data)->depth[ii] <= 0)
    {(*data)->S[ii] = 0;}
  }
}

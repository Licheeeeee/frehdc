#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

// -----------------------------------------------------------------------------
// This file contains the functions to solve the 2D depth integrated Navier
// Stokes equations (including the free surface equations).
// - ZhiLi 2017-05-05 -
// -----------------------------------------------------------------------------

#include "configuration.h"
#include "bathymetry.h"
#include "map.h"
#include "fileio.h"
#include "mpifunctions.h"
#include "nsfunctions.h"
#include "scalar.h"
#include "subgrid.h"
#include "utilities.h"

#include "laspack/errhandl.h"
#include "laspack/vector.h"
#include "laspack/matrix.h"
#include "laspack/qmatrix.h"
#include "laspack/operats.h"
#include "laspack/factor.h"
#include "laspack/precond.h"
#include "laspack/eigenval.h"
#include "laspack/rtc.h"
#include "laspack/itersolv.h"
#include "laspack/mlsolv.h"
#include "laspack/version.h"
#include "laspack/copyrght.h"

void SolveAll(Data **data, Sub **sub, Maps *map, BC *bc, Bath *bath, Config *setting, int irank, int nrank);
void oneCompleteStep(Data **data, Sub **sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, \
  QMatrix A, Vector x, Vector z, int itank, int nrank);
void solveSourceTerm(Data **data, Sub *sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, \
  int irank, int nrank);
void solveFreeSurface(Data **data, Sub *sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, \
  QMatrix A, Vector x, Vector z, int irank, int nrank);
void updateData(Data **data, Sub **sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, \
  int irank, int nrank);

// ========== The top level code of solve ==========
void SolveAll(Data **data, Sub **sub, Maps *map, BC *bc, Bath *bath, Config *setting, int irank, int nrank)
{
  int tt = 0, writeInd, root = 0;
  float t0, t1, tstep;
  // output initial condition
  if (setting->isRestart == 1)
  {writeInd = tt+setting->ttRestart;}
  else
  {writeInd = tt;}
  DataOutput(data, bath, setting, writeInd, root, irank);
  if (setting->useSubgrid == 1) {SubOutput(sub, setting, writeInd, root, irank);}
  // begin time stepping
  if (irank == 0)
  {printf("Ready for time stepping ...\n");}
  for (tt = 1; tt <= setting->Nt; tt++)
  {
    //printf("tide is %lf\n",bc->tideP[tt]);
    if (setting->isRestart == 1)
    {writeInd = tt+setting->ttRestart;}
    else
    {writeInd = tt;}
    QMatrix A;
    Q_Constr(&A, "A", setting->N2ci, False, Rowws, Normal, True);
    Vector z;
    V_Constr(&z, "z", setting->N2ci, Normal, True);
    Vector x;
    V_Constr(&x, "x", setting->N2ci, Normal, True);
    if (irank == 0)
    {t0 = clock();}
    // Solve the Navier Stokes equations for the current time step
    oneCompleteStep(data, sub, map, bc, bath, setting, tt, A, x, z, irank, nrank);
    Q_Destr(&A);
    V_Destr(&x);
    V_Destr(&z);
    // Monitor for stability. Give warning if the maximum CFL number is
    // greater than 2 (which meant possible wetting of 3 cells in 1 step)
    monitorCFL(data, irank, setting, tt);
    // output data
    if (tt % setting->OutItvl == 0)
    {
        DataOutput(data, bath, setting, writeInd, root, irank);
        if (setting->useSubgrid == 1) {SubOutput(sub, setting, writeInd, root, irank);}
    }
    // report time
    if (irank == 0)
    {
      t1 = clock();
      tstep = (t1 - t0)/(float)CLOCKS_PER_SEC;
      printf("---------- Step %d executed with time %.4f sec... ----------\n",writeInd,tstep);
    }
  }
  if (irank == 0)
  {printf("Writing restart file...\n");}
  writeRestartFile(data, bath, setting, writeInd, root, irank);
  if (irank == 0)
  {printf("Time stepping completed...\n");}
}

// ========== Function to solve 1 time step ==========
void oneCompleteStep(Data **data, Sub **sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, \
  QMatrix A, Vector x, Vector z, int irank, int nrank)
{
  solveSourceTerm(data, *sub, map, bc, bath, setting, tt, irank, nrank);
  solveFreeSurface(data, *sub, map, bc, bath, setting, tt, A, x, z, irank, nrank);
  updateData(data, sub, map, bc, bath, setting, tt, irank, nrank);
  if (setting->useScalar == 1)
  {scalarTransport(data, map, bath, bc, *sub, setting, tt, irank, nrank);}
}

// ========== Compute the explicit source term ==========
void solveSourceTerm(Data **data, Sub *sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, int irank, int nrank)
{
  /*if (setting->useSubgrid == 1)
  {computeVolumeSub(data, sub, setting);}
  else
  {computeVolume(data, setting);}*/
  int ii = 33540-400, jj = ii+400, kk = jj + 400, pp = jj-1, qq = jj+1;
  explicitVelocity(data, map, setting);
  //if (irank == 3){
  //printf("EnXP base = %lf,%lf,%lf\n",(*data)->EnXP[pp],(*data)->EnXP[jj],(*data)->EnXP[qq]);
  //printf("EnYP base = %lf,%lf,%lf\n",(*data)->EnYP[ii],(*data)->EnYP[jj],(*data)->EnYP[kk]);}
  if (setting->useSubgrid == 1)
  {
    //advectionTermSub(data, map, sub, setting);
    advectionTerm(data, map, setting);
      //if (irank == 5){
      //printf("EnXP adv = %lf,%lf,%lf\n",(*data)->EnXP[pp],(*data)->EnXP[jj],(*data)->EnXP[qq]);
      //printf("EnYP adv = %lf,%lf,%lf\n",(*data)->EnYP[ii],(*data)->EnYP[jj],(*data)->EnYP[kk]);}
    diffusionTermSub(data, map, sub, setting);
    //windTermSub(data, bc, sub, setting, tt);
    windTerm(data, bc, setting, tt);
      //if (irank == 5){
      //printf("EnXP wind = %lf,%lf,%lf\n",(*data)->EnXP[pp],(*data)->EnXP[jj],(*data)->EnXP[qq]);
      //printf("EnYP wind = %lf,%lf,%lf\n",(*data)->EnYP[ii],(*data)->EnYP[jj],(*data)->EnYP[kk]);}
    //dragTerm(data, map, setting);
    dragTermSub(data, map, sub, setting);
  }
  else
  {
    advectionTerm(data, map, setting);
    //printf("EnXP adv = %lf,%lf,%lf\n",(*data)->EnXP[pp],(*data)->EnXP[jj],(*data)->EnXP[qq]);
    //printf("EnYP adv = %lf,%lf,%lf\n",(*data)->EnYP[ii],(*data)->EnYP[jj],(*data)->EnYP[kk]);
    diffusionTerm(data, map, setting);
    //printf("EnXP diff = %lf,%lf,%lf\n",(*data)->EnXP[pp],(*data)->EnXP[jj],(*data)->EnXP[qq]);
    //printf("EnYP diff = %lf,%lf,%lf\n",(*data)->EnYP[ii],(*data)->EnYP[jj],(*data)->EnYP[kk]);
    windTerm(data, bc, setting, tt);
    //printf("EnXP wind = %lf,%lf,%lf\n",(*data)->EnXP[pp],(*data)->EnXP[jj],(*data)->EnXP[qq]);
    //printf("EnYP wind = %lf,%lf,%lf\n",(*data)->EnYP[ii],(*data)->EnYP[jj],(*data)->EnYP[kk]);
    dragTerm(data, map, setting);
  }
  dragInversion(data, map, setting);
    //if (irank == 3){
  //printf("EnXP drag = %lf,%lf,%lf\n",(*data)->EnXP[pp],(*data)->EnXP[jj],(*data)->EnXP[qq]);
  //printf("EnYP drag = %lf,%lf,%lf\n",(*data)->EnYP[ii],(*data)->EnYP[jj],(*data)->EnYP[kk]);}
  if (setting->useMPI == 1)
  {
    mpiexchange((*data)->EnXP, map, setting, irank, nrank);
    mpiexchange((*data)->EnYP, map, setting, irank, nrank);
  }
  if (setting->useSubgrid != 0)
  {matrixSourceTermSub(data, map, bc, sub, setting, tt, irank);}
  else
  {matrixSourceTerm(data, map, bc, setting, tt, irank);}
  /*int ii,jj;
  for (ii = 38; ii < 48; ii++)
  {
    jj = ii;
    printf("ii,depthXP,depthYP = %d,%lf,%lf\n",ii,(*data)->depthXP[jj],(*data)->depthYP[jj]);
  }*/
}

// ========== Compute the free surface elevation ==========
void solveFreeSurface(Data **data, Sub *sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, \
  QMatrix A, Vector x, Vector z, int irank, int nrank)
{
  if (setting->useSubgrid == 1)
  {matrixCoeffSub(data, map, sub, setting);}
  else
  {matrixCoeff(data, map, setting);}
  adjustMatrixBoundary(data, bc, map, setting, irank, nrank);
  adjustTidalBoundary(data, bc, bath, setting, tt, irank, nrank);
  setupMatrix(*data, map, setting, A);
  solveMatrix(*data, setting, A, x, z);
  getFreeSurface(data, map, setting, x);
  int ii = 2602-60, jj = ii+60, kk = jj + 60, pp = jj-1, qq = jj+1;
    //if (irank == 0){
  //printf("Face areaN = %lf,%lf,%lf\n",sub->Nx[ii],sub->Nx[jj],sub->Nx[kk]);
  //printf("Face areaO = %lf,%lf,%lf\n",sub->Oy[ii],sub->Oy[jj],sub->Oy[kk]);
  //printf("Matrix coeff 1 = %lf,%lf,%lf,%lf,%lf\n",(*data)->GnCt[ii],\
    (*data)->GnXM[ii],(*data)->GnXP[ii],(*data)->GnYM[ii],(*data)->GnYP[ii]);
  //printf("Matrix coeff 2 = %lf,%lf,%lf,%lf,%lf\n",10000*(*data)->GnCt[jj],\
    10000*(*data)->GnXM[jj],10000*(*data)->GnXP[jj],10000*(*data)->GnYM[jj],10000*(*data)->GnYP[jj]);
  //printf("Matrix coeff 3 = %lf,%lf,%lf,%lf,%lf\n",(*data)->GnCt[kk],\
    (*data)->GnXM[kk],(*data)->GnXP[kk],(*data)->GnYM[kk],(*data)->GnYP[kk]);
  //printf("Source term = %lf,%lf,%lf\n",(*data)->z[map->sprt[ii]],(*data)->z[map->sprt[jj]],(*data)->z[map->sprt[kk]]);
  //printf("raw free surface = %lf,%lf,%lf,%lf,%lf\n",(*data)->surf[ii],(*data)->surf[jj],(*data)->surf[kk],(*data)->surf[pp],(*data)->surf[qq]);}

  // compute evaporation and rainfall
  if (setting->useEvap == 1)
  {evaporation(data, bc, tt, setting);}
  if (setting->useRain == 1)
  {rainfall(data, bc, tt, setting);}
  zeroNegSurf(data, bath, setting);
  enforceFreeSurfaceBC(data, map, setting);
  enforceTidalBC(data, bath, bc, setting, tt, irank, nrank);

  //printf("SurfaceX = %lf,%lf,%lf\n",(*data)->surf[pp],(*data)->surf[jj],(*data)->surf[qq]);
  //printf("SurfaceY = %lf,%lf,%lf\n",(*data)->surf[ii],(*data)->surf[jj],(*data)->surf[kk]);
}

// =============== Update depth and velocity ===============
void updateData(Data **data, Sub **sub, Maps *map, BC *bc, Bath *bath, Config *setting, int tt, \
  int irank, int nrank)
{
  int ii = 10999-70, jj = ii+70, kk = jj+70, pp = jj-1, qq = jj+1;
  // ----- remove high velocity or small depth cells -----
  limiterCFL(data, bath, map, setting);
  if (setting->useMPI == 1)
  {
    mpiexchange((*data)->surf, map, setting, irank, nrank);
  }
  // ----- compute grid cell volume and scalar mass -----
  if (setting->useSubgrid != 0)
  {
    computeVolumeSub(data, *sub, setting);
    // here is where the subgrid variables are updated
    computeSubArea(sub, *data, bath, map, setting);
    if (setting->useMPI == 1)
    {
      mpiexchangeInt((*sub)->ind, map, setting, irank, nrank);
      mpiexchange((*sub)->Op, map, setting, irank, nrank);
      mpiexchange((*sub)->Om, map, setting, irank, nrank);
      mpiexchange((*sub)->Vyp, map, setting, irank, nrank);
      mpiexchange((*sub)->Vym, map, setting, irank, nrank);
      mpiexchange((*sub)->Zyp, map, setting, irank, nrank);
      mpiexchange((*sub)->Zym, map, setting, irank, nrank);
    }
    combineSubArea(sub, *data, bath, map, setting, irank, nrank);
    if (setting->useMPI == 1)
    {
      mpiexchange((*sub)->Nx, map, setting, irank, nrank);
      mpiexchange((*sub)->Oy, map, setting, irank, nrank);
      mpiexchange((*sub)->N, map, setting, irank, nrank);
      mpiexchange((*sub)->O, map, setting, irank, nrank);
      mpiexchange((*sub)->V, map, setting, irank, nrank);
      mpiexchange((*sub)->Z, map, setting, irank, nrank);
      mpiexchange((*sub)->Vx, map, setting, irank, nrank);
      mpiexchange((*sub)->Vy, map, setting, irank, nrank);
      mpiexchange((*sub)->Zx, map, setting, irank, nrank);
      mpiexchange((*sub)->Zy, map, setting, irank, nrank);
    }
  }
  else
  {computeVolume(data, setting);}
  // ----- update grid cell depth -----
  updateDepth(data, map, bath, setting);
  updateAllFaceDepth(data, bath, map, setting, irank, nrank);
  if (setting->useMPI == 1)
  {
    mpiexchange((*data)->depthXP, map, setting, irank, nrank);
    mpiexchange((*data)->depthYP, map, setting, irank, nrank);
  }
  // ----- update cell face velocities -----
  if (setting->useSubgrid != 0)
  {
    detectWaterfallLocationSub(data, *sub, bath, map, setting);
      //if (irank == 0){
    //printf("GnXM, GnYM = %lf,%lf\n",(*data)->GnXM[pp],(*data)->GnYM[pp]);
    //printf("free surface = %lf,%lf,%lf,%lf,%lf\n",(*data)->surf[ii],(*data)->surf[jj],(*data)->surf[kk],(*data)->surf[pp],(*data)->surf[qq]);
    //printf("depth = %lf,%lf,%lf,%lf,%lf\n",(*data)->depth[ii],(*data)->depth[jj],(*data)->depth[kk],(*data)->depth[pp],(*data)->depth[qq]);
    //printf("Bathymetry = %lf,%lf,%lf,%lf,%lf\n",bath->bottomZ[ii],bath->bottomZ[jj],bath->bottomZ[kk],bath->bottomZ[pp],bath->bottomZ[qq]);
    //printf("Bathymetry raw = %lf,%lf,%lf,%lf,%lf\n",bath->bottomZ[ii] - bath->offset[0],bath->bottomZ[jj]- bath->offset[0],bath->bottomZ[kk]- bath->offset[0],bath->bottomZ[pp]- bath->offset[0],bath->bottomZ[qq]- bath->offset[0]);
    //printf("Edge depth = %lf,%lf,%lf\n",(*data)->depthYP[ii],(*data)->depthYP[jj],(*data)->depthYP[kk]);
    //printf("Edge bath = %lf,%lf\n",bath->bottomXP[ii],bath->bottomXP[jj]);
    //printf("New Face areaN = %lf,%lf,%lf\n",(*sub)->Nx[pp],(*sub)->Nx[jj],(*sub)->Nx[qq]);
    //printf("New Face areaO = %lf,%lf,%lf\n",(*sub)->Oy[ii],(*sub)->Oy[jj],(*sub)->Oy[kk]);
    //printf("New Face volumeN = %lf,%lf,%lf\n",(*sub)->Vx[pp],(*sub)->Vx[jj],(*sub)->Vx[qq]);
    //printf("New Face volumeO = %lf,%lf,%lf\n",(*sub)->Vy[ii],(*sub)->Vy[jj],(*sub)->Vy[kk]);}
    updateVelocitySub(data, map, *sub, setting);
    //printf("VelocityU = %lf,%lf,%lf\n",(*data)->uuXP[ii],(*data)->uuXP[jj],(*data)->uuXP[kk]);
    //printf("EnYP = %lf,%lf,%lf\n",(*data)->EnYP[ii],(*data)->EnYP[jj],(*data)->EnYP[kk]);
    //printf("VelocityV = %lf,%lf,%lf\n",(*data)->vvYP[ii],(*data)->vvYP[jj],(*data)->vvYP[pp]);
    waterfallCorrectionSub(data, *sub, bath, map, setting);
    //printf("Waterfall VelocityU = %lf,%lf,%lf\n",(*data)->uuXP[ii],(*data)->uuXP[jj],(*data)->uuXP[kk]);
    //printf("Waterfall VelocityV = %lf,%lf,%lf\n",(*data)->vvYP[ii],(*data)->vvYP[jj],(*data)->vvYP[pp]);
  }
  else
  {
    detectWaterfallLocation(data, bath, map, setting);
    updateVelocity(data, map, setting);
    waterfallCorrection(data, bath, map, setting);
  }
  enforceVelocityBC(data, map, setting);
  if (setting->useMPI == 1)
  {
    mpiexchange((*data)->uuXP, map, setting, irank, nrank);
    mpiexchange((*data)->vvYP, map, setting, irank, nrank);
  }
  // ----- update cell face flow rates -----
  if (setting->useSubgrid != 0)
  {
    computeFaceFlowRateSub(data, map, *sub, setting);
    adjustTidalVelocitySub(data, map, *sub, setting, irank, nrank);
    //velocityInterpSub(data, *sub, setting);
    velocityInterp(data, setting);
    volumeByFluxSub(data, *sub, map, bc, setting, tt, irank);
  }
  else
  {
    computeFaceFlowRate(data, map, setting);
    adjustTidalVelocity(data, map, setting, irank, nrank);
    velocityInterp(data, setting);
    volumeByFlux(data, map, bc, setting, tt, irank);
  }
  if (setting->useMPI == 1)
  {
    mpiexchange((*data)->Fuu, map, setting, irank, nrank);
    mpiexchange((*data)->Fvv, map, setting, irank, nrank);
  }

  //if (irank == 3){
  //printf("EnXP = %lf,%lf,%lf\n",(*data)->EnXP[pp],(*data)->EnXP[jj],(*data)->EnXP[qq]);
  //printf("EnYP = %lf,%lf,%lf\n",(*data)->EnYP[ii],(*data)->EnYP[jj],(*data)->EnYP[kk]);
  //printf("VelocityU = %lf,%lf,%lf\n",(*data)->uuXP[pp],(*data)->uuXP[jj],(*data)->uuXP[qq]);
  //printf("VelocityV = %lf,%lf,%lf\n",(*data)->vvYP[ii],(*data)->vvYP[jj],(*data)->vvYP[kk]);
  //printf("Flux uu = %lf,%lf,%lf\n",(*data)->Fuu[pp],(*data)->Fuu[jj],(*data)->Fuu[qq]);
  //printf("Flux vv = %lf,%lf,%lf\n",(*data)->Fvv[ii],(*data)->Fvv[jj],(*data)->Fvv[kk]);}
  //printf("Volume by flux = %lf,%lf,%lf\n",(*data)->cellV[ii],(*data)->cellV[jj],(*data)->cellV[jj+1]);
  //printf("Subgrid volume = %lf,%lf,%lf\n",(*sub)->V[ii],(*sub)->V[jj],(*sub)->V[jj+1]);
  if (setting->useSubgrid != 0 & setting->useSubDrag != 0)
  {
    //printf("bath offset = %lf\n",bath->offset[0]);
    //printf("free surface = %lf,%lf,%lf,%lf,%lf\n",(*data)->surf[ii],(*data)->surf[jj],(*data)->surf[kk],(*data)->surf[pp],(*data)->surf[qq]);
    //printf("Yh = %lf,%lf,%lf,%lf,%lf\n",(*sub)->Yh[ii],(*sub)->Yh[jj],(*sub)->Yh[kk],(*sub)->Yh[pp],(*sub)->Yh[qq]);
    updateCDSub(data, *sub, setting);
    //printf("CDXP = %lf,%lf,%lf,%lf,%lf\n",(*data)->CDXP[ii],(*data)->CDXP[jj],(*data)->CDXP[kk],(*data)->CDXP[pp],(*data)->CDXP[qq]);
  }
  else
  {
    updateCD(data, setting);
  }
  thinLayerDrag(data, setting);
  // printf("CDXP = %lf,%lf,%lf,%lf,%lf\n",(*data)->CDXP[ii],(*data)->CDXP[jj],(*data)->CDXP[kk],(*data)->CDXP[pp],(*data)->CDXP[qq]);
}

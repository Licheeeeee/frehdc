#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

// -----------------------------------------------------------------------------
// Functions defining the grid maps of the model
// - 2017-05-04 by Zhi Li -
// -----------------------------------------------------------------------------

#include "configuration.h"
#include "map.h"

void createMaps(Maps **map, Config *setting)
{
  int ii, jj, kk, col;
  *map = malloc(sizeof(Maps));
  (*map)->cntr = malloc(setting->N2ci*sizeof(int));
  (*map)->trps = malloc(setting->N2ci*sizeof(int));
  (*map)->sprt = malloc(setting->N2ci*sizeof(int));
  (*map)->iPjc = malloc(setting->N2ci*sizeof(int));
  (*map)->iMjc = malloc(setting->N2ci*sizeof(int));
  (*map)->icjP = malloc(setting->N2ci*sizeof(int));
  (*map)->icjM = malloc(setting->N2ci*sizeof(int));
  (*map)->iPjP = malloc(1*sizeof(int));
  (*map)->iPjM = malloc(1*sizeof(int));
  (*map)->iMjP = malloc(1*sizeof(int));
  (*map)->iMjM = malloc(1*sizeof(int));
  (*map)->ii2d = malloc(setting->N2ci*sizeof(int));
  (*map)->jj2d = malloc(setting->N2ci*sizeof(int));
  (*map)->iPbd = malloc(setting->ny*sizeof(int));
  (*map)->iPgt = malloc(setting->ny*sizeof(int));
  (*map)->iMbd = malloc(setting->ny*sizeof(int));
  (*map)->iMgt = malloc(setting->ny*sizeof(int));
  (*map)->jPbd = malloc(setting->nx*sizeof(int));
  (*map)->jPgt = malloc(setting->nx*sizeof(int));
  (*map)->jMbd = malloc(setting->nx*sizeof(int));
  (*map)->jMgt = malloc(setting->nx*sizeof(int));
  // set index for the center map
  for (ii = 0; ii < setting->N2ci; ii++)
  {(*map)->cntr[ii] = ii;}
  // set index for the transposed map
  for (kk = 0; kk < setting->N2ci; kk++)
  {
    ii = floor(kk / setting->ny);
    jj = kk % setting->ny;
    (*map)->trps[kk] = jj*setting->nx + ii;
    ii = floor(kk / setting->nx);
    jj = kk % setting->nx;
    (*map)->sprt[kk] = jj*setting->ny + ii;
  }
  // set index for the jP and jM maps
  for (ii = 0; ii < setting->N2ci-setting->nx; ii++)
  {
    (*map)->icjP[ii] = ii + setting->nx;
    (*map)->icjM[ii+setting->nx] = ii;
  }
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*map)->icjP[ii+setting->N2ci-setting->nx] = ii + setting->N2ci;
    (*map)->icjM[ii] = ii + setting->N2ci + setting->nx;
  }
  // set index for the iP and iM maps
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    col = floor(ii/setting->nx);
    if (ii % setting->nx == setting->nx-1)
    {(*map)->iPjc[ii] = setting->N2ci + 2*setting->nx + col; (*map)->iMjc[ii] = ii - 1;}
    else if (ii % setting->nx == 0)
    {(*map)->iMjc[ii] = setting->N2ci + 2*setting->nx + setting->ny + col; (*map)->iPjc[ii] = ii + 1;}
    else
    {(*map)->iPjc[ii] = ii + 1;   (*map)->iMjc[ii] = ii - 1;}
  }
  // set corner maps
  (*map)->iMjM[0] = setting->N2ct - 1;
  (*map)->iPjM[0] = setting->N2ct - 2;
  (*map)->iPjP[0] = setting->N2ct - 3;
  (*map)->iMjP[0] = setting->N2ct - 4;
  // set the 2d (ii,jj) maps
  for (ii = 0; ii < setting->N2ci; ii++)
  {
    col = floor(ii/setting->nx);
    (*map)->ii2d[ii] = ii - col*setting->nx;
    (*map)->jj2d[ii] = col;
  }
  // set the maps for the jP and jM ghost cells
  for (ii = 0; ii < setting->nx; ii++)
  {
    (*map)->jPbd[ii] = (*map)->cntr[setting->N2ci-setting->nx+ii];
    (*map)->jMbd[ii] = (*map)->cntr[ii];
    (*map)->jPgt[ii] = (*map)->icjP[setting->N2ci-setting->nx+ii];
    (*map)->jMgt[ii] = (*map)->icjM[ii];
  }
  // set the maps for the iP and iM ghost cells
  for (ii = 0; ii < setting->ny; ii++)
  {
    (*map)->iPbd[ii] = (*map)->cntr[((ii+1)*setting->nx)-1];
    (*map)->iMbd[ii] = (*map)->cntr[ii*setting->nx];
    (*map)->iPgt[ii] = (*map)->iPjc[((ii+1)*setting->nx)-1];
    (*map)->iMgt[ii] = (*map)->iMjc[ii*setting->nx];
  }
}

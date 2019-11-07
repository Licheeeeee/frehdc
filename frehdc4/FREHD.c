#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<mpi.h>

// -----------------------------------------------------------------------------
// This is the top level script of the SubFREHD-C model.
// - Zhi Li 2017-05-10 -
// -----------------------------------------------------------------------------

#include "bathymetry.h"
#include "configuration.h"
#include "fileio.h"
#include "initialize.h"
#include "map.h"
#include "mpifunctions.h"
#include "nsfunctions.h"
#include "nssolve.h"
#include "subgrid.h"
#include "utilities.h"

int main(int argc, char *argv[])
{
  Config *setting;
  Bath *bath;
  Maps *map;
  Data *data;
  IC *ic;
  BC *bc;
  Sub *sub;
  int irank = 0, nrank = 1;
  double ts, te;

  // timer
  if (irank == 0) {ts = clock();}
  // read user settings
  writeText("----- STATUS ----- Reading user settings...", irank);
  ReadUserSettings(&setting);
  // start MPI
  if (setting->useMPI == 1)
  {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
  }
  // read bathymetry
  writeText("----- STATUS ----- Reading bathymetry...", irank);
  ReadBathymetry(&bath, setting, irank);
  InitBathymetry(bath, setting, irank);
  writeText("----- STATUS ----- Creating maps...", irank);
  createMaps(&map, setting);
  // Initialize the model
  writeText("----- STATUS ----- Initializing...", irank);
  Init(&data, &map, &ic, &bc, &sub, bath, setting, irank, nrank);
  // Solve the Navier Stokes equations
  writeText("----- STATUS ----- Begin to solve...", irank);
  SolveAll(&data, &sub, map, bc, bath, setting, irank, nrank);
  // end MPI
  if (setting->useMPI == 1)
  { MPI_Finalize();}
  writeText("----- STATUS ----- Completing...", irank);
  // timer
  if (irank == 0) {te = clock(); printf("Total time = %lf min\n",(te-ts)/60/(float)CLOCKS_PER_SEC);}
  return 0;
}

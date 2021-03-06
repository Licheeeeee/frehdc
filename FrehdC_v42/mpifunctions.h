

#include "configuration.h"
#include "map.h"

void combineAllRanks(double *ally, double *y, Config *setting, int root);
void mpiexchange(double *S, Maps *map, Config *setting, int irank, int nrank);
void mpiexchangeInt(int *S, Maps *map, Config *setting, int irank, int nrank);

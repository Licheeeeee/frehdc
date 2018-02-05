

#ifndef MAP_H
#define MAP_H

typedef struct Maps
{
  int *cntr, *trps, *sprt, *iPjc, *iMjc, *icjP, *icjM, *ii2d, *jj2d;
  int *iPjP, *iPjM, *iMjP, *iMjM;
  int *iPbd, *iPgt, *iMbd, *iMgt, *jPbd, *jPgt, *jMbd, *jMgt;
}Maps;

#endif

void createMaps(Maps **map, Config *setting);

/*
 * Demonstrates how to use staggered grids with boundary conditions
 */

#include <bout.hxx>
#include <boutmain.hxx>

#include <derivs.hxx>
#include <interpolation.hxx>

Field3D n, v, n_stag;

int physics_init(bool restart) {
  
  v.setLocation(CELL_YLOW); // Staggered relative to n
  n_stag.setLocation(CELL_YLOW); // Staggered relative to n
  
  SOLVE_FOR2(n, v);

  return 0;
}

int physics_run(BoutReal time) {
  mesh->communicate(n, v);
  n_stag = interp_to(n,CELL_YLOW);

  //ddt(n) = -Div_par_flux(v, n, CELL_CENTRE);
  ddt(n) = -n*Grad_par(v, CELL_CENTRE) - Vpar_Grad_par(v, n, CELL_CENTRE);
  
  //ddt(v) = -Grad_par(n, CELL_YLOW);
  ddt(v) = -Grad_par(n_stag);
 
  // Have to manually apply the lower Y boundary region, using a width of 3
  for( RangeIterator rlow = mesh->iterateBndryLowerY(); !rlow.isDone(); rlow++)
    for(int y=2;y>=0;y--) 
      for(int z=0;z<mesh->LocalNz;z++) {
        ddt(v)(rlow.ind,y,z) = ddt(v)(rlow.ind,y+1,z);
      }
  
  return 0;
}


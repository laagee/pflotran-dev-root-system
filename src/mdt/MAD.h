#ifndef MAD_H_
#define MAD_H_

#include <string.h>

#include "petscsys.h"

#include "Globals.h"
#include "AsciiGrid.h"
#include "BoundarySet.h"
#include "Connection.h"
#include "Grid.h"
#include "Polygon.h"

class MAD {
  
public:
  MAD(Grid **grid);
  virtual ~MAD();

  void computeTopBoundary(Grid *grid, PetscInt complete);
  void computeBottomBoundary(Grid *grid, PetscInt complete);
  void computeNorthBoundary(Grid *grid, PetscInt complete);
  void computeSouthBoundary(Grid *grid, PetscInt complete);
  void computeEastBoundary(Grid *grid, PetscInt complete);
  void computeWestBoundary(Grid *grid, PetscInt complete);
  void computeIFCBoundary(Grid *grid, Polygon *p);
  void computeUnSatSPPDomain(Grid *grid, Polygon *p);
  void computeSatSPPDomain(Grid *grid, Polygon *p);
  void flagGridCells(Grid *grid);

private:

  Polygon *ifc_polygon;
  Polygon *spp_polygon;
  Polygon *river_polygon;
  AsciiGrid **ascii_grids;

  void setMaterialIdBasedOnNaturalId(PetscInt natural_id, PetscInt material_id,
                                     Grid *grid);
  void setActiveBasedOnNaturalId(PetscInt natural_id, PetscInt active,
                                 Grid *grid);

};

#endif /*MAD_H_*/

#include "IFC.h"

#ifndef MAX
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#endif
#ifndef MIN
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#endif

IFC::IFC(Grid **grid_) {

  ifc_polygon = NULL;
  river_polygon = NULL;
  ascii_grids = NULL;
  spp_polygon = NULL;

  PetscReal mx = 1.;
  PetscReal my = 1.;
  PetscReal mz = 1.;

  PetscInt nx, ny, nz;

  char filename[1024];
  PetscBool option_found;
  strcpy(filename,"mdt.in");
  PetscOptionsGetString(PETSC_NULL,"-mdtin",filename,1024,&option_found);

  FileIO *file = new FileIO(filename);
  file->getLine();
  file->readInt(&nx);
  file->readInt(&ny);
  file->readInt(&nz);
  delete file;

  PetscReal dx;
  PetscReal dy;
  PetscReal dz;// */

#if 0
  PetscReal len_x = 500.;
  PetscReal len_y = 700.;
#else
  PetscReal len_x = 850.;
  PetscReal len_y = 1000.;
#endif
  PetscReal len_z = 20.;

  dx = len_x/(PetscReal)nx;
  dy = len_y/(PetscReal)ny;
  dz = len_z/(PetscReal)nz;

  PetscInt n = nx*ny*nz;

  PetscPrintf(PETSC_COMM_WORLD,"nx = %d, dx = %f, lenx = %f\n",nx,dx,nx*dx);
  PetscPrintf(PETSC_COMM_WORLD,"ny = %d, dy = %f, leny = %f\n",ny,dy,ny*dy);
  PetscPrintf(PETSC_COMM_WORLD,"nz = %d, dz = %f, lenz = %f\n",nz,dz,nz*dz);
  *grid_ = new Grid(nx,ny,nz);
  Grid *grid = *grid_;

// grid spacing with a bias
#if 0
  PetscReal sum_x = 0.;
  PetscReal sum_y = 0.;
  dx = 0.8470329472543
  PetscReal *dx_array = new double[nx];
  PetscReal *dy_array = new double[ny];
  PetscReal *dz_array = new double[nz];
  for (int i=0; i<nx; i++)
    dx_array[i] = 10.;
  for (int i=0; i<ny; i++)
    dy_array[i] = 10.;
  for (int i=0; i<nz; i++)
    dz_array[i] = 0.25;

  for (int i=11; i<19; i++) {
    dx_array[i] = 10.*pow(1.30242241518419,(double)(10-i));
    sum_x += dx_array[i];
  }

  for (int i=19; i<89; i++) {
    dx_array[i] = 1.;
    sum_x += dx_array[i];
  }

  for (int i=89; i<97; i++) {
    dx_array[i] = 10.*pow(1.30242241518419,(double)(i-97));
    sum_x += dx_array[i];
  }

  for (int i=97; i<9; i++) {
    dy_array[110+i] = 10.*pow(1.353088,i+1.);
    dy_array[9-i] = 10.*pow(1.353088,i+1.);
    sum_y += dy_array[9-i];
  }
  grid->setGridSpacing(dx_array,dy_array,dz_array);
#else
  grid->setGridSpacing(dx,dy,dz);
#endif

  grid->setRotation(14.); // must come before ->setOrigin()

#if 0
  grid->setOrigin(594200.,115550.,90.);
#else
  grid->setOrigin(593875.,115360.,90.);
#endif

//  grid->computeCoordinates();
//  grid->computeConnectivity();
  grid->computeCellMapping();
  grid->setUpCells();
  grid->computeVertexMapping();
  grid->setUpVertices();
  grid->mapVerticesToCells();

  // don't know that we need this since we can use conductance of river.
//  river_polygon = new Polygon();
//  river_polygon->createRiverEdgePolygon();

  ifc_polygon = new Polygon();
  ifc_polygon->createIFCPolygon();
  spp_polygon = new Polygon();
  spp_polygon->createSPPPolygon();

#if 1
#if 0
  char ascii_filename[1024];
  strcpy(ascii_filename,"test.asc");
  AsciiGrid **ascii_grid = new AsciiGrid*[2];
  ascii_grid[0] = new AsciiGrid(ascii_filename);
  ascii_grid[0]->setMaterialId(1);
  ascii_grid[1] = new AsciiGrid("default",2,2,-1.e10,-1.e10,1.e10,-9999.,
                                1.e20,0); 
#else

  AsciiGrid::nasciigrids = 6;
  string *grid_filenames = new string[AsciiGrid::nasciigrids];
#if 1
  grid_filenames[0].append("./basalt_300area.asc");
  grid_filenames[1].append("./u9_300area.asc");
  grid_filenames[2].append("./u8_300area.asc");
  grid_filenames[3].append("./u5gravel_300area.asc");
  grid_filenames[4].append("./u5silt_300area.asc");
  grid_filenames[5].append("./newbath_10mDEM_grid.ascii");
#else
  grid_filenames[0].append("../basalt_300area.asc");
  grid_filenames[1].append("../u9_300area.asc");
  grid_filenames[2].append("../u8_300area.asc");
  grid_filenames[3].append("../u5gravel_300area.asc");
  grid_filenames[4].append("../u5silt_300area.asc");
  grid_filenames[5].append("../newbath_10mDEM_grid.ascii");
#endif

  ascii_grids = new AsciiGrid*[AsciiGrid::nasciigrids];
  for (PetscInt i=0; i<AsciiGrid::nasciigrids; i++) {
    char filename[32];
    strcpy(filename,grid_filenames[i].c_str());
    ascii_grids[i] = new AsciiGrid(filename);
  }
  ascii_grids[0]->setMaterialId(10);
  ascii_grids[1]->setMaterialId(9);
  ascii_grids[2]->setMaterialId(8);
  ascii_grids[3]->setMaterialId(6);
  ascii_grids[4]->setMaterialId(5);
  ascii_grids[5]->setMaterialId(1);

  PetscInt mod = grid->num_cells_ghosted/10;
  for (PetscInt i=0; i<grid->num_cells_ghosted; i++) {
    PetscInt material_id = 0;
    PetscReal x = grid->cells[i].getX();
    PetscReal y = grid->cells[i].getY();
    PetscReal z = grid->cells[i].getZ();
    for (PetscInt ilayer=0; ilayer<AsciiGrid::nasciigrids; ilayer++) {
      PetscReal zlayer = ascii_grids[ilayer]->computeElevationFromCoordinate(x,y);
      if (zlayer > ascii_grids[ilayer]->nodata && zlayer >= z) {
        material_id = ascii_grids[ilayer]->getMaterialId();
        break;
      }
    }
    if (material_id == 0) grid->cells[i].setActive(0);
    grid->cells[i].setMaterialId(material_id);
    if (river_polygon) {
      if (!river_polygon->pointInPolygon(x,y)) {
        grid->cells[i].setActive(0);
        grid->cells[i].negateMaterialId();
      }
    }
    if (i%mod == 0) {
      PetscPrintf(PETSC_COMM_WORLD,"%d of %d cells mapped with materials and activity.\n",
                  i,grid->num_cells_ghosted);
    }
  }

#endif
#endif

  flagGridCells(grid);

//  setEastBoundaryMaterialTo2(grid);

  computeEastBoundary(grid,1);
  computeWestBoundary(grid,1);
  computeNorthBoundary(grid,0);
  computeSouthBoundary(grid,0);
  computeTopBoundary(grid,0);

  computeIFCBoundary(grid,ifc_polygon);
  computeUnSatSPPDomain(grid,spp_polygon);
  computeSatSPPDomain(grid,spp_polygon);
  computeTopSPPDomain(grid,spp_polygon);

  BoundarySet *river = grid->getBoundarySet("East");
  BoundarySet *west = grid->getBoundarySet("West");
  BoundarySet *north = grid->getBoundarySet("North");
  BoundarySet *south = grid->getBoundarySet("South");
  BoundarySet *recharge = grid->getBoundarySet("Top");

/*
  Condition *new_condition = new Condition("river.bc");
  river->condition = new_condition;
  new_condition = new Condition("west.bc");
  west->condition = new_condition;
  new_condition = new Condition("recharge.bc");
  recharge->condition = new_condition;
  new_condition = NULL;
*/  

}

void IFC::setEastBoundaryMaterialTo2(Grid *grid) {

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].getMaterialId() == 1 &&
        (grid->cells[i].flag & EAST_DIR_EAST_FACE || 
         grid->cells[i].flag & EAST_DIR_SOUTH_FACE || 
         grid->cells[i].flag & EAST_DIR_NORTH_FACE || 
         grid->cells[i].flag & EAST_DIR_BOTTOM_FACE || 
         grid->cells[i].flag & EAST_DIR_TOP_FACE)) {
        grid->cells[i].setMaterialId(2);
        if (grid->cells[i-1].getActive() &&
            grid->cells[i-1].getMaterialId() == 1) grid->cells[i-1].setMaterialId(2);
      }
    }
  }

}

void IFC::computeWestBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *west = new BoundarySet("West");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & WEST_DIR_WEST_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(WEST,vertex_list);
        west->addConnection(new Connection(local_id,vertex_list,WEST));
      }
      if (complete) {
        if (grid->cells[i].flag & WEST_DIR_SOUTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list,SOUTH));
        }
        if (grid->cells[i].flag & WEST_DIR_NORTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list,NORTH));
        }
        if (grid->cells[i].flag & WEST_DIR_BOTTOM_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list,BOTTOM));
        }
        if (grid->cells[i].flag & WEST_DIR_TOP_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          west->addConnection(new Connection(local_id,vertex_list,TOP));
        }
      }
    }
  }

  grid->addBoundarySet(west);
  west = NULL;

}

void IFC::computeEastBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *east = new BoundarySet("East");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & EAST_DIR_EAST_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(EAST,vertex_list);
        east->addConnection(new Connection(local_id,vertex_list,EAST));
      }
      if (complete) {
        if (grid->cells[i].flag & EAST_DIR_SOUTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list,SOUTH));
        }
        if (grid->cells[i].flag & EAST_DIR_NORTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list,NORTH));
        }
        if (grid->cells[i].flag & EAST_DIR_BOTTOM_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list,BOTTOM));
        }
        if (grid->cells[i].flag & EAST_DIR_TOP_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          east->addConnection(new Connection(local_id,vertex_list,TOP));
        }
      }
    }
  }

  grid->addBoundarySet(east);
  east = NULL;

}

void IFC::computeSouthBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *south = new BoundarySet("South");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & SOUTH_DIR_SOUTH_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
        south->addConnection(new Connection(local_id,vertex_list,SOUTH));
      }
      if (complete) {
        if (grid->cells[i].flag & SOUTH_DIR_WEST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(WEST,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list,WEST));
        }
        if (grid->cells[i].flag & SOUTH_DIR_EAST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(EAST,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list,EAST));
        }
        if (grid->cells[i].flag & SOUTH_DIR_BOTTOM_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list,BOTTOM));
        }
        if (grid->cells[i].flag & SOUTH_DIR_TOP_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          south->addConnection(new Connection(local_id,vertex_list,TOP));
        }
      }
    }
  }

  grid->addBoundarySet(south);
  south = NULL;
}

void IFC::computeNorthBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *north = new BoundarySet("North");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & NORTH_DIR_NORTH_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
        north->addConnection(new Connection(local_id,vertex_list,NORTH));
      }
      if (complete) {
        if (grid->cells[i].flag & NORTH_DIR_WEST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(WEST,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list,WEST));
        }
        if (grid->cells[i].flag & NORTH_DIR_EAST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(EAST,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list,EAST));
        }
        if (grid->cells[i].flag & NORTH_DIR_BOTTOM_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list,BOTTOM));
        }
        if (grid->cells[i].flag & NORTH_DIR_TOP_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(TOP,vertex_list);
          north->addConnection(new Connection(local_id,vertex_list,TOP));
        }
      }
    }
  }

  grid->addBoundarySet(north);
  north = NULL;

}

void IFC::computeBottomBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *bottom = new BoundarySet("Bottom");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & BOTTOM_DIR_BOTTOM_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(BOTTOM,vertex_list);
        bottom->addConnection(new Connection(local_id,vertex_list,BOTTOM));
      }
      if (complete) {
        if (grid->cells[i].flag & BOTTOM_DIR_WEST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(WEST,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list,WEST));
        }
        if (grid->cells[i].flag & BOTTOM_DIR_EAST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(EAST,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list,EAST));
        }
        if (grid->cells[i].flag & BOTTOM_DIR_SOUTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list,SOUTH));
        }
        if (grid->cells[i].flag & BOTTOM_DIR_NORTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
          bottom->addConnection(new Connection(local_id,vertex_list,NORTH));
        }
      }
    }
  }

  grid->addBoundarySet(bottom);
  bottom = NULL;

}

void IFC::computeTopBoundary(Grid *grid, PetscInt complete) {

  BoundarySet *top = new BoundarySet("Top");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & TOP_DIR_TOP_FACE) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(TOP,vertex_list);
        top->addConnection(new Connection(local_id,vertex_list,TOP));
      }
      if (complete) {
#if 0  
! not enough space to define these
        if (grid->cells[i].flag & TOP_DIR_WEST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(WEST,vertex_list);
          top->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & TOP_DIR_EAST_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(EAST,vertex_list);
          top->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & TOP_DIR_SOUTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(SOUTH,vertex_list);
          top->addConnection(new Connection(local_id,vertex_list));
        }
        if (grid->cells[i].flag & TOP_DIR_NORTH_FACE) {
          PetscInt vertex_list[5] = {4,0,0,0,0};
          grid->cells[i].getHexFaceVertices(NORTH,vertex_list);
          top->addConnection(new Connection(local_id,vertex_list));
        }
#endif
      }
    }
  }

  grid->addBoundarySet(top);
  top = NULL;

}

void IFC::computeIFCBoundary(Grid *grid, Polygon *p) {

  BoundarySet *plume = new BoundarySet("IFC_Boundary");

  PetscInt count = 0;
  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
//      if (grid->cells[i].flag & TOP_DIR_TOP_FACE &&
      if (grid->cells[i].getZ() >= 105. && grid->cells[i].getZ() <= 108.5 &&
          p->pointInPolygon(grid->cells[i].getX(),
                            grid->cells[i].getY())) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(TOP,vertex_list);
        plume->addConnection(new Connection(local_id,vertex_list,TOP));
        count++;
      }
    }
  }

  PetscPrintf(PETSC_COMM_WORLD,"%d cells mapped to Plume Surface.\n",
                  count);

  grid->addBoundarySet(plume);
  plume = NULL;

}

void IFC::computeTopSPPDomain(Grid *grid, Polygon *p) {

  BoundarySet *top = new BoundarySet("SPP_Top");

  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
      if (grid->cells[i].flag & TOP_DIR_TOP_FACE &&
          p->pointInPolygon(grid->cells[i].getX(),
                            grid->cells[i].getY())) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(TOP,vertex_list);
        top->addConnection(new Connection(local_id,vertex_list,TOP));
      }
    }
  }

  grid->addBoundarySet(top);
  top = NULL;

}

void IFC::computeSatSPPDomain(Grid *grid, Polygon *p) {

  BoundarySet *spp = new BoundarySet("SPP_Saturated");

  PetscInt count = 0;
  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
//      if (grid->cells[i].flag & TOP_DIR_TOP_FACE &&
      if (grid->cells[i].getZ() >= 103. && grid->cells[i].getZ() <= 105. &&
          p->pointInPolygon(grid->cells[i].getX(),
                            grid->cells[i].getY())) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(TOP,vertex_list);
        spp->addConnection(new Connection(local_id,vertex_list,TOP));
        count++;
      }
    }
  }

  PetscPrintf(PETSC_COMM_WORLD,"%d cells mapped to SPP_Saturated Domain.\n",
                  count);

  grid->addBoundarySet(spp);
  spp = NULL;

}

void IFC::computeUnSatSPPDomain(Grid *grid, Polygon *p) {

  BoundarySet *spp = new BoundarySet("SPP_Unsaturated");

  PetscInt count = 0;
  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) {
    PetscInt local_id = grid->cells[i].getIdLocal();
    if (local_id > -1) {
//      if (grid->cells[i].flag & TOP_DIR_TOP_FACE &&
      if (grid->cells[i].getZ() > 105. &&
          p->pointInPolygon(grid->cells[i].getX(),
                            grid->cells[i].getY())) {
        PetscInt vertex_list[5] = {4,0,0,0,0};
        grid->cells[i].getHexFaceVertices(TOP,vertex_list);
        spp->addConnection(new Connection(local_id,vertex_list,TOP));
        count++;
      }
    }
  }

  PetscPrintf(PETSC_COMM_WORLD,"%d cells mapped to SPP_Unsaturated Domain.\n",
                  count);

  grid->addBoundarySet(spp);
  spp = NULL;

}

void IFC::flagGridCells(Grid *grid) {

  PetscReal top_stage_max = 1.e20;
  PetscReal top_stage_min = 108.;
  PetscReal bottom_stage_max = 1.e20;
  PetscReal bottom_stage_min = -1.e20;
  PetscReal north_stage_max = 1.e20;
  PetscReal north_stage_min = -1.e20;
  PetscReal south_stage_max = 1.e20;
  PetscReal south_stage_min = -1.e20;
  PetscReal west_stage_max = 1.e20;
  PetscReal west_stage_min = -1.e20;
  PetscReal east_stage_max = 107.5;
  PetscReal east_stage_min = -1.e20;

  PetscInt nx = grid->getNx();
  PetscInt ny = grid->getNy();
  PetscInt nz = grid->getNz();

  PetscInt lnx,lny,lnz,lxs,lys,lzs,lxe,lye,lze;
  PetscInt gnx,gny,gnz,gxs,gys,gzs,gxe,gye,gze;

  grid->getCorners(&lxs,&lys,&lzs,&lnx,&lny,&lnz);
  grid->getGhostCorners(&gxs,&gys,&gzs,&gnx,&gny,&gnz);

  lxe = lxs+lnx;
  lye = lys+lny;
  lze = lzs+lnz;
  gxe = gxs+gnx;
  gye = gys+gny;
  gze = gzs+gnz;

  PetscInt gnxXny = gnx*gny;

  PetscInt istart = lxs-gxs;
  PetscInt jstart = lys-gys;
  PetscInt kstart = lzs-gzs;
  PetscInt iend = istart+lnx;
  PetscInt jend = jstart+lny;
  PetscInt kend = kstart+lnz;


  grid->zeroGridCellFlags();

  PetscInt *matrix = NULL;

  // west
  matrix = new PetscInt[ny*nz];
  for (PetscInt i=0; i<ny*nz; i++)
    matrix[i] = -1;

  for (PetscInt kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (PetscInt jglobal=0; jglobal<ny; jglobal++) {
        if (lys <= jglobal && jglobal < lye) {
          PetscInt j = jglobal-gys;
          PetscInt k = kglobal-gzs;

          PetscInt flag[3];
          flag[0] = 0;
          flag[1] = 0;
          flag[2] = 0;
          if (lxs > 0) grid->receiveFlag(flag,WEST);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt i=0; i<gnx; i++) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  grid->cells[ghosted_id].getZ() <= west_stage_max &&
                  grid->cells[ghosted_id].getZ() >= west_stage_min) {
                grid->cells[ghosted_id].flag |= WEST_DIR_WEST_FACE;
                flag[0] = 1;
                matrix[jglobal+kglobal*ny] = i+gxs;
                break;
              }
            }
          }
          if (lxe < nx) grid->sendFlag(flag,EAST);
        }
      }
    }
  }
  PetscMPIInt nyXnz = ny*nz;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nyXnz,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);

  for (PetscInt k=0; k<gnz; k++) {
    PetscInt kglobal = k+gzs;
    for (PetscInt j=0; j<gny; j++) {
      PetscInt jglobal = j+gys;
      // y-direction
      if (j < gny-1) {
        PetscInt iup = (PetscInt)matrix[jglobal+kglobal*ny];
        PetscInt idown = (PetscInt)matrix[jglobal+1+kglobal*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (PetscInt i=idown; i<iup; i++)
              grid->cells[i+(j+1)*gnx+k*gnxXny].flag |= WEST_DIR_SOUTH_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (PetscInt i=iup; i<idown; i++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= WEST_DIR_NORTH_FACE;     
          }
        }
      }
      // z-direction
      if (k < gnz-1) {
        PetscInt iup = (PetscInt)matrix[jglobal+kglobal*ny];
        PetscInt idown = (PetscInt)matrix[jglobal+(kglobal+1)*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (PetscInt i=idown; i<iup; i++)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= WEST_DIR_BOTTOM_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (PetscInt i=iup; i<idown; i++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= WEST_DIR_TOP_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // east
  matrix = new PetscInt[ny*nz];
  for (PetscInt i=0; i<ny*nz; i++)
    matrix[i] = -1;

  for (PetscInt kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (PetscInt jglobal=0; jglobal<ny; jglobal++) {
        if (lys <= jglobal && jglobal < lye) {
          PetscInt j = jglobal-gys;
          PetscInt k = kglobal-gzs;

          PetscInt flag[3] = {0,0,0};

          if (lxe < nx) grid->receiveFlag(flag,EAST);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt i=gnx-1; i>=0; i--) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  grid->cells[ghosted_id].getZ() <= east_stage_max &&
                  grid->cells[ghosted_id].getZ() >= east_stage_min) {
                grid->cells[ghosted_id].flag |= EAST_DIR_EAST_FACE;
                flag[0] = 1;
                matrix[jglobal+kglobal*ny] = i+gxs;
                break;
              }
            }
          }
          if (lxs > 0) grid->sendFlag(flag,WEST);
        }
      }
    }
  }
  nyXnz = ny*nz;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nyXnz,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);

  for (PetscInt k=0; k<gnz; k++) {
    PetscInt kglobal = k+gzs;
    for (PetscInt j=0; j<gny; j++) {
      PetscInt jglobal = j+gys;
      // y-direction
      if (j < gny-1) {
        PetscInt iup = (PetscInt)matrix[jglobal+kglobal*ny];
        PetscInt idown = (PetscInt)matrix[jglobal+1+kglobal*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (PetscInt i=iup; i>idown; i--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= EAST_DIR_NORTH_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (PetscInt i=idown; i>iup; i--)
              grid->cells[i+(j+1)*gnx+k*gnxXny].flag |= EAST_DIR_SOUTH_FACE;     
          }
        }
      }
      // z-direction
      if (k < gnz-1) {
        PetscInt iup = (PetscInt)matrix[jglobal+kglobal*ny];
        PetscInt idown = (PetscInt)matrix[jglobal+(kglobal+1)*ny];
        if (iup>-1 && idown>-1 && iup != idown) {
          if (iup > idown) {                        
            iup = MIN(iup-gxs,gnx-1);
            idown = MAX(idown-gxs,0);
            for (PetscInt i=iup; i>idown; i--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= EAST_DIR_TOP_FACE;     
          }                                         
          else {                                    
            idown = MIN(idown-gxs,gnx-1);
            iup = MAX(iup-gxs,0);
            for (PetscInt i=idown; i>iup; i--)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= EAST_DIR_BOTTOM_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // south
  matrix = new PetscInt[nx*nz];
  for (PetscInt i=0; i<nx*nz; i++)
    matrix[i] = -1;

  for (PetscInt kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (PetscInt iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          PetscInt k = kglobal-gzs;
          PetscInt i = iglobal-gxs;

          PetscInt flag[3] = {0,0,0};

          if (lys > 0) grid->receiveFlag(flag,SOUTH);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt j=0; j<gny; j++) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  grid->cells[ghosted_id].getZ() <= south_stage_max &&
                  grid->cells[ghosted_id].getZ() >= south_stage_min) {
                grid->cells[ghosted_id].flag |= SOUTH_DIR_SOUTH_FACE;
                flag[0] = 1;
                matrix[iglobal+kglobal*nx] = j+gys;
                break;
              }
            }
          }
          if (lye < ny) grid->sendFlag(flag,NORTH);
        }
      }
    }
  }
  PetscMPIInt nxXnz = nx*nz;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nxXnz,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);

  for (PetscInt k=0; k<gnz; k++) {
    PetscInt kglobal = k+gzs;
    for (PetscInt i=0; i<gnx; i++) {
      PetscInt iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        PetscInt jup = (PetscInt)matrix[iglobal+kglobal*nx];
        PetscInt jdown = (PetscInt)matrix[iglobal+1+kglobal*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (PetscInt j=jdown; j<jup; j++)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= SOUTH_DIR_WEST_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (PetscInt j=jup; j<jdown; j++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= SOUTH_DIR_EAST_FACE;     
          }                                         
        }
      }
      // z-direction
      if (k < gnz-1) {
        PetscInt jup = (PetscInt)matrix[iglobal+kglobal*nx];
        PetscInt jdown = (PetscInt)matrix[iglobal+(kglobal+1)*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (PetscInt j=jdown; j<jup; j++)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= SOUTH_DIR_BOTTOM_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (PetscInt j=jup; j<jdown; j++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= SOUTH_DIR_TOP_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;
      
  // north
  matrix = new PetscInt[nx*nz];
  for (PetscInt i=0; i<nx*nz; i++)
    matrix[i] = -1;

  for (PetscInt kglobal=0; kglobal<nz; kglobal++) {
    if (lzs <= kglobal && kglobal < lze) {
      for (PetscInt iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          PetscInt k = kglobal-gzs;
          PetscInt i = iglobal-gxs;

          PetscInt flag[3] = {0,0,0};

          if (lye < ny) grid->receiveFlag(flag,NORTH);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt j=gny-1; j>=0; j--) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  // added to prevent connections within 100m of east boundary
//                  grid->cells[ghosted_id].getXLocal() < 400. &&
                  grid->cells[ghosted_id].getXLocal() < 750. &&
                  grid->cells[ghosted_id].getZ() <= north_stage_max &&
                  grid->cells[ghosted_id].getZ() >= north_stage_min) {
                grid->cells[ghosted_id].flag |= NORTH_DIR_NORTH_FACE;
                flag[0] = 1;
                matrix[iglobal+kglobal*nx] = j+gys;
                break;
              }
            }
          }
          if (lys > 0) grid->sendFlag(flag,SOUTH);
        }
      }
    }
  }
  nxXnz = nx*nz;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nxXnz,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);

  for (PetscInt k=0; k<gnz; k++) {
    PetscInt kglobal = k+gzs;
    for (PetscInt i=0; i<gnx; i++) {
      PetscInt iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        PetscInt jup = (PetscInt)matrix[iglobal+kglobal*nx];
        PetscInt jdown = (PetscInt)matrix[iglobal+1+kglobal*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (PetscInt j=jup; j>jdown; j--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= NORTH_DIR_EAST_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (PetscInt j=jdown; j>jup; j--)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= NORTH_DIR_WEST_FACE;     
          }                                         
        }
      }
      // z-direction
      if (k < gnz-1) {
        PetscInt jup = (PetscInt)matrix[iglobal+kglobal*nx];
        PetscInt jdown = (PetscInt)matrix[iglobal+(kglobal+1)*nx];
        if (jup>-1 && jdown>-1 && jup != jdown) {
          if (jup > jdown) {                        
            jup = MIN(jup-gys,gny-1);
            jdown = MAX(jdown-gys,0);
            for (PetscInt j=jup; j>jdown; j--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= NORTH_DIR_TOP_FACE;     
          }                                         
          else {                                    
            jdown = MIN(jdown-gys,gny-1);
            jup = MAX(jup-gys,0);
            for (PetscInt j=jdown-1; j>=jup; j--)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= NORTH_DIR_BOTTOM_FACE;
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // bottom
  matrix = new PetscInt[nx*ny];
  for (PetscInt i=0; i<nx*ny; i++)
    matrix[i] = -1;

  for (PetscInt jglobal=0; jglobal<ny; jglobal++) {
    if (lys <= jglobal && jglobal < lye) {
      for (PetscInt iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          PetscInt j = jglobal-gys;
          PetscInt i = iglobal-gxs;

          PetscInt flag[3] = {0,0,0};

          if (lzs > 0) grid->receiveFlag(flag,BOTTOM);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt k=0; k<gnz; k++) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  grid->cells[ghosted_id].getZ() <= bottom_stage_max &&
                  grid->cells[ghosted_id].getZ() >= bottom_stage_min) {
                grid->cells[ghosted_id].flag |= BOTTOM_DIR_BOTTOM_FACE;
                flag[0] = 1;
                matrix[iglobal+jglobal*nx] = k+gzs;
                break;
              }
            }
          }
          if (lze < nz) grid->sendFlag(flag,TOP);
        }
      }
    }
  }
  PetscMPIInt nxXny = nx*ny;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nxXny,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);
 
  for (PetscInt j=0; j<gny; j++) {
    PetscInt jglobal = j+gys;
    for (PetscInt i=0; i<gnx; i++) {
      PetscInt iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        PetscInt kup = (PetscInt)matrix[iglobal+jglobal*nx];
        PetscInt kdown = (PetscInt)matrix[iglobal+1+jglobal*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (PetscInt k=kdown; k<kup; k++)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= BOTTOM_DIR_WEST_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (PetscInt k=kup; k<kdown; k++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= BOTTOM_DIR_EAST_FACE;     
          }                                         
        }
      }
      // y-direction
      if (j < gny-1) {
        PetscInt kup = (PetscInt)matrix[iglobal+jglobal*nx];
        PetscInt kdown = (PetscInt)matrix[iglobal+(jglobal+1)*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (PetscInt k=kdown; k<kup; k++)
              grid->cells[i+j*gnx+(k+1)*gnxXny].flag |= BOTTOM_DIR_SOUTH_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (PetscInt k=kup; k<kdown; k++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= BOTTOM_DIR_NORTH_FACE;     
          }                                         
        }
      }
    }
  }
  delete [] matrix;
  matrix = NULL;

  // top
  matrix = new PetscInt[nx*ny];
  for (PetscInt i=0; i<nx*ny; i++)
    matrix[i] = -1;

  for (PetscInt jglobal=0; jglobal<ny; jglobal++) {
    if (lys <= jglobal && jglobal < lye) {
      for (PetscInt iglobal=0; iglobal<nx; iglobal++) {
        if (lxs <= iglobal && iglobal < lxe) {
          PetscInt j = jglobal-gys;
          PetscInt i = iglobal-gxs;

          PetscInt flag[3] = {0,0,0};

          if (lze < nz) grid->receiveFlag(flag,TOP);
          if (flag[0] == 0) {

          // loop over local cells
            for (PetscInt k=gnz-1; k>=0; k--) {
              PetscInt ghosted_id = i+j*gnx+k*gnxXny;
              if (grid->cells[ghosted_id].getIdLocal() > -1 &&
                  grid->cells[ghosted_id].getActive() && 
                  grid->cells[ghosted_id].getZ() <= top_stage_max &&
                  grid->cells[ghosted_id].getZ() >= top_stage_min) {
                grid->cells[ghosted_id].flag |= TOP_DIR_TOP_FACE;
                flag[0] = 1;
                matrix[iglobal+jglobal*nx] = k+gzs;
                break;
              }
            }
          }
          if (lzs > 0) grid->sendFlag(flag,BOTTOM);
        }
      }
    }
  }
#if 0
  nxXny = nx*ny;
  MPI_Allreduce(MPI_IN_PLACE,matrix,nxXny,MPIU_INT,MPI_MAX,
                PETSC_COMM_WORLD);
 
  for (PetscInt j=0; j<gny; j++) {
    PetscInt jglobal = j+gys;
    for (PetscInt i=0; i<gnx; i++) {
      PetscInt iglobal = i+gxs;
      // x-direction
      if (i < gnx-1) {
        PetscInt kup = (PetscInt)matrix[iglobal+jglobal*nx];
        PetscInt kdown = (PetscInt)matrix[iglobal+1+jglobal*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (PetscInt k=kup; k>kdown; k--)
              grid->cells[i+j*gnx+k*gnxXny].flag |= TOP_DIR_EAST_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (PetscInt k=kdown; k>kup; k--)
              grid->cells[i+1+j*gnx+k*gnxXny].flag |= TOP_DIR_WEST_FACE;     
          }                                         
        }
      }
      // y-direction
      if (j < gny-1) {
        PetscInt kup = (PetscInt)matrix[iglobal+jglobal*nx];
        PetscInt kdown = (PetscInt)matrix[iglobal+(jglobal+1)*nx];
        if (kup>-1 && kdown>-1 && kup != kdown) {
          if (kup > kdown) {                        
            kup = MIN(kup-gzs,gnz-1);
            kdown = MAX(kdown-gzs,0);
            for (PetscInt k=kdown; k<kup; k++)
              grid->cells[i+j*gnx+k*gnxXny].flag |= TOP_DIR_NORTH_FACE;     
          }                                         
          else {                                    
            kdown = MIN(kdown-gzs,gnz-1);
            kup = MAX(kup-gzs,0);
            for (PetscInt k=kup; k<kdown; k++)
              grid->cells[i+(j+1)*gnx+k*gnxXny].flag |= TOP_DIR_SOUTH_FACE;     
          }                                         
        }
      }
    }
  }
#endif
  delete [] matrix;
  matrix = NULL;

}

void IFC::setMaterialIdBasedOnNaturalId(PetscInt natural_id, PetscInt material_id,
                                             Grid *grid) {
  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) 
    if (grid->cells[i].getIdNatural() == natural_id)//  {
      grid->cells[i].setMaterialId(material_id);
//prPetscIntf("%d %d %d %d\n",myrank,i,grid->cells[i].getIdNatural(),grid->cells[i].getMaterialId());
//}
}

void IFC::setActiveBasedOnNaturalId(PetscInt natural_id, PetscInt active,
                                         Grid *grid) {
  for (PetscInt i=0; i<grid->getNumberOfCellsGhosted(); i++) 
    if (grid->cells[i].getIdNatural() == natural_id) 
      grid->cells[i].setActive(active);
}


IFC::~IFC() {
  if (ascii_grids) {
    for (PetscInt i=0; i<AsciiGrid::nasciigrids; i++)
      delete ascii_grids[i];
    delete [] ascii_grids;
  }
  if (river_polygon) delete river_polygon;
}

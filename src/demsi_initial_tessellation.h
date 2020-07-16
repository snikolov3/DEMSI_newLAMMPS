/*
  DEMSI initial tessellation class

 */

#ifndef DEMSI_INITIAL_TESSELLATION_H
#define DEMSI_INITIAL_TESSELLATION_H

#include "demsi_typedefs.h"
#include "demsi_polygon.h"
#include "demsi_partition.h"
#include "demsi_logging.h"
#include "demsi_configs.h"
#include "demsi_configs.h"
#include "demsi_lmp_instance.h"
#include "demsi_grid.h"
#include "demsi_particles.h"

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class Tessellation
    \brief
*/
//------------------------------------------------------------------------------
class Tessellation{
public:

  /*! \brief Constructor for Tessellation class */
  Tessellation(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Configs* configsIn, DEMSI::LammpsInstance* lammpsInstanceIn, DEMSI::Grid* gridIn, DEMSI::Particles* particlesIn);

  /*! \brief Default destructor for Tessellation class */
  ~Tessellation() = default;

  //----------------------------------------------------------------------------
  // Cells of initial tessellation
  //----------------------------------------------------------------------------

  /*! Total number of particles in initial distribution fully covering domain */
  int nParticlesInit;

  /*! Total number of particles in initial distribution in halo region */
  int nParticlesHaloInit;

  /*! Total init particles owned and halo */
  int nParticlesAllInit;

  /*! Maximum number of vertices on a cell in the initial tessellation */
  int maxVerticesInit;

  /*! Global index of initial particle */
  int *globalIndexInit;

  /*! Type of initial particle */
  int *typeInit;

  /*! x position of initial particles */
  double *xInit;

  /*! y position of initial particles */
  double *yInit;

  /*! Radius of initial particles */
  double *radiusInit;

  /*! Latitude position of initial particles */
  double *latitudeInit;

  /*! Longitude position of initial particles */
  double *longitudeInit;

  /*! Number of particles in initial tessellation on each processor */
  int* nParticlesInitOnAllProcs;

  /*! Number of particles in initial tessellation on previous (iProc < thisProc) processor */
  int* nParticlesInitOnPrevProc;

  /*! Map from local to global particle indices for initial tessellation */
  std::vector<int> localToGlobalIndicesInit;

  /*! Maxmimum size, vertex to vertex, of any cell in initial distribution */
  double maxCellSizeInit;

  /*! Maxmimum size, center to vertex, of any cell in initial distribution */
  double maxHalfCellSizeInit;

  //----------------------------------------------------------------------------
  // Vertices of cells of initial tessellation
  //----------------------------------------------------------------------------

  /*! Polygon tesselation of initial particles */
  std::vector<DEMSI::Polygon> initPolygons;

  /*! Number of vertices on a particular cell */
  int *nVerticesOnCellInit;

  /*! Maximum number of vertices on a particular cell */
  int nVerticesOnCellInitMax;

  /*! x position of vertices of initial tessellation */
  double** xVertexInit;

  /*! y position of vertices of initial tessellation */
  double** yVertexInit;

  //----------------------------------------------------------------------------
  // Edges of cells of initial tessellation
  //----------------------------------------------------------------------------

  /* Number of owned edges in initial tessellation */
  int nEdgesInit;

  /* Number of halo edges in the initial tessellation */
  int nEdgesHaloInit;

  /* Total number of edges in the initial tessellation */
  int nEdgesAllInit;

  /*! Number of edges in initial tessellation on each processor */
  int* nEdgesInitOnAllProcs;

  /*! Number of edges in initial tessellation on previous (iProc < thisProc) processor */
  int* nEdgesInitOnPrevProc;

  //----------------------------------------------------------------------------
  // Connectivity
  //----------------------------------------------------------------------------

  /*! Global indices of cells neighbouring cells */
  int** cellsOnCellGlobal;

  /*! Local indices of cells neighbouring cells */
  int** cellsOnCellLocal;

  /*! Global indices of edges neighbouring cells */
  int** edgesOnCellGlobal;

  /*! Local indices of edges neighbouring cells */
  int** edgesOnCellLocal;

  /*! Global indices of cells neighbouring edges */
  int** cellsOnEdgeGlobal;

  /*! Local indices of cells neighbouring edges */
  int** cellsOnEdgeLocal;

  //----------------------------------------------------------------------------
  // Halo
  //----------------------------------------------------------------------------

  /*! Local indices of distant halo cells mapped by proc ID - these are local owned cells */
  std::map<int, std::vector<int>> localOwnedParticlesIndicesExchInit;

  /*! local indices of distant owned cells mapped by processor ID - these are local halo cells */
  std::map<int, std::vector<int>> localHaloParticlesIndicesExchInit;

  /*! Local indices of distant halo edges mapped by proc ID - these are local owned edges */
  std::map<int, std::vector<int>> localOwnedEdgesIndicesExchInit;

  /*! local indices of distant owned edges mapped by processor ID - these are local halo edges */
  std::map<int, std::vector<int>> localHaloEdgesIndicesExchInit;

  /*! Perform halo exchange for particle integer vector */
  void update_init_halo_particles(std::vector<int> &array);

  /*! Perform halo exchange for particle double vector */
  void update_init_halo_particles(std::vector<double> &array);

  /*! Perform halo exchange for particle integer array */
  void update_init_halo_particles(int *array);

  /*! Perform halo exchange for particle double array */
  void update_init_halo_particles(double* array);

  /*! Perform halo exchange for edge integer vector */
  void update_init_halo_edges(std::vector<int> &array);

  /*! Perform halo exchange for edge double vector */
  void update_init_halo_edges(std::vector<double> &array);

  /*! Perform halo exchange for edge integer array */
  void update_init_halo_edges(int *array);

  /*! Perform halo exchange for edge double array */
  void update_init_halo_edges(double* array);

  //----------------------------------------------------------------------------
  // IO
  //----------------------------------------------------------------------------

  /*! Vector of fields registered for possible output */
  std::vector<std::pair<double*,std::string>> fieldsWrite;

  /*! Register fields for output */
  void register_output_field(double* field, const std::string fieldName);

private:

  /*! Pointer to partition object */
  DEMSI::Partition* partition;

  /*! Pointer to Log object */
  DEMSI::Log* log;

  /*! Pointer to configs object */
  DEMSI::Configs* configs;

  /*! Pointer to lammps instance */
  DEMSI::LammpsInstance* lammpsInstance;

  /*! Pointer to grid object */
  DEMSI::Grid* grid;

  /*! Pointer to particles object */
  DEMSI::Particles* particles;

  /*! Setup initial particle distribution */
  void init_initial_distribution(void);

  /*! Read in the initial particle distribution */
  void read_initial_distribution(std::vector<int> &initHaloIndices);

  /*! Initialize halo exchange for initial tessellation particles */
  void init_particles_halo(std::vector<int> initHaloIndices);

  /*! Create the polygon objects of the initial tessellation */
  void set_init_particles_polygons(void);

  /*! Determine local cellsOnCell index from global */
  void set_local_cells_on_cell(void);

  /*! Set the edge to cell connectivity */
  void set_edge_connectivity(void);

  /*! Initialize halo exchange for initial tessellation edges */
  void init_edge_halo(void);

  /*! Write out diagnostics for the initial particle distribution */
  void init_distribution_diagnostics(void);

  /*! Set mapping from local to global indices for initial tessellation particles. */
  void set_local_to_global_index_map_init(void);

  /*! Output particle array to test file for debugging */
  void output_init_particles(const std::string filenamePrefix, double* array);

  /*! Output edge array to test file for debugging */
  void output_init_edges(const std::string filenamePrefix, double* array);

};

}

#endif

#ifndef DEMSI_COLUMN_TRACERS_H_
#define DEMSI_COLUMN_TRACERS_H_

/*! \file   demsi_column_tracers.h
    \brief  Header file for the DEMSI::ColumnTracers class
*/

#include "demsi_configs.h"
#include "demsi_column_variables.h"

namespace DEMSI {

/*! \class ColumnTracers
 \brief Class describing the hierarchical relationship of the column tracers
*/
class ColumnTracers {

public:

  //-----------------------------------------------------------------------
  // tracer parameters
  //-----------------------------------------------------------------------

  /*! If true use ice age tracer */
  bool useIceAge;
  /*! If true use first year ice tracer */
  bool useFirstYearIce;
  /*! If true use level ice tracer */
  bool useLevelIce;
  /*! If true use melt pond tracer */
  bool useMeltPonds;
  /*! If true use CESM melt pond tracer */
  bool useCesmMeltponds;
  /*! If true use level ice meltpond tracer */
  bool useLevelMeltponds;
  /*! If true use the topo meltpond tracer */
  bool useTopoMeltponds;
  /*! If true use aerosol tracers */
  bool useAerosols;

  //-----------------------------------------------------------------------
  // base tracer object
  //-----------------------------------------------------------------------

  /*! length of tracer array */
  int nTracers; //ntrcr

  /*! number of base tracers */
  int nBaseTracers = 3;

  /*! maximum number of ancestor tracers */
  int nMaxAncestorTracers = 2;

  /*! category tracer array: trcrn */
  double* tracerArrayCategory;

  /*! cell tracer array: trcr */
  double* tracerArray;

  /*! index of the parent tracer: trcr_depend */
  int* parentIndex;

  /*! first ancestor type mask: trcr_base */
  double* firstAncestorMask;

  /*! indices of ancestor tracers excluding base tracer: nt_strata */
  int* ancestorIndices;

  /*! number of ancestor tracers excluding base tracer: n_trcr_strata */
  int* ancestorNumber;

  //-----------------------------------------------------------------------
  // physics
  //-----------------------------------------------------------------------

  // indexes of physics tracers in tracer array
  /*! tracer index for surface temperature: nt_Tsfc */
  int indexSurfaceTemperature;
  /*! tracer index for ice enthalpy: nt_qice */
  int indexIceEnthalpy;
  /*! tracer index for snow enthalpy: nt_qsno */
  int indexSnowEnthalpy;
  /*! tracer index for ice salinity: nt_sice */
  int indexIceSalinity;
  /*! tracer index for ice age: nt_iage */
  int indexIceAge;
  /*! tracer index for first year ice tracer: nt_FY */
  int indexFirstYearIceArea;
  /*! tracer index for level ice area: nt_alvl */
  int indexLevelIceArea;
  /*! tracer index for level ice volume: nt_vlvl */
  int indexLevelIceVolume;
  /*! tracer index for melt pond area: nt_apnd */
  int indexPondArea;
  /*! tracer index for melt pond depth: nt_hpnd */
  int indexPondDepth;
  /*! tracer index for melt pond refrozen lid thickness: nt_ipnd */
  int indexPondLidThickness;
  /*! tracer index for aerosols: nt_aero */
  int indexAerosols;

  //-----------------------------------------------------------------------
  // public member functions
  //-----------------------------------------------------------------------

  /*! \brief ColumnTracers class constructor
      \param columnDimensionsIn Pointer to the column dimensions object.
      \param configsIn Pointer to the configs object
  */
  ColumnTracers(DEMSI::ColumnDimensions* columnDimensionsIn, DEMSI::Configs* configsIn);

  /*! \brief default destructor */
  ~ColumnTracers() = default;

  /*! \brief Print the tracer object */
  void print(void);

  /*! \brief Set pointers to the column tracer variables
      \param iceAreaCategoryIn Pointer to the ice area tracer
      \param iceVolumeCategoryIn Pointer to the ice volume tracer
      \param snowVolumeCategoryIn Pointer to the snow volume tracer
      \param surfaceTemperatureIn Pointer to the surface temperature tracer
      \param iceEnthalpyIn Pointer to the ice enthalpy tracer
      \param iceSalinityIn Pointer to the ice salinity tracer
      \param snowEnthalpyIn Pointer to the snow enthalpy tracer
      \param iceAgeIn Pointer to the ice age tracer
      \param firstYearIceAreaIn Pointer to the first year ice tracer
      \param levelIceAreaIn Pointer to the level ice area tracer
      \param levelIceVolumeIn Pointer to the level ice volume tracer
      \param pondAreaIn Pointer to the meltpond area tracer
      \param pondDepthIn Pointer to the meltpond depth tracer
      \param pondLidThicknessIn Pointer to the meltpond refrozen lid thickness tracer
      \param snowScatteringAerosolIn Pointer to the snow scattering aerosol tracer
      \param snowBodyAerosolIn Pointer to the snow body aerosol tracer
      \param iceScatteringAerosolIn Pointer to the ice scattering aerosol tracer
      \param iceBodyAerosolIn Pointer to the ice body aerosol tracer
  */
  void set_tracer_variable_pointers(
	DEMSI::ColumnVariable<double>* iceAreaCategoryIn,
	DEMSI::ColumnVariable<double>* iceVolumeCategoryIn,
	DEMSI::ColumnVariable<double>* snowVolumeCategoryIn,
	DEMSI::ColumnVariable<double>* surfaceTemperatureIn,
	DEMSI::ColumnVariable<double>* iceEnthalpyIn,
	DEMSI::ColumnVariable<double>* iceSalinityIn,
	DEMSI::ColumnVariable<double>* snowEnthalpyIn,
	DEMSI::ColumnVariable<double>* iceAgeIn,
	DEMSI::ColumnVariable<double>* firstYearIceAreaIn,
	DEMSI::ColumnVariable<double>* levelIceAreaIn,
	DEMSI::ColumnVariable<double>* levelIceVolumeIn,
	DEMSI::ColumnVariable<double>* pondAreaIn,
	DEMSI::ColumnVariable<double>* pondDepthIn,
	DEMSI::ColumnVariable<double>* pondLidThicknessIn,
	DEMSI::ColumnVariable<double>* snowScatteringAerosolIn,
	DEMSI::ColumnVariable<double>* snowBodyAerosolIn,
	DEMSI::ColumnVariable<double>* iceScatteringAerosolIn,
	DEMSI::ColumnVariable<double>* iceBodyAerosolIn);

  /*! \brief Set pointers to the column aggregated tracer variables
      \param iceAreaCellIn Pointer to the ice area aggregated tracer
      \param iceVolumeCellIn Pointer to the ice volume aggregated tracer
      \param snowVolumeCellIn Pointer to the snow volume aggregated tracer
      \param surfaceTemperatureCellIn Pointer to the surface temperature aggregated tracer
      \param iceEnthalpyCellIn Pointer to the ice enthalpy aggregated tracer
      \param iceSalinityCellIn Pointer to the ice salinity aggregated tracer
      \param snowEnthalpyCellIn Pointer to the snow enthalpy aggregated tracer
      \param iceAgeCellIn Pointer to the ice age aggregated tracer
      \param firstYearIceAreaCellIn Pointer to the first year aggregated ice tracer
      \param levelIceAreaCellIn Pointer to the level ice area aggregated tracer
      \param levelIceVolumeCellIn Pointer to the level ice volume aggregated tracer
      \param pondAreaCellIn Pointer to the meltpond area aggregated tracer
      \param pondDepthCellIn Pointer to the meltpond depth aggregated tracer
      \param pondLidThicknessCellIn Pointer to the meltpond refrozen lid thickness aggregated tracer
      \param snowScatteringAerosolCellIn Pointer to the snow scattering aerosol aggregated tracer
      \param snowBodyAerosolCellIn Pointer to the snow body aerosol aggregated tracer
      \param iceScatteringAerosolCellIn Pointer to the ice scattering aeroso aggregatedl tracer
      \param iceBodyAerosolCellIn Pointer to the ice body aerosol aggregated tracer
  */
  void set_aggregated_tracer_variable_pointers(
	DEMSI::ColumnVariable<double>* iceAreaCellIn,
	DEMSI::ColumnVariable<double>* iceVolumeCellIn,
	DEMSI::ColumnVariable<double>* snowVolumeCellIn,
	DEMSI::ColumnVariable<double>* surfaceTemperatureCellIn,
	DEMSI::ColumnVariable<double>* iceEnthalpyCellIn,
	DEMSI::ColumnVariable<double>* iceSalinityCellIn,
	DEMSI::ColumnVariable<double>* snowEnthalpyCellIn,
	DEMSI::ColumnVariable<double>* iceAgeCellIn,
	DEMSI::ColumnVariable<double>* firstYearIceAreaCellIn,
	DEMSI::ColumnVariable<double>* levelIceAreaCellIn,
	DEMSI::ColumnVariable<double>* levelIceVolumeCellIn,
	DEMSI::ColumnVariable<double>* pondAreaCellIn,
	DEMSI::ColumnVariable<double>* pondDepthCellIn,
	DEMSI::ColumnVariable<double>* pondLidThicknessCellIn,
	DEMSI::ColumnVariable<double>* snowScatteringAerosolCellIn,
	DEMSI::ColumnVariable<double>* snowBodyAerosolCellIn,
	DEMSI::ColumnVariable<double>* iceScatteringAerosolCellIn,
	DEMSI::ColumnVariable<double>* iceBodyAerosolCellIn);

  /*! \brief Set the tracer category arrays for a particle
      \param iParticle Index of particle to set array for
  */
  void set_icepack_tracer_array_category(const int iParticle);

  /*! \brief Get the tracer category arrays for a particle
      \param iParticle Index of particle to get array for
  */
  void get_icepack_tracer_array_category(const int iParticle);

  /*! \brief Set the tracer particle arrays for a particle
      \param iParticle Index of particle to set array for
  */
  void set_icepack_tracer_array(const int iParticle);

  /*! \brief Get the tracer particle arrays for a particle
      \param iParticle Index of particle to get array for
  */
  void get_icepack_tracer_array(const int iParticle);

private:

  /*! Pointers to tracers */
  DEMSI::ColumnVariable<double>* iceAreaCategory;
  DEMSI::ColumnVariable<double>* iceVolumeCategory;
  DEMSI::ColumnVariable<double>* snowVolumeCategory;
  DEMSI::ColumnVariable<double>* surfaceTemperature;
  DEMSI::ColumnVariable<double>* iceEnthalpy;
  DEMSI::ColumnVariable<double>* iceSalinity;
  DEMSI::ColumnVariable<double>* snowEnthalpy;
  DEMSI::ColumnVariable<double>* iceAge;
  DEMSI::ColumnVariable<double>* firstYearIceArea;
  DEMSI::ColumnVariable<double>* levelIceArea;
  DEMSI::ColumnVariable<double>* levelIceVolume;
  DEMSI::ColumnVariable<double>* pondArea;
  DEMSI::ColumnVariable<double>* pondDepth;
  DEMSI::ColumnVariable<double>* pondLidThickness;
  DEMSI::ColumnVariable<double>* snowScatteringAerosol;
  DEMSI::ColumnVariable<double>* snowBodyAerosol;
  DEMSI::ColumnVariable<double>* iceScatteringAerosol;
  DEMSI::ColumnVariable<double>* iceBodyAerosol;

  /*! Pointers to aggregated tracers */
  DEMSI::ColumnVariable<double>* iceAreaCell;
  DEMSI::ColumnVariable<double>* iceVolumeCell;
  DEMSI::ColumnVariable<double>* snowVolumeCell;
  DEMSI::ColumnVariable<double>* surfaceTemperatureCell;
  DEMSI::ColumnVariable<double>* iceEnthalpyCell;
  DEMSI::ColumnVariable<double>* iceSalinityCell;
  DEMSI::ColumnVariable<double>* snowEnthalpyCell;
  DEMSI::ColumnVariable<double>* iceAgeCell;
  DEMSI::ColumnVariable<double>* firstYearIceAreaCell;
  DEMSI::ColumnVariable<double>* levelIceAreaCell;
  DEMSI::ColumnVariable<double>* levelIceVolumeCell;
  DEMSI::ColumnVariable<double>* pondAreaCell;
  DEMSI::ColumnVariable<double>* pondDepthCell;
  DEMSI::ColumnVariable<double>* pondLidThicknessCell;
  DEMSI::ColumnVariable<double>* snowScatteringAerosolCell;
  DEMSI::ColumnVariable<double>* snowBodyAerosolCell;
  DEMSI::ColumnVariable<double>* iceScatteringAerosolCell;
  DEMSI::ColumnVariable<double>* iceBodyAerosolCell;

  /*! Pointer to the configs object */
  DEMSI::Configs* configs;

  /*! Pointer to the column dimensions */
  DEMSI::ColumnDimensions* columnDimensions;

  // routines to init the various parts of the tracer object
  // and initialize tracers quantities in Icepack

  /*! \brief Set which tracers to use in the tracer object */
  void get_use_tracers(void);

  /*! \brief Set the number of tracers in the tracer object */
  void get_tracer_number(void);

  /*! \brief Allocate the tracer info arrays in the tracer object */
  void allocate_arrays(void);

  /*! \brief Init the child tracer indices in the tracer object */
  void init_child_indices(void);

  /*! \brief Init the parent tracer indices in the tracer object */
  void init_parent_indices(void);

  /*! \brief Set the first ancestor mask in the tracer object */
  void set_first_ancestor_mask(void);

  /*! \brief Set the ancester indices in the tracer object */
  void set_ancester_indices(void);

  /*! \brief Set the use tracer flags in icepack */
  void set_icepack_tracer_flags(void);

  /*! \brief Set the number of tracers in icepack */
  void set_icepack_tracer_numbers(void);

  /*! \brief Set the tracer indices in icepack */
  void set_icepack_tracer_indices(void);

}; // ColumnTracers class

/*! \class TracerTree
 \brief Class defining a node of the tracer hierarchy tree
*/
class TracerTree {

public:

  /*! \brief TracerTree constructor
      \param tracerIn Pointer to tracer variable to create the node from
  */
  TracerTree(DEMSI::ColumnVariable<double>* tracerIn);

  /*! \brief TracerTree constructor
      \param tracerIn Pointer to tracer variable to create the node from
      \param parentIn Pointer to the parent tracer variable of the tracer variable
  */
  TracerTree(DEMSI::ColumnVariable<double>* tracerIn, DEMSI::ColumnVariable<double>* parentIn);

  /*! \brief default destructor */
  ~TracerTree() = default;

  /*! \brief Add a tracer node to the tracer tree
      \param tracer Tracer variable to add
      \param parent Parent tracer for tracer to add
  */
  void add(DEMSI::ColumnVariable<double>* tracer, DEMSI::ColumnVariable<double>* parent);

  /*! \brief Print out the tracers in the tracer tree */
  void print(int level = 0);

  /*! \brief Get a vector of all the tracer variable pointers in the tree
      \return A vector of all the tracer variable pointers in the tree
  */
  std::vector<DEMSI::ColumnVariable<double>*> get_tracer_ptrs(void);

  /*! \brief Set the previous pointers of the tracer tree
      \param prev Pointer to the previous tracer, starts as NULL
      \return pointer to last tracer in order
  */
  DEMSI::TracerTree* set_prev_pointer(DEMSI::TracerTree* prev);

  /*! \brief Return the name of the underlying tracer
      \return Name of the underlying tracer
  */
  std::string name(void);

  /*! Pointer to the underlying tracer variable */
  DEMSI::ColumnVariable<double>* tracer;

  /*! Pointer to the parent variable of the tracer node */
  DEMSI::ColumnVariable<double>* parent = NULL;

  /*! Pointer to next tracer in order */
  DEMSI::TracerTree* next;

  /*! Pointer to prev tracer in reverse order */
  DEMSI::TracerTree* prev;

private:

  /*! Children of this node of the tracer tree */
  std::vector<DEMSI::TracerTree*> children;

}; // tracersTree class

} // namespace DEMSI

#endif /* COLUMN_TRACERS_H_ */

#include "demsi_interpolation.h"
#include "demsi_logging.h"
#include <iostream>
#include <fstream>
#include <cmath>

#include <lammps/atom.h>
#include <lammps/math_const.h>

namespace DEMSI {

  Interpolation::Interpolation(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Configs* configsIn, 
                               DEMSI::Grid* gridIn, DEMSI::Particles* particlesIn, DEMSI::Column* columnIn){


    partition = partitionIn;
    log = logIn;
    configs = configsIn;
    grid = gridIn;
    particles = particlesIn;
    column = columnIn;

    gridPartitions = grid->getGridPartitions();

    interpToGrid = 0;

    // find fields that need interpolation
    if (configs->exists({"ConfigGroup:griddedOutput","Array:griddedOutputStreams"})) {

    // get gridded output stream names
    std::vector<std::string> outputStreamNames = configs->get_array({"ConfigGroup:griddedOutput","Array:griddedOutputStreams"}, "GriddedOutputStream", "name");

       // iterate over output streams
       for (int iStream=0 ; iStream<outputStreamNames.size() ; iStream++) {

           std::list<std::string> configLocationRoot;
           configLocationRoot.push_back("ConfigGroup:griddedOutput");
           configLocationRoot.push_back("Array:griddedOutputStreams");
           configLocationRoot.push_back("GriddedOutputStream:" + outputStreamNames[iStream]);

           // get output field names
           std::list<std::string> outputFieldsLocation = configLocationRoot;
           outputFieldsLocation.push_back("Array:OutputFields");
           std::vector<std::string> outputFieldNames = configs->get_array(outputFieldsLocation, "OutputField");
 
           // loop over output fields
           for (int ifield = 0; ifield<outputFieldNames.size(); ifield++) {
              if (outputFieldNames[ifield] == "iceThicknessGrid"){
              //    iceThicknessGrid = new  Kokkos::View<double **>("iceThicknessGrid",grid->nx(),grid->ny());
                  iceThicknessGrid = Kokkos::View<double **>("iceThicknessGrid",grid->nx(),grid->ny());
                  gridFieldsInterp.push_back(iceThicknessGrid);
                  partFieldsInterp.push_back(column->iceVolumeCell);
                  grid->register_field_for_write(&iceThicknessGrid);
                  interpToGrid = 1;
              }
              if (outputFieldNames[ifield] == "iceConcentrationGrid"){
                  iceConcentrationGrid = Kokkos::View<double **>("iceConcentrationGrid",grid->nx(),grid->ny());
                  gridFieldsInterp.push_back(iceConcentrationGrid);
                  partFieldsInterp.push_back(column->iceAreaCell);
                  grid->register_field_for_write(&iceConcentrationGrid);
                  interpToGrid = 1;
              }
           }
        }
     }
  }

  /*===========================================================================================
   * Interpolate array of column variables to grid
   *===========================================================================================*/
  void Interpolation::particle_fields_to_grid() {

    // compute moment matrix
    Kokkos::View<double ***> momentMat = Kokkos::View<double ***> ("momentMatrix",grid->nx(),grid->ny(),7);
    Kokkos::deep_copy(momentMat, 0);
    Kokkos::fence();
    get_moment_matrix(momentMat);

    // communicate boundary values
    exchange_moment_matrix(momentMat);

    // compute moment matrix inverse (only need the first row)
    Kokkos::View<double ***> momentMatInv = Kokkos::View<double ***> ("momentMatrixInverse",grid->nx(),grid->ny(),3);
    Kokkos::deep_copy(momentMatInv, 0);
    Kokkos::fence();
    get_moment_matrix_inv(momentMat, momentMatInv);

    // interpolate fields
    for (int ifield = 0; ifield < partFieldsInterp.size(); ifield++){

        Kokkos::View<double **> minField = Kokkos::View<double **>("gridFieldMin",grid->nx(),grid->ny());
        Kokkos::View<double **> maxField = Kokkos::View<double **>("gridFieldMax",grid->nx(),grid->ny());
        get_field_bounds(partFieldsInterp[ifield], minField, maxField);
        exchange_grid_field_min(minField);
        exchange_grid_field_max(maxField);
  
        Kokkos::View<double **> gridFieldLO = Kokkos::View<double **>("gridFieldLO",grid->nx(),grid->ny());
        Kokkos::deep_copy(gridFieldLO, 0);
        interpolate_fields(momentMat, momentMatInv, partFieldsInterp[ifield], gridFieldsInterp[ifield],
                           gridFieldLO);

        exchange_grid_field(gridFieldsInterp[ifield]);
        exchange_grid_field(gridFieldLO);

        limit_grid_field(gridFieldsInterp[ifield],gridFieldLO,minField,maxField);
    }

  }

  /*===========================================================================================
   * Interpolate fields from particle to grid
   *===========================================================================================*/
   void Interpolation::interpolate_fields(Kokkos::View<double***> & momentMat, Kokkos::View<double***> & momentMatInv,
                                          DEMSI::ColumnVariable<double>* partField, Kokkos::View<double**>  gridField,
                                          Kokkos::View<double**>  gridFieldLO) {


    auto particles_x = particles->x;
    auto particles_type = particles->type;
    auto particles_id = particles->globalID;
    auto grid_dx = grid->dx();
    auto grid_dy = grid->dy();
    auto grid_nx = grid->nx();
    auto grid_ny = grid->ny();
    auto total_get_i0 = -gridPartitions[partition->this_proc()]->ilocal_total(0);
    auto total_get_j0 = -gridPartitions[partition->this_proc()]->jlocal_total(0);

    Kokkos::deep_copy(gridField, 0);
    Kokkos::deep_copy(gridFieldLO, 0);
    Kokkos::fence();

    // loop over particles
    Kokkos::parallel_for("Interpolation::particles_to_grid", Kokkos::RangePolicy<device_execution_space>
        (0, *(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

      if (particles_type(iParticle) == 1){

      //Get particle position
      double x = particles_x(iParticle,0);
      double y = particles_x(iParticle,1);

      //Locate particle on grid
      int iglobal,jglobal;
      double xi,eta;
      locate_particle(x,y,grid_dx,grid_dy,iglobal,jglobal,xi,eta);

      //Get basis values
      double N[4];
      grid_basis(xi,eta,N);

      //Loop over cell nodes
      int i[4];
      int j[4];
      i[0] = iglobal;
      j[0] = jglobal;
      i[1] = i[0]+1; i[2] = i[0]+1; i[3] = i[0];
      j[1] = j[0];   j[2] = j[0]+1; j[3] = j[0]+1;

      double partval = (*partField)(iParticle);

      for (int k=0; k<4; ++k){

        int iloc = i[k] - total_get_i0;
        int jloc = j[k] - total_get_j0;
        double xnode = i[k]*grid_dx;
        double ynode = j[k]*grid_dy;

        // compute interpolation coefficient 
        double psi = momentMatInv(iloc,jloc,0) + momentMatInv(iloc,jloc,1)*(xnode - x) 
                     + momentMatInv(iloc,jloc,2)*(ynode - y);  
        
        // evaluate grid fields
         //gridField(iloc,jloc) += N[k]*psi*partval;
         Kokkos::atomic_add(&gridField(iloc,jloc), N[k]*psi*partval);
         Kokkos::atomic_add(&gridFieldLO(iloc,jloc), N[k]*partval/momentMat(iloc,jloc,0));
      }
     }
    }); // end loop over particles
    Kokkos::fence();

  }

  /*===========================================================================================
   * Get local particle field bounds at grid points
   *===========================================================================================*/
   void Interpolation::get_field_bounds(DEMSI::ColumnVariable<double>* partField, Kokkos::View<double**>  & minField,
                                          Kokkos::View<double**>  & maxField) {


    auto particles_x = particles->x;
    auto particles_type = particles->type;
    auto grid_dx = grid->dx();
    auto grid_dy = grid->dy();
    auto total_get_i0 = -gridPartitions[partition->this_proc()]->ilocal_total(0);
    auto total_get_j0 = -gridPartitions[partition->this_proc()]->jlocal_total(0);

    Kokkos::deep_copy(minField, 1.0e3);
    Kokkos::deep_copy(maxField, -1.0e3);
    //Kokkos::deep_copy(minField, 0);
    //Kokkos::deep_copy(maxField, 0);
    Kokkos::fence();

    // loop over particles
    Kokkos::parallel_for("Interpolation::get_field_bounds", Kokkos::RangePolicy<device_execution_space>
        (0, *(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

      if (particles_type(iParticle) == 1){

      //Get particle position
      double x = particles_x(iParticle,0);
      double y = particles_x(iParticle,1);

      //Locate particle on grid
      int iglobal,jglobal;
      double xi,eta;
      locate_particle(x,y,grid_dx,grid_dy,iglobal,jglobal,xi,eta);

      //Loop over cell nodes
      int i[4];
      int j[4];
      i[0] = iglobal;
      j[0] = jglobal;
      i[1] = i[0]+1; i[2] = i[0]+1; i[3] = i[0];
      j[1] = j[0];   j[2] = j[0]+1; j[3] = j[0]+1;

      double partval = (*partField)(iParticle);
      for (int k=0; k<4; ++k){

        int iloc = i[k] - total_get_i0;
        int jloc = j[k] - total_get_j0;
        double xnode = i[k]*grid_dx;
        double ynode = j[k]*grid_dy;

        // get min/max
      /* 
        if (minField(iloc,jloc) == 0) {
            minField(iloc,jloc) = partval;
        }
        else {
            minField(iloc,jloc) = std::min(partval,minField(iloc,jloc));
        }
        if (maxField(iloc,jloc) == 0) {
            maxField(iloc,jloc) = partval;
        }
        else {
            maxField(iloc,jloc) = std::max(partval,maxField(iloc,jloc));
        }
       */
         minField(iloc,jloc) = std::min(partval,minField(iloc,jloc));
         maxField(iloc,jloc) = std::max(partval,maxField(iloc,jloc));
      }
     }
    }); // end loop over particles
    Kokkos::fence();

  }

  /*===========================================================================================
   * Get local particle field bounds at grid points
   *===========================================================================================*/
   void Interpolation::limit_grid_field(Kokkos::View<double**>  gridField, Kokkos::View<double**>  gridFieldLO,
                                        Kokkos::View<double**>  & minField, Kokkos::View<double**>  & maxField) {


    typedef typename Kokkos::Experimental::MDRangePolicy< Kokkos::Experimental::Rank<2> > MDPolicyType_2D;
    MDPolicyType_2D mdpolicy_2d( {{0,0}}, {{ grid->nx(), grid->ny() }} );

    Kokkos::parallel_for( "ForcingField::limit_grid_field", mdpolicy_2d,
    KOKKOS_LAMBDA( const int i, const int j ) {
      
        if (gridField(i,j) > maxField(i,j)) 
            gridField(i,j) = gridFieldLO(i,j);
        if (gridField(i,j) < minField(i,j)) 
            gridField(i,j) = gridFieldLO(i,j);
      
    });
    Kokkos::fence();

   }

  /*===========================================================================================
   * Interpolate scalar particle field to Kokkos view grid field
   *===========================================================================================*/
  void Interpolation::particles_to_grid(Kokkos::View<double *> &partField, Kokkos::View<double **> &gridField) {

    // compute moment matrix
    Kokkos::View<double ***> momentMat = Kokkos::View<double ***> ("momentMatrix",grid->nx(),grid->ny(),7);
    Kokkos::deep_copy(momentMat, 0);
    Kokkos::fence();
    get_moment_matrix(momentMat);

    // communicate boundary values
    exchange_moment_matrix(momentMat);

    // compute moment matrix inverse (only need the first row)
    Kokkos::View<double ***> momentMatInv = Kokkos::View<double ***> ("momentMatrixInverse",grid->nx(),grid->ny(),3);
    Kokkos::fence();
    Kokkos::deep_copy(momentMatInv, 0);
    get_moment_matrix_inv(momentMat, momentMatInv);

    // ensure grid values are zero to start
    Kokkos::deep_copy(gridField, 0);
    Kokkos::fence();

    auto particles_x = particles->x;
    auto particles_type = particles->type;
    auto grid_dx = grid->dx();
    auto grid_dy = grid->dy();
    auto grid_nx = grid->nx();
    auto grid_ny = grid->ny();
    auto total_get_i0 = -gridPartitions[partition->this_proc()]->ilocal_total(0);
    auto total_get_j0 = -gridPartitions[partition->this_proc()]->jlocal_total(0);
    // loop over particles
    Kokkos::parallel_for("Interpolation::particles_to_grid", Kokkos::RangePolicy<device_execution_space>
        (0, *(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

      if (particles_type(iParticle) == 1) {

      //Get particle position
      double x = particles_x(iParticle,0);
      double y = particles_x(iParticle,1);

      //Locate particle on grid
      int iglobal,jglobal;
      double xi,eta;
      locate_particle(x,y,grid_dx,grid_dy,iglobal,jglobal,xi,eta);

      //Get basis values
      double N[4];
      grid_basis(xi,eta,N);

      //Loop over cell nodes
      int i[4];
      int j[4];
      i[0] = iglobal;
      j[0] = jglobal;
      i[1] = i[0]+1; i[2] = i[0]+1; i[3] = i[0];
      j[1] = j[0];   j[2] = j[0]+1; j[3] = j[0]+1;
      for (int k=0; k<4; ++k){

        int iloc = i[k] - total_get_i0;
        int jloc = j[k] - total_get_j0;
        double xnode = i[k]*grid_dx;
        double ynode = j[k]*grid_dy;

        // compute interpolation coefficient 
        double psi = momentMatInv(iloc,jloc,0) + momentMatInv(iloc,jloc,1)*(xnode - x) 
                     + momentMatInv(iloc,jloc,2)*(ynode - y);  
        
        // evaluate grid field
        Kokkos::atomic_add(&gridField(iloc,jloc), N[k]*psi*partField(iParticle));

      }
    }
    }); // end loop over particles
    Kokkos::fence();

    exchange_grid_field(gridField);
  }

  /*===========================================================================================
   * Interpolate Kokkos view vector particle field to Kokkos view grid field
   *===========================================================================================*/
  void Interpolation::particles_to_grid(Kokkos::View<double **> &partField, Kokkos::View<double ***> &gridField) {

    // compute moment matrix
    Kokkos::View<double ***> momentMat = Kokkos::View<double ***> ("momentMatrix",grid->nx(),grid->ny(),7);
    Kokkos::deep_copy(momentMat, 0);
    Kokkos::fence();
    get_moment_matrix(momentMat);

    // communicate boundary values
    exchange_moment_matrix(momentMat);

    // compute moment matrix inverse (only need the first row)
    Kokkos::View<double ***> momentMatInv = Kokkos::View<double ***> ("momentMatrixInverse",grid->nx(),grid->ny(),3);
    Kokkos::deep_copy(momentMatInv, 0);
    Kokkos::fence();
    get_moment_matrix_inv(momentMat, momentMatInv);

    // ensure grid values are zero to start
    Kokkos::deep_copy(gridField, 0);
    Kokkos::fence();

    auto particles_x = particles->x;
    auto particles_type = particles->type;
    auto grid_dx = grid->dx();
    auto grid_dy = grid->dy();
    auto grid_nx = grid->nx();
    auto grid_ny = grid->ny();
    auto total_get_i0 = -gridPartitions[partition->this_proc()]->ilocal_total(0);
    auto total_get_j0 = -gridPartitions[partition->this_proc()]->jlocal_total(0);
    // loop over particles
    Kokkos::parallel_for("Interpolation::particles_to_grid", Kokkos::RangePolicy<device_execution_space>
        (0, *(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

      if (particles_type(iParticle) == 1) {

      //Get particle position
      double x = particles_x(iParticle,0);
      double y = particles_x(iParticle,1);

      //Locate particle on grid
      int iglobal,jglobal;
      double xi,eta;
      Interpolation::locate_particle(x,y,grid_dx,grid_dy,iglobal,jglobal,xi,eta);

      //Get basis values
      double N[4];
      grid_basis(xi,eta,N);

      //Loop over cell nodes
      int i[4];
      int j[4];
      i[0] = iglobal;
      j[0] = jglobal;
      i[1] = i[0]+1; i[2] = i[0]+1; i[3] = i[0];
      j[1] = j[0];   j[2] = j[0]+1; j[3] = j[0]+1;
      for (int k=0; k<4; ++k){

        int iloc = i[k] - total_get_i0;
        int jloc = j[k] - total_get_j0;
        double xnode = i[k]*grid_dx;
        double ynode = j[k]*grid_dy;

        // compute interpolation coefficient 
        double psi = momentMatInv(iloc,jloc,0) + momentMatInv(iloc,jloc,1)*(xnode - x) 
                     + momentMatInv(iloc,jloc,2)*(ynode - y);  
        
        // evaluate grid field - for now component by component
        Kokkos::atomic_add(&gridField(iloc,jloc,0), N[k]*psi*partField(iParticle,0));
        Kokkos::atomic_add(&gridField(iloc,jloc,1), N[k]*psi*partField(iParticle,1));
       }
      }
    }); // end loop over particles
    Kokkos::fence();

  }

  /*===========================================================================================
   * Interpolate column variable to Kokkos view grid field.
   *===========================================================================================*/
  void Interpolation::particles_to_grid(DEMSI::ColumnVariable<double>* partField, Kokkos::View<double **> &gridField) {

    // compute moment matrix
    Kokkos::View<double ***> momentMat = Kokkos::View<double ***> ("momentMatrix",grid->nx(),grid->ny(),7);
    Kokkos::deep_copy(momentMat, 0);
    Kokkos::fence();
    get_moment_matrix(momentMat);

    // communicate boundary values
    exchange_moment_matrix(momentMat);

    // compute moment matrix inverse (only need the first row)
    Kokkos::View<double ***> momentMatInv = Kokkos::View<double ***> ("momentMatrixInverse",grid->nx(),grid->ny(),3);
    Kokkos::deep_copy(momentMatInv, 0);
    Kokkos::fence();
    get_moment_matrix_inv(momentMat, momentMatInv);

    // ensure grid values are zero to start
    Kokkos::deep_copy(gridField, 0);
    Kokkos::fence();

    auto particles_x = particles->x;
    auto particles_type = particles->type;
    auto grid_dx = grid->dx();
    auto grid_dy = grid->dy();
    auto grid_nx = grid->nx();
    auto grid_ny = grid->ny();
    auto total_get_i0 = -gridPartitions[partition->this_proc()]->ilocal_total(0);
    auto total_get_j0 = -gridPartitions[partition->this_proc()]->jlocal_total(0);

    // loop over particles
    Kokkos::parallel_for("Interpolation::particles_to_grid", Kokkos::RangePolicy<device_execution_space>
        (0, *(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

      if (particles_type(iParticle) == 1) {

      //Get particle position
      double x = particles_x(iParticle,0);
      double y = particles_x(iParticle,1);

      //Locate particle on grid
      int iglobal,jglobal;
      double xi,eta;
      locate_particle(x,y,grid_dx,grid_dy,iglobal,jglobal,xi,eta);

      //Get basis values
      double N[4];
      grid_basis(xi,eta,N);

      //Loop over cell nodes
      int i[4];
      int j[4];
      i[0] = iglobal;
      j[0] = jglobal;
      i[1] = i[0]+1; i[2] = i[0]+1; i[3] = i[0];
      j[1] = j[0];   j[2] = j[0]+1; j[3] = j[0]+1;
      for (int k=0; k<4; ++k){

        int iloc = i[k] - total_get_i0;
        int jloc = j[k] - total_get_j0;
        double xnode = i[k]*grid_dx;
        double ynode = j[k]*grid_dy;

        // compute interpolation coefficient 
        double psi = momentMatInv(iloc,jloc,0) + momentMatInv(iloc,jloc,1)*(xnode - x) 
                     + momentMatInv(iloc,jloc,2)*(ynode - y);  
        
        // evaluate grid field 
        double partval = (*partField)(iParticle);
         //gridField(iloc,jloc) += N[k]*psi*partval;
        Kokkos::atomic_add(&gridField(iloc,jloc), N[k]*psi*partval);
      }
    }
    }); // end loop over particles
    Kokkos::fence();
    exchange_grid_field(gridField);
  }

  /*===========================================================================================
   * Evaluate moment matrix for use in interpolating particle values to grid.
   *===========================================================================================*/
  void Interpolation::get_moment_matrix(Kokkos::View<double***> & momentMat) {

    auto particles_x = particles->x;
    auto particles_type = particles->type;
    auto grid_dx = grid->dx();
    auto grid_dy = grid->dy();

    auto total_get_i0 = -gridPartitions[partition->this_proc()]->ilocal_total(0);
    auto total_get_j0 = -gridPartitions[partition->this_proc()]->jlocal_total(0);

    // loop over particles
    Kokkos::parallel_for("Interpolation::get_moment_matrix", Kokkos::RangePolicy<device_execution_space>
        (0, *(particles->nParticles)), KOKKOS_LAMBDA (const int iParticle) {

     if (particles_type(iParticle) == 1){
      //Get particle position
      double x = particles_x(iParticle,0);
      double y = particles_x(iParticle,1);

      //Locate particle on grid
      double xi,eta;
      int iglobal,jglobal;
      locate_particle(x,y,grid_dx,grid_dy,iglobal,jglobal,xi,eta);

      //Get basis values
      double N[4];
      grid_basis(xi,eta,N);

      //Loop over cell nodes
      int i[4];
      int j[4];
      i[0] = iglobal;
      j[0] = jglobal;
      i[1] = i[0]+1; i[2] = i[0]+1; i[3] = i[0];
      j[1] = j[0];   j[2] = j[0]+1; j[3] = j[0]+1;

      for (int k=0; k<4; ++k){
         // need to get what this ilocal_total and jlocal_total are doing and 
         // then just do the calculation on the vector
         int iloc = i[k] - total_get_i0;
         int jloc = j[k] - total_get_j0;
         double xnode = i[k]*grid_dx;
         double ynode = j[k]*grid_dy;
         Kokkos::atomic_add(&momentMat(iloc,jloc,0), N[k]);
         Kokkos::atomic_add(&momentMat(iloc,jloc,1), N[k]*(xnode - x));
         Kokkos::atomic_add(&momentMat(iloc,jloc,2), N[k]*(ynode - y));
         Kokkos::atomic_add(&momentMat(iloc,jloc,3), N[k]*(xnode - x)*(xnode - x));
         Kokkos::atomic_add(&momentMat(iloc,jloc,4), N[k]*(xnode - x)*(ynode - y));
         Kokkos::atomic_add(&momentMat(iloc,jloc,5), N[k]*(ynode - y)*(ynode - y));
         Kokkos::atomic_add(&momentMat(iloc,jloc,6), 1.0);
      }
     } // if particle type = 1
    }); // end loop over particles
    Kokkos::fence();

  }

  /*===========================================================================================
   * Evaluate inverse of moment matrix 
   *===========================================================================================*/
  void Interpolation::get_moment_matrix_inv(Kokkos::View<double***> & momentMat,
                                            Kokkos::View<double***> & momentMatInv) { 

    // loop over local nodes
    auto grid_nx = grid->nx();
    auto grid_ny = grid->ny();
    Kokkos::parallel_for("Interpolation::get_moment_matrix_inv", Kokkos::RangePolicy<device_execution_space>
        (0, grid_nx), KOKKOS_LAMBDA (const int iloc) {

      for (int jloc = 0; jloc < grid_ny; jloc++ ) {

        // if enough particles contributed to moment matrix compute determinant
        if (momentMat(iloc,jloc,6) >= 3) {

          // compute determinant of moment matrix
          double Mdet = momentMat(iloc,jloc,0)*(momentMat(iloc,jloc,3)*momentMat(iloc,jloc,5)
                         - momentMat(iloc,jloc,4)*momentMat(iloc,jloc,4)) + 
                        momentMat(iloc,jloc,1)*(momentMat(iloc,jloc,4)*momentMat(iloc,jloc,2)
                         - momentMat(iloc,jloc,1)*momentMat(iloc,jloc,5)) + 
                        momentMat(iloc,jloc,2)*(momentMat(iloc,jloc,1)*momentMat(iloc,jloc,4)
                         - momentMat(iloc,jloc,2)*momentMat(iloc,jloc,3)); 
          Mdet = std::abs(Mdet);

          // check for determinant > 0, could be 0 if particles are colinear
          if (Mdet > 0.0) { 

            momentMatInv(iloc,jloc,0) = (momentMat(iloc,jloc,3)*momentMat(iloc,jloc,5)-
                                         momentMat(iloc,jloc,4)*momentMat(iloc,jloc,4))/Mdet;
            momentMatInv(iloc,jloc,1) = (momentMat(iloc,jloc,2)*momentMat(iloc,jloc,4)-
                                         momentMat(iloc,jloc,1)*momentMat(iloc,jloc,5))/Mdet;
            momentMatInv(iloc,jloc,2) = (momentMat(iloc,jloc,1)*momentMat(iloc,jloc,4)-
                                         momentMat(iloc,jloc,2)*momentMat(iloc,jloc,3))/Mdet;

          }
          else { // default to first order
            momentMatInv(iloc,jloc,0) = 1.0/(momentMat(iloc,jloc,0));
            momentMatInv(iloc,jloc,1) = 0.0;
            momentMatInv(iloc,jloc,2) = 0.0;
          }

        }
        else { // default to first order
          momentMatInv(iloc,jloc,0) = 1.0/(momentMat(iloc,jloc,0));
          momentMatInv(iloc,jloc,1) = 0.0;
          momentMatInv(iloc,jloc,2) = 0.0;
        }
      }
    });
    Kokkos::fence();
  }

  /*===========================================================================================
   * Communicate grid field boundary values and min
   *===========================================================================================*/
  void Interpolation::exchange_grid_field_min(Kokkos::View<double**> & field){

    auto grid_nx = grid->nx();
    auto grid_ny = grid->ny();

    auto iglobal0 = gridPartitions[partition->this_proc()]->iglobal_total(0);
    auto iglobal1 = gridPartitions[partition->this_proc()]->iglobal_total(grid_nx-1);
    auto jglobal0 = gridPartitions[partition->this_proc()]->jglobal_total(0);
    auto jglobal1 = gridPartitions[partition->this_proc()]->jglobal_total(grid_ny-1);

    int nx_owned = gridPartitions[partition->this_proc()]->size_owned(0);
    int ny_owned = gridPartitions[partition->this_proc()]->size_owned(1);

    auto igowned0 = gridPartitions[partition->this_proc()]->iglobal_owned(0);
    auto igowned1 = gridPartitions[partition->this_proc()]->iglobal_owned(nx_owned-1);
    auto jgowned0 = gridPartitions[partition->this_proc()]->jglobal_owned(0);
    auto jgowned1 = gridPartitions[partition->this_proc()]->jglobal_owned(ny_owned-1);

    auto nx0 = igowned0-iglobal0;
    auto nx1 = iglobal1-igowned1;
    auto ny0 = jgowned0-jglobal0;
    auto ny1 = jglobal1-jgowned1;

    auto eneigh = grid->neighbor_east();
    auto wneigh = grid->neighbor_west();
    auto sneigh = grid->neighbor_south();
    auto nneigh = grid->neighbor_north();

    char mpiErrBuffer[MPI_MAX_ERROR_STRING];
    int mpiErrLen;
    int err;
    MPI_Request request;

    double *recvbuf; 
    double *sendbuf; 

   // Receive West
    if (nx0 > 0) {
       recvbuf = new double[grid_ny*2];
       err = MPI_Irecv(&recvbuf[0], grid_ny*2, MPI_DOUBLE, wneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }

   // Send east
    if (nx1 > 0) {
       sendbuf = new double[grid_ny*2];
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          sendbuf[ind] = field(grid_nx-nx1,i);
          ind++;
          sendbuf[ind] = field(grid_nx-nx1-1,i);
          ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_ny*2, MPI_DOUBLE, eneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

    if (nx0 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          field(1,i) = std::min(field(1,i),recvbuf[ind]);
          ind++;
          field(0,i) = std::min(field(0,i),recvbuf[ind]);
          ind++;
       }
    } 

    // Receive south
    if (ny0 > 0) {
       recvbuf = new double[grid_nx*2];

       err = MPI_Irecv(&recvbuf[0], grid_nx*2, MPI_DOUBLE, sneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

    // Send north
    if (ny1 > 0) {
       sendbuf = new double[grid_nx*2];
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          sendbuf[ind] = field(i,grid_ny-ny1);
          ind++;
          sendbuf[ind] = field(i,grid_ny-ny1-1);
          ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_nx*2, MPI_DOUBLE, nneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }

     if (ny0 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          field(i,1) = std::min(field(i,1),recvbuf[ind]);
          ind++;
          field(i,0) = std::min(field(i,0),recvbuf[ind]);
          ind++;
       }
    }

    MPI_Barrier(partition->comm());

   // Receive east
    if (nx1 > 0) {
       recvbuf = new double[grid_ny*2];
       err = MPI_Irecv(&recvbuf[0], grid_ny*2, MPI_DOUBLE, eneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

   // Send west
    if (nx0 > 0) {
       sendbuf = new double[grid_ny*2];
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          sendbuf[ind] = field(1,i);
          ind++;
          sendbuf[ind] = field(0,i);
          ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_ny*2, MPI_DOUBLE, wneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }
   
    if (nx1 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          field(grid_nx-nx1,i) = std::min(field(grid_nx-nx1,i),recvbuf[ind]);
          ind++;
          field(grid_nx-nx1-1,i) = std::min(field(grid_nx-nx1-1,i),recvbuf[ind]);
          ind++;
       }
     } 


    // Receive north
    if (ny1>0) {
       recvbuf = new double[grid_nx*2];
       err = MPI_Irecv(&recvbuf[0], grid_nx*2, MPI_DOUBLE, nneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

     // Send south
     if (ny0 > 0) {
       sendbuf = new double[grid_nx*2];
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          sendbuf[ind] = field(i,1);
          ind++;
          sendbuf[ind] = field(i,0);
          ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_nx*2, MPI_DOUBLE, sneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
      }

    if (ny1 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          field(i,grid_ny-ny1) = std::min(field(i,grid_ny-ny1),recvbuf[ind]);
          ind++;
          field(i,grid_ny-ny1-1) = std::min(field(i,grid_ny-ny1-1),recvbuf[ind]);
          ind++;
       }
    }
  }

  /*===========================================================================================
   * Communicate grid field boundary values and max
   *===========================================================================================*/
  void Interpolation::exchange_grid_field_max(Kokkos::View<double**> & field){

    auto grid_nx = grid->nx();
    auto grid_ny = grid->ny();

    auto iglobal0 = gridPartitions[partition->this_proc()]->iglobal_total(0);
    auto iglobal1 = gridPartitions[partition->this_proc()]->iglobal_total(grid_nx-1);
    auto jglobal0 = gridPartitions[partition->this_proc()]->jglobal_total(0);
    auto jglobal1 = gridPartitions[partition->this_proc()]->jglobal_total(grid_ny-1);

    int nx_owned = gridPartitions[partition->this_proc()]->size_owned(0);
    int ny_owned = gridPartitions[partition->this_proc()]->size_owned(1);

    auto igowned0 = gridPartitions[partition->this_proc()]->iglobal_owned(0);
    auto igowned1 = gridPartitions[partition->this_proc()]->iglobal_owned(nx_owned-1);
    auto jgowned0 = gridPartitions[partition->this_proc()]->jglobal_owned(0);
    auto jgowned1 = gridPartitions[partition->this_proc()]->jglobal_owned(ny_owned-1);

    auto nx0 = igowned0-iglobal0;
    auto nx1 = iglobal1-igowned1;
    auto ny0 = jgowned0-jglobal0;
    auto ny1 = jglobal1-jgowned1;

    auto eneigh = grid->neighbor_east();
    auto wneigh = grid->neighbor_west();
    auto sneigh = grid->neighbor_south();
    auto nneigh = grid->neighbor_north();

    char mpiErrBuffer[MPI_MAX_ERROR_STRING];
    int mpiErrLen;
    int err;
    MPI_Request request;

    double *recvbuf; 
    double *sendbuf; 

   // Receive West
    if (nx0 > 0) {
       recvbuf = new double[grid_ny*2];
       err = MPI_Irecv(&recvbuf[0], grid_ny*2, MPI_DOUBLE, wneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }

   // Send east
    if (nx1 > 0) {
       sendbuf = new double[grid_ny*2];
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          sendbuf[ind] = field(grid_nx-nx1,i);
          ind++;
          sendbuf[ind] = field(grid_nx-nx1-1,i);
          ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_ny*2, MPI_DOUBLE, eneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

    if (nx0 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          field(1,i) = std::max(field(1,i),recvbuf[ind]);
          ind++;
          field(0,i) = std::max(field(0,i),recvbuf[ind]);
          ind++;
       }
    } 

    // Receive south
    if (ny0 > 0) {
       recvbuf = new double[grid_nx*2];

       err = MPI_Irecv(&recvbuf[0], grid_nx*2, MPI_DOUBLE, sneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

    // Send north
    if (ny1 > 0) {
       sendbuf = new double[grid_nx*2];
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
            sendbuf[ind] = field(i,grid_ny-ny1);
            ind++;
            sendbuf[ind] = field(i,grid_ny-ny1-1);
            ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_nx*2, MPI_DOUBLE, nneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }

     if (ny0 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);

       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          field(i,1) = std::max(field(i,1),recvbuf[ind]);
          ind++;
          field(i,0) = std::max(field(i,0),recvbuf[ind]);
          ind++;
       }
    }

    MPI_Barrier(partition->comm());

   // Receive east
    if (nx1 > 0) {
       recvbuf = new double[grid_ny*2];
       err = MPI_Irecv(&recvbuf[0], grid_ny*2, MPI_DOUBLE, eneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

   // Send west
    if (nx0 > 0) {
       sendbuf = new double[grid_ny*2];
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          sendbuf[ind] = field(1,i);
          ind++;
          sendbuf[ind] = field(0,i);
          ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_ny*2, MPI_DOUBLE, wneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }
   
    if (nx1 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          field(grid_nx-nx1,i) = std::max(field(grid_nx-nx1,i),recvbuf[ind]);
          ind++;
          field(grid_nx-nx1-1,i) = std::max(field(grid_nx-nx1-1,i),recvbuf[ind]);
          ind++;
       }
     } 

    // Receive north
    if (ny1>0) {
       recvbuf = new double[grid_nx*2];
       err = MPI_Irecv(&recvbuf[0], grid_nx*2, MPI_DOUBLE, nneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

    // Send south
     if (ny0 > 0) {
       sendbuf = new double[grid_nx*2];
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          sendbuf[ind] = field(i,1);
          ind++;
          sendbuf[ind] = field(i,0);
          ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_nx*2, MPI_DOUBLE, sneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
      }

    if (ny1 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          field(i,grid_ny-ny1) = std::max(field(i,grid_ny-ny1),recvbuf[ind]);
          ind++;
          field(i,grid_ny-ny1-1) = std::max(field(i,grid_ny-ny1-1),recvbuf[ind]);
          ind++;
       }
    }
  }

  /*===========================================================================================
   * Communicate grid field boundary values and sum
   *===========================================================================================*/
  void Interpolation::exchange_grid_field(Kokkos::View<double**> & field){

    auto grid_nx = grid->nx();
    auto grid_ny = grid->ny();

    auto iglobal0 = gridPartitions[partition->this_proc()]->iglobal_total(0);
    auto iglobal1 = gridPartitions[partition->this_proc()]->iglobal_total(grid_nx-1);
    auto jglobal0 = gridPartitions[partition->this_proc()]->jglobal_total(0);
    auto jglobal1 = gridPartitions[partition->this_proc()]->jglobal_total(grid_ny-1);

    int nx_owned = gridPartitions[partition->this_proc()]->size_owned(0);
    int ny_owned = gridPartitions[partition->this_proc()]->size_owned(1);

    auto igowned0 = gridPartitions[partition->this_proc()]->iglobal_owned(0);
    auto igowned1 = gridPartitions[partition->this_proc()]->iglobal_owned(nx_owned-1);
    auto jgowned0 = gridPartitions[partition->this_proc()]->jglobal_owned(0);
    auto jgowned1 = gridPartitions[partition->this_proc()]->jglobal_owned(ny_owned-1);

    auto nx0 = igowned0-iglobal0;
    auto nx1 = iglobal1-igowned1;
    auto ny0 = jgowned0-jglobal0;
    auto ny1 = jglobal1-jgowned1;

    auto eneigh = grid->neighbor_east();
    auto wneigh = grid->neighbor_west();
    auto sneigh = grid->neighbor_south();
    auto nneigh = grid->neighbor_north();

    char mpiErrBuffer[MPI_MAX_ERROR_STRING];
    int mpiErrLen;
    int err;
    MPI_Request request;

    double *recvbuf; 
    double *sendbuf; 

   // Receive West
    if (nx0 > 0) {
       recvbuf = new double[grid_ny*2];
       err = MPI_Irecv(&recvbuf[0], grid_ny*2, MPI_DOUBLE, wneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }

   // Send east
    if (nx1 > 0) {
       sendbuf = new double[grid_ny*2];
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          sendbuf[ind] = field(grid_nx-nx1,i);
          ind++;
          sendbuf[ind] = field(grid_nx-nx1-1,i);
          ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_ny*2, MPI_DOUBLE, eneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

    if (nx0 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          field(1,i) = field(1,i) + recvbuf[ind];
          ind++;
          field(0,i) = field(0,i) + recvbuf[ind];
          ind++;
       }
    } 

    // Receive south
    if (ny0 > 0) {
       recvbuf = new double[grid_nx*2];
       err = MPI_Irecv(&recvbuf[0], grid_nx*2, MPI_DOUBLE, sneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

    // Send north
    if (ny1 > 0) {
       sendbuf = new double[grid_nx*2];
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          sendbuf[ind] = field(i,grid_ny-ny1);
          ind++;
          sendbuf[ind] = field(i,grid_ny-ny1-1);
          ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_nx*2, MPI_DOUBLE, nneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }

     if (ny0 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          field(i,1) = field(i,1) + recvbuf[ind];
          ind++;
          field(i,0) = field(i,0) + recvbuf[ind];
          ind++;
       }
    }

    MPI_Barrier(partition->comm());

   // Receive east
    if (nx1 > 0) {
       recvbuf = new double[grid_ny*2];
       err = MPI_Irecv(&recvbuf[0], grid_ny*2, MPI_DOUBLE, eneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

   // Send west
    if (nx0 > 0) {
       sendbuf = new double[grid_ny*2];
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          sendbuf[ind] = field(1,i);
          ind++;
          sendbuf[ind] = field(0,i);
          ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_ny*2, MPI_DOUBLE, wneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }
   
    if (nx1 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);

       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          field(grid_nx-nx1,i) = recvbuf[ind];
          ind++;
          field(grid_nx-nx1-1,i) = recvbuf[ind];
          ind++;
       }
     } 


    // Receive north
    if (ny1>0) {
       recvbuf = new double[grid_nx*2];
       err = MPI_Irecv(&recvbuf[0], grid_nx*2, MPI_DOUBLE, nneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

     // Send south
     if (ny0 > 0) {
       sendbuf = new double[grid_nx*2];
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          sendbuf[ind] = field(i,1);
          ind++;
          sendbuf[ind] = field(i,0);
          ind++;
       }
       err = MPI_Send(&sendbuf[0], grid_nx*2, MPI_DOUBLE, sneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
      }

    if (ny1 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);

       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          field(i,grid_ny-ny1) = recvbuf[ind];
          ind++;
          field(i,grid_ny-ny1-1) = recvbuf[ind];
          ind++;
       }
    }
  }

  /*===========================================================================================
   * Communicate moment matrix boundary values and sum
   *===========================================================================================*/
  void Interpolation::exchange_moment_matrix(Kokkos::View<double***> & momentMat){

    auto grid_nx = grid->nx();
    auto grid_ny = grid->ny();

    auto iglobal0 = gridPartitions[partition->this_proc()]->iglobal_total(0);
    auto iglobal1 = gridPartitions[partition->this_proc()]->iglobal_total(grid_nx-1);
    auto jglobal0 = gridPartitions[partition->this_proc()]->jglobal_total(0);
    auto jglobal1 = gridPartitions[partition->this_proc()]->jglobal_total(grid_ny-1);

    int nx_owned = gridPartitions[partition->this_proc()]->size_owned(0);
    int ny_owned = gridPartitions[partition->this_proc()]->size_owned(1);

    auto igowned0 = gridPartitions[partition->this_proc()]->iglobal_owned(0);
    auto igowned1 = gridPartitions[partition->this_proc()]->iglobal_owned(nx_owned-1);
    auto jgowned0 = gridPartitions[partition->this_proc()]->jglobal_owned(0);
    auto jgowned1 = gridPartitions[partition->this_proc()]->jglobal_owned(ny_owned-1);

    auto nx0 = igowned0-iglobal0;
    auto nx1 = iglobal1-igowned1;
    auto ny0 = jgowned0-jglobal0;
    auto ny1 = jglobal1-jgowned1;

    auto eneigh = grid->neighbor_east();
    auto wneigh = grid->neighbor_west();
    auto sneigh = grid->neighbor_south();
    auto nneigh = grid->neighbor_north();

    char mpiErrBuffer[MPI_MAX_ERROR_STRING];
    int mpiErrLen;
    int err;

    MPI_Request request;

    double *recvbuf; 
    double *sendbuf; 

   // Receive West
    if (nx0 > 0) {
       recvbuf = new double[grid_ny*7*2];
       err = MPI_Irecv(&recvbuf[0], grid_ny*7*2, MPI_DOUBLE, wneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }

   // Send east
    if (nx1 > 0) {
       sendbuf = new double[grid_ny*7*2];
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          for (int j=0; j<7; j++) {
            sendbuf[ind] = momentMat(grid_nx-nx1,i,j);
            ind++;
            sendbuf[ind] = momentMat(grid_nx-nx1-1,i,j);
            ind++;
          }
       }
       err = MPI_Send(&sendbuf[0], grid_ny*7*2, MPI_DOUBLE, eneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

    if (nx0 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          for (int j=0; j<7; j++) {
            momentMat(1,i,j) = momentMat(1,i,j) + recvbuf[ind];
            ind++;
            momentMat(0,i,j) = momentMat(0,i,j) + recvbuf[ind];
            ind++;
          }
       }
    } 

    // Receive south
    if (ny0 > 0) {
       recvbuf = new double[grid_nx*7*2];
       err = MPI_Irecv(&recvbuf[0], grid_nx*7*2, MPI_DOUBLE, sneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

    // Send north
    if (ny1 > 0) {
       sendbuf = new double[grid_nx*7*2];
       int ind = 0;
       int jtemp = jglobal1-jgowned1;
       for (int i=0; i < grid_nx; i++){
          for (int j=0; j<7; j++) {
            sendbuf[ind] = momentMat(i,grid_ny-ny1,j);
            ind++;
            sendbuf[ind] = momentMat(i,grid_ny-ny1-1,j);
            ind++;
          }
       }
       err = MPI_Send(&sendbuf[0], grid_nx*7*2, MPI_DOUBLE, nneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }

     if (ny0 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          for (int j=0; j<7; j++) {
            momentMat(i,1,j) = momentMat(i,1,j) + recvbuf[ind];
            ind++;
            momentMat(i,0,j) = momentMat(i,0,j) + recvbuf[ind];
            ind++;
          }
       }
    }

    MPI_Barrier(partition->comm());

   // Receive east
    if (nx1 > 0) {
       recvbuf = new double[grid_ny*7*2];
       err = MPI_Irecv(&recvbuf[0], grid_ny*7*2, MPI_DOUBLE, eneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

   // Send west
    if (nx0 > 0) {
       sendbuf = new double[grid_ny*7*2];
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          for (int j=0; j<7; j++) {
            sendbuf[ind] = momentMat(1,i,j);
            ind++;
            sendbuf[ind] = momentMat(0,i,j);
            ind++;
          }
       }
       err = MPI_Send(&sendbuf[0], grid_ny*7*2, MPI_DOUBLE, wneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
     }
   
    if (nx1 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_ny; i++){
          for (int j=0; j<7; j++) {
            momentMat(grid_nx-nx1,i,j) = recvbuf[ind];
            ind++;
            momentMat(grid_nx-nx1-1,i,j) = recvbuf[ind];
            ind++;
          }
       }
     } 

    // Receive north
    if (ny1>0) {
       auto nneigh = grid->neighbor_north();
       recvbuf = new double[grid_nx*7*2];
       err = MPI_Irecv(&recvbuf[0], grid_nx*7*2, MPI_DOUBLE, nneigh, 0, partition->comm(), &request);
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
    }

     // Send south
     if (ny0 > 0) {
       sendbuf = new double[grid_nx*7*2];
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          for (int j=0; j<7; j++) {
            sendbuf[ind] = momentMat(i,1,j);
            ind++;
            sendbuf[ind] = momentMat(i,0,j);
            ind++;
          }
       }
       err = MPI_Send(&sendbuf[0], grid_nx*7*2, MPI_DOUBLE, sneigh, 0, partition->comm());
       MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
       log->check(err == MPI_SUCCESS, ": ", (std::string) mpiErrBuffer);
      }

    if (ny1 > 0) {
       MPI_Status status;
       MPI_Wait(&request,&status);
       int ind = 0;
       for (int i=0; i < grid_nx; i++){
          for (int j=0; j<7; j++) {
            momentMat(i,grid_ny-ny1,j) = recvbuf[ind];
            ind++;
            momentMat(i,grid_ny-ny1-1,j) = recvbuf[ind];
            ind++;
          }
       }
    }
  }

  /*===========================================================================================
   * Locate particle on structured grid
   *===========================================================================================*/
  KOKKOS_INLINE_FUNCTION
  void Interpolation::locate_particle(double x, double y, double grid_dx, double grid_dy,
                                      int & iglobal, int & jglobal, double & xi, double & eta) {

       // grid resolution
       double dx = grid_dx;
       double dy = grid_dy;
       
       // global index of cell where particle is located
       iglobal = std::floor( x / dx );
       jglobal = std::floor( y / dy );

       // reference coordinates for [-1,1]x[-1,1] cell
       double x1 = dx*iglobal + dx;
       double y1 = dy*jglobal + dy;
       xi  = 1.0 - 2.0*(x1 - x)/dx;
       eta = 1.0 - 2.0*(y1 - y)/dy;

  }

  /*===========================================================================================
   * Evaluate basis function at particle location
   *===========================================================================================*/
  KOKKOS_INLINE_FUNCTION
  void Interpolation::grid_basis(double xi, double eta, double basis[4]) {

      // bilinear basis functions
       basis[0] = 0.25*(1.0-xi)*(1.0-eta);
       basis[1] = 0.25*(1.0+xi)*(1.0-eta);
       basis[2] = 0.25*(1.0+xi)*(1.0+eta);
       basis[3] = 0.25*(1.0-xi)*(1.0+eta);
  }

} // namespace DEMSI

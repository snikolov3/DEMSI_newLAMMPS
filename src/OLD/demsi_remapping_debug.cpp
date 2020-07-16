#include "demsi_remapping.h"
#include "demsi_communication.h"
#include "demsi_column.h"

#include <netcdf.h>

namespace DEMSI {

//------------------------------------------------------------------------------
// Remapping routines for debugging
//------------------------------------------------------------------------------

void Remapping::output_init_particle_array(const std::string filenameOut, double array[]) {

  int err, ncID;
  int dimidParticles, dimidMaxVertices;
  int varidNVerticesOnCell, varidType, varidX, varidY, varidArray;

  err = nc_create(filenameOut.c_str(), NC_CLOBBER, &ncID);
  log->check(err == NC_NOERR, "Problem creating output file: ", nc_strerror(err), ", for ", filenameOut);

  err = nc_def_dim(ncID, "nParticles", tessellation->nParticlesAllInit, &dimidParticles);
  log->check(err == NC_NOERR, "Problem defining nParticles dimension: ", nc_strerror(err), ", for ", filenameOut);

  err = nc_def_dim(ncID, "maxVertices", tessellation->maxVerticesInit, &dimidMaxVertices);
  log->check(err == NC_NOERR, "Problem defining nParticles dimension: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsNVerticesOnCell[1] = {dimidParticles};
  err = nc_def_var(ncID, "nVerticesOnCell", NC_INT, 1, dimidsNVerticesOnCell, &varidNVerticesOnCell);
  log->check(err == NC_NOERR, "Problem defining nVerticesOnCell variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsType[1] = {dimidParticles};
  err = nc_def_var(ncID, "type", NC_INT, 1, dimidsType, &varidType);
  log->check(err == NC_NOERR, "Problem defining type variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsX[2] = {dimidParticles,dimidMaxVertices};
  err = nc_def_var(ncID, "xVertex", NC_DOUBLE, 2, dimidsX, &varidX);
  log->check(err == NC_NOERR, "Problem defining x variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsY[2] = {dimidParticles,dimidMaxVertices};
  err = nc_def_var(ncID, "yVertex", NC_DOUBLE, 2, dimidsY, &varidY);
  log->check(err == NC_NOERR, "Problem defining y variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsArray[1] = {dimidParticles};
  err = nc_def_var(ncID, "array", NC_DOUBLE, 1, dimidsArray, &varidArray);
  log->check(err == NC_NOERR, "Problem defining array variable: ", nc_strerror(err), ", for ", filenameOut);

  err = nc_enddef(ncID);
  log->check(err == NC_NOERR, "Problem with enddef: ", nc_strerror(err), ", for ", filenameOut);

  size_t start1D[1];
  size_t count1D[1];
  start1D[0] = 0;
  count1D[0] = (size_t) tessellation->nParticlesAllInit;

  size_t start2D[2];
  size_t count2D[2];
  start2D[0] = 0;
  count2D[0] = (size_t) tessellation->nParticlesAllInit;
  start2D[1] = 0;
  count2D[1] = (size_t) tessellation->maxVerticesInit;

  double xx[tessellation->nParticlesAllInit * tessellation->maxVerticesInit];
  double yy[tessellation->nParticlesAllInit * tessellation->maxVerticesInit];

  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    for (int iVertex = 0 ; iVertex < tessellation->maxVerticesInit ; iVertex++) {
      int ij = tessellation->maxVerticesInit * iParticleInit + iVertex;
      xx[ij] = tessellation->xVertexInit[iParticleInit][iVertex];
      yy[ij] = tessellation->yVertexInit[iParticleInit][iVertex];
    } // iVertex
  } // iParticle

  err = nc_put_vara_int(ncID, varidType, start1D, count1D, &tessellation->typeInit[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable type for file ", filenameOut);

  err = nc_put_vara_int(ncID, varidNVerticesOnCell, start1D, count1D, &tessellation->nVerticesOnCellInit[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable nVerticesOnCell for file ", filenameOut);

  err = nc_put_vara_double(ncID, varidX, start2D, count2D, &xx[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable x for file ", filenameOut);

  err = nc_put_vara_double(ncID, varidY, start2D, count2D, &yy[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable y for file ", filenameOut);

  err = nc_put_vara_double(ncID, varidArray, start1D, count1D, &array[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable array for file ", filenameOut);

  err = nc_close(ncID);
  log->check(err == NC_NOERR, "Problem closing file: ", nc_strerror(err), ", for file ");

} // Remapping::output_init_particle_array

void Remapping::output_init_particle_array(const std::string filenameOut, int array[]) {

  int err, ncID;
  int dimidParticles, dimidMaxVertices;
  int varidNVerticesOnCell, varidType, varidX, varidY, varidArray;

  err = nc_create(filenameOut.c_str(), NC_CLOBBER, &ncID);
  log->check(err == NC_NOERR, "Problem creating output file: ", nc_strerror(err), ", for ", filenameOut);

  err = nc_def_dim(ncID, "nParticles", tessellation->nParticlesAllInit, &dimidParticles);
  log->check(err == NC_NOERR, "Problem defining nParticles dimension: ", nc_strerror(err), ", for ", filenameOut);

  err = nc_def_dim(ncID, "maxVertices", tessellation->maxVerticesInit, &dimidMaxVertices);
  log->check(err == NC_NOERR, "Problem defining nParticles dimension: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsNVerticesOnCell[1] = {dimidParticles};
  err = nc_def_var(ncID, "nVerticesOnCell", NC_INT, 1, dimidsNVerticesOnCell, &varidNVerticesOnCell);
  log->check(err == NC_NOERR, "Problem defining nVerticesOnCell variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsType[1] = {dimidParticles};
  err = nc_def_var(ncID, "type", NC_INT, 1, dimidsType, &varidType);
  log->check(err == NC_NOERR, "Problem defining type variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsX[2] = {dimidParticles,dimidMaxVertices};
  err = nc_def_var(ncID, "xVertex", NC_DOUBLE, 2, dimidsX, &varidX);
  log->check(err == NC_NOERR, "Problem defining x variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsY[2] = {dimidParticles,dimidMaxVertices};
  err = nc_def_var(ncID, "yVertex", NC_DOUBLE, 2, dimidsY, &varidY);
  log->check(err == NC_NOERR, "Problem defining y variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsArray[1] = {dimidParticles};
  err = nc_def_var(ncID, "array", NC_INT, 1, dimidsArray, &varidArray);
  log->check(err == NC_NOERR, "Problem defining array variable: ", nc_strerror(err), ", for ", filenameOut);

  err = nc_enddef(ncID);
  log->check(err == NC_NOERR, "Problem with enddef: ", nc_strerror(err), ", for ", filenameOut);

  size_t start1D[1];
  size_t count1D[1];
  start1D[0] = 0;
  count1D[0] = (size_t) tessellation->nParticlesAllInit;

  size_t start2D[2];
  size_t count2D[2];
  start2D[0] = 0;
  count2D[0] = (size_t) tessellation->nParticlesAllInit;
  start2D[1] = 0;
  count2D[1] = (size_t) tessellation->maxVerticesInit;

  double xx[tessellation->nParticlesAllInit * tessellation->maxVerticesInit];
  double yy[tessellation->nParticlesAllInit * tessellation->maxVerticesInit];

  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    for (int iVertex = 0 ; iVertex < tessellation->maxVerticesInit ; iVertex++) {
      int ij = tessellation->maxVerticesInit * iParticleInit + iVertex;
      xx[ij] = tessellation->xVertexInit[iParticleInit][iVertex];
      yy[ij] = tessellation->yVertexInit[iParticleInit][iVertex];
    } // iVertex
  } // iParticle

  err = nc_put_vara_int(ncID, varidType, start1D, count1D, &tessellation->typeInit[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable type for file ", filenameOut);

  err = nc_put_vara_int(ncID, varidNVerticesOnCell, start1D, count1D, &tessellation->nVerticesOnCellInit[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable nVerticesOnCell for file ", filenameOut);

  err = nc_put_vara_double(ncID, varidX, start2D, count2D, &xx[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable x for file ", filenameOut);

  err = nc_put_vara_double(ncID, varidY, start2D, count2D, &yy[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable y for file ", filenameOut);

  err = nc_put_vara_int(ncID, varidArray, start1D, count1D, &array[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable array for file ", filenameOut);

  err = nc_close(ncID);
  log->check(err == NC_NOERR, "Problem closing file: ", nc_strerror(err), ", for file ");

} // Remapping::output_init_particle_array

// output areas
void Remapping::output_areas(void) {

  int err, ncID;
  int dimidParticles;
  int varidType, varidX, varidY, varidR, varidAreaNew, varidAreaOld;

  err = nc_create("areas_remap.nc", NC_CLOBBER, &ncID);
  log->check(err == NC_NOERR, "Problem creating output file: ", nc_strerror(err), ", for output_areas");

  err = nc_def_dim(ncID, "nParticles", *(particles->nParticles), &dimidParticles);
  log->check(err == NC_NOERR, "Problem defining nParticles dimension: ", nc_strerror(err), ", for output_areas");

  int dimidsType[1] = {dimidParticles};
  err = nc_def_var(ncID, "type", NC_INT, 1, dimidsType, &varidType);
  log->check(err == NC_NOERR, "Problem defining type variable: ", nc_strerror(err), ", for output_areas");

  int dimidsX[1] = {dimidParticles};
  err = nc_def_var(ncID, "x", NC_DOUBLE, 1, dimidsX, &varidX);
  log->check(err == NC_NOERR, "Problem defining x variable: ", nc_strerror(err), ", for output_areas");

  int dimidsY[1] = {dimidParticles};
  err = nc_def_var(ncID, "y", NC_DOUBLE, 1, dimidsY, &varidY);
  log->check(err == NC_NOERR, "Problem defining y variable: ", nc_strerror(err), ", for output_areas");

  int dimidsR[1] = {dimidParticles};
  err = nc_def_var(ncID, "radius", NC_DOUBLE, 1, dimidsR, &varidR);
  log->check(err == NC_NOERR, "Problem defining radius variable: ", nc_strerror(err), ", for output_areas");

  int dimidsAreaNew[1] = {dimidParticles};
  err = nc_def_var(ncID, "areaNew", NC_DOUBLE, 1, dimidsAreaNew, &varidAreaNew);
  log->check(err == NC_NOERR, "Problem defining areaNew variable: ", nc_strerror(err), ", for output_areas");

  int dimidsAreaOld[1] = {dimidParticles};
  err = nc_def_var(ncID, "areaOld", NC_DOUBLE, 1, dimidsAreaOld, &varidAreaOld);
  log->check(err == NC_NOERR, "Problem defining areaOld variable: ", nc_strerror(err), ", for output_areas");

  err = nc_enddef(ncID);
  log->check(err == NC_NOERR, "Problem with enddef: ", nc_strerror(err), ", for output_areas");

  size_t start[1];
  size_t count[1];
  start[0] = 0;
  count[0] = (size_t) *(particles->nParticles);

  int tt[*(particles->nParticles)];
  double xx[*(particles->nParticles)];
  double yy[*(particles->nParticles)];
  double rr[*(particles->nParticles)];
  double areaNew[*(particles->nParticles)];
  double areaOld[*(particles->nParticles)];

  for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {
    tt[iParticle] = particles->type(iParticle);
    xx[iParticle] = particles->x(iParticle,0);
    yy[iParticle] = particles->x(iParticle,1);
    rr[iParticle] = particles->radius(iParticle);
    areaNew[iParticle] = (*column->effectiveElementArea)(iParticle);
    areaOld[iParticle] = (*column->areaInit)(iParticle);
  } // iParticle

  err = nc_put_vara_int(ncID, varidType, start, count, &tt[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable type");

  err = nc_put_vara_double(ncID, varidX, start, count, &xx[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable x");

  err = nc_put_vara_double(ncID, varidY, start, count, &yy[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable y");

  err = nc_put_vara_double(ncID, varidR, start, count, &rr[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable radius");

  err = nc_put_vara_double(ncID, varidAreaNew, start, count, &areaNew[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable areaNew");

  err = nc_put_vara_double(ncID, varidAreaOld, start, count, &areaOld[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable areaOld");

  err = nc_close(ncID);
  log->check(err == NC_NOERR, "Problem closing file: ", nc_strerror(err), ", for file output_areas");

  exit(0);

} // Remapping::output_areas

// output column variable mid-remapping
void Remapping::output_post_remapping(const std::string filenameOut, DEMSI::ColumnVariable<double>* columnTracer, const int nParticlesNew, const std::vector<int> iParticleInitFromNew) {

  int err, ncID;
  int dimidParticles;
  int varidT, varidX, varidY, varidR, varidArray;

  err = nc_create(filenameOut.c_str(), NC_CLOBBER, &ncID);
  log->check(err == NC_NOERR, "Problem creating output file: ", nc_strerror(err), ", for ", filenameOut);

  err = nc_def_dim(ncID, "nParticles", nParticlesNew, &dimidParticles);
  log->check(err == NC_NOERR, "Problem defining nParticles dimension: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsT[1] = {dimidParticles};
  err = nc_def_var(ncID, "type", NC_INT, 1, dimidsT, &varidT);
  log->check(err == NC_NOERR, "Problem defining type variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsX[1] = {dimidParticles};
  err = nc_def_var(ncID, "x", NC_DOUBLE, 1, dimidsX, &varidX);
  log->check(err == NC_NOERR, "Problem defining x variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsY[1] = {dimidParticles};
  err = nc_def_var(ncID, "y", NC_DOUBLE, 1, dimidsY, &varidY);
  log->check(err == NC_NOERR, "Problem defining y variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsR[1] = {dimidParticles};
  err = nc_def_var(ncID, "radius", NC_DOUBLE, 1, dimidsR, &varidR);
  log->check(err == NC_NOERR, "Problem defining radius variable: ", nc_strerror(err), ", for ", filenameOut);

  int dimidsArray[1] = {dimidParticles};
  err = nc_def_var(ncID, "array", NC_DOUBLE, 1, dimidsArray, &varidArray);
  log->check(err == NC_NOERR, "Problem defining array variable: ", nc_strerror(err), ", for ", filenameOut);

  err = nc_enddef(ncID);
  log->check(err == NC_NOERR, "Problem with enddef: ", nc_strerror(err), ", for ", filenameOut);

  size_t start1D[1];
  size_t count1D[1];
  start1D[0] = 0;
  count1D[0] = (size_t) nParticlesNew;

  double xx[nParticlesNew];
  double yy[nParticlesNew];
  double rr[nParticlesNew];
  int tt[nParticlesNew];
  double array[nParticlesNew];

  for (int iParticle = 0 ; iParticle < nParticlesNew ; iParticle++) {
    int iParticleInit = iParticleInitFromNew[iParticle];
    xx[iParticle] = tessellation->xInit[iParticleInit];
    yy[iParticle] = tessellation->yInit[iParticleInit];
    rr[iParticle] = tessellation->radiusInit[iParticleInit];
    tt[iParticle] = tessellation->typeInit[iParticleInit];
    if (tt[iParticle] == 0) tt[iParticle] = 1;
    array[iParticle] = (*columnTracer)(iParticle);
  } // iParticle

  err = nc_put_vara_int(ncID, varidT, start1D, count1D, &tt[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable type for file ", filenameOut);

  err = nc_put_vara_double(ncID, varidX, start1D, count1D, &xx[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable x for file ", filenameOut);

  err = nc_put_vara_double(ncID, varidY, start1D, count1D, &yy[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable y for file ", filenameOut);

  err = nc_put_vara_double(ncID, varidR, start1D, count1D, &rr[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable radius for file ", filenameOut);

  err = nc_put_vara_double(ncID, varidArray, start1D, count1D, &array[0]);
  log->check(err == NC_NOERR, "Problem writing variable: ", nc_strerror(err), ", for variable array for file ", filenameOut);

  err = nc_close(ncID);
  log->check(err == NC_NOERR, "Problem closing file: ", nc_strerror(err), ", for file ");

} // Remapping::output_post_remapping

} // namespace DEMSI

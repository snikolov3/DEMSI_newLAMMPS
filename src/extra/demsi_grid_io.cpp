#include "demsi_grid_io.h"

#include "demsi_grid.h"
#include "demsi_logging.h"
#include "demsi_file_utils.h"

#include <netcdf.h>

namespace DEMSI {

//------------------------------------------------------------------------------
// class GridWrite
//------------------------------------------------------------------------------

GridWrite::GridWrite(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Grid* gridIn, const std::string filenameTemplateIn, DEMSI::Clock* clockIn, DEMSI::Time alarmStartTime, DEMSI::TimeInterval alarmInterval, const bool clobberIn, std::vector<std::string> outputFieldNamesIn) {

  partition = partitionIn;
  log = logIn;
  grid = gridIn;
  filenameTemplate = filenameTemplateIn;
  clobber = clobberIn;
  outputFieldNames = outputFieldNamesIn;

  // check if filename is none or variant
  std::string filenameTemplateUpper = filenameTemplate;
  std::transform(filenameTemplateUpper.begin(), filenameTemplateUpper.end(), filenameTemplateUpper.begin(), toupper);
  if (filenameTemplateUpper == "NONE") {
    active = false;
  } else {
    active = true;
    outputAlarm.set(clockIn, alarmStartTime, alarmInterval);
  }

  clock = clockIn;

}

// Write a single Eulerian field to a netcdf file.
void GridWrite::write_field(Kokkos::View<double**>* field, const int nc_id, const size_t iTime) {

  int err;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int sendCount = grid->ownedDimensions[0]*grid->ownedDimensions[1];
  double* sendBuffer = new double [grid->ownedDimensions[0]*grid->ownedDimensions[1]];

  int ij = 0;
  for (int i=0 ; i<grid->ownedDimensions[0] ; i++) {
    for (int j=0 ; j<grid->ownedDimensions[1] ; j++) {

      int iGlobal = grid->gridPartitions[partition->this_proc()]->iglobal_owned(i);
      int jGlobal = grid->gridPartitions[partition->this_proc()]->jglobal_owned(j);

      int iLocalTot = grid->gridPartitions[partition->this_proc()]->ilocal_total(iGlobal);
      int jLocalTot = grid->gridPartitions[partition->this_proc()]->jlocal_total(jGlobal);

      sendBuffer[ij] = (*field)(iLocalTot,jLocalTot);
      ij++;
    }
  }

  int* recvCount;
  int* displacements;
  double* recvBuffer;
  int recvCountTotal;

  if (partition->on_master_proc()) {

    recvCount = new int [partition->nprocs()];
    displacements = new int [partition->nprocs()];

    recvCountTotal = 0;
    for (int iProc = 0 ; iProc<partition->nprocs() ; iProc++) {
      recvCount[iProc] = grid->gridPartitions[iProc]->size_owned();
      recvCountTotal = recvCountTotal + recvCount[iProc];
    }
    displacements[0] = 0;
    for (int iProc = 1 ; iProc<partition->nprocs() ; iProc++) {
      displacements[iProc] = displacements[iProc-1] + recvCount[iProc-1];
    }

    recvBuffer = new double [recvCountTotal];

  }

  err = MPI_Gatherv(sendBuffer, sendCount, MPI_DOUBLE, recvBuffer, recvCount, displacements, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "GridWrite::write_field: Problem gathering field: ", (std::string) mpiErrBuffer);

  if (partition->on_master_proc()) {

    int varid;

    std::string fieldName = field->label();
    err = nc_inq_varid(nc_id, fieldName.c_str(), &varid);
    log->check(err == NC_NOERR, "GridWrite::write_field: Problem getting field varid: ", fieldName, ": ", nc_strerror(err));

    double* arrayOut = new double [grid->globalDimensions[0]*grid->globalDimensions[1]];

    int ij = 0;
    for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {

      for (int i=0 ; i<grid->gridPartitions[iProc]->size_owned(0) ; i++) {
	for (int j=0 ; j<grid->gridPartitions[iProc]->size_owned(1) ; j++) {

	  int ijGlobal = grid->gridPartitions[iProc]->iglobal_owned(i) * grid->globalDimensions[1] + grid->gridPartitions[iProc]->jglobal_owned(j);
	  arrayOut[ijGlobal] = recvBuffer[ij];
	  ij++;

	} // j
      } // i

    } // iProc

    size_t start[3] = {iTime, 0, 0};
    size_t count[3] = {1, (size_t) grid->globalDimensions[0], (size_t) grid->globalDimensions[1]};

    err = nc_put_vara_double(nc_id, varid, start, count, &arrayOut[0]);
    log->check(err == NC_NOERR, "GridWrite::write_field: Problem writing field: ", nc_strerror(err), ", for field: ", fieldName);

    delete [] arrayOut;
    delete [] recvCount;
    delete [] displacements;
    delete [] recvBuffer;

  } // master proc

  delete [] sendBuffer;

}

void GridWrite::write(void) {

  int nc_id;

  // get the current time from the simulation clock
  Time writeTime = clock->get_time();

  // check if we should write an output
  if (active and outputAlarm.is_ringing()) {

    int err;
    char mpiErrBuffer[MPI_MAX_ERROR_STRING];
    int mpiErrLen;

    int dimidTime, dimidStrLen, dimidnx, dimidny, dimidTWO;
    size_t nTimes;

    if (partition->on_master_proc()) {

      // expand filename template
      std::string filenameOut = expand_filename_template(filenameTemplate, writeTime);

      // try to create new file
      err = nc_create(filenameOut.c_str(), NC_NOCLOBBER, &nc_id);
      if (err == NC_NOERR) {

	nTimes = 0;

	// define dimensions
	err = nc_def_dim(nc_id, "Time", NC_UNLIMITED, &dimidTime);
	log->check(err == NC_NOERR, "GridWrite::write: Problem defining Time dimension: ", nc_strerror(err));

	err = nc_def_dim(nc_id, "strLen", 64, &dimidStrLen);
	log->check(err == NC_NOERR, "GridWrite::write: Problem defining strLen dimension: ", nc_strerror(err));

	err = nc_def_dim(nc_id, "nx", grid->globalDimensions[0], &dimidnx);
	log->check(err == NC_NOERR, "GridWrite::write: Problem defining nx dimension: ", nc_strerror(err));

	err = nc_def_dim(nc_id, "ny", grid->globalDimensions[1], &dimidny);
	log->check(err == NC_NOERR, "GridWrite::write: Problem defining ny dimension: ", nc_strerror(err));

	err = nc_def_dim(nc_id, "TWO", 2, &dimidTWO);
	log->check(err == NC_NOERR, "GridWrite::write: Problem defining TWO dimension: ", nc_strerror(err));

	int varidTime;
	int dimidsTime[2] = {dimidTime,dimidStrLen};
	err = nc_def_var(nc_id, "Time", NC_CHAR, 2, dimidsTime, &varidTime);
	log->check(err == NC_NOERR, "GridWrite::write: Problem defining Time variable: ", nc_strerror(err), ", for file: ", filenameOut);

	int varidResolution;
	int dimidsResolution[1] = {dimidTWO};
	err = nc_def_var(nc_id, "resolution", NC_DOUBLE, 1, dimidsResolution, &varidResolution);
	log->check(err == NC_NOERR, "GridWrite::write: Problem defining resolution variable: ", nc_strerror(err));

	for (int iFieldOut = 0; iFieldOut < outputFieldNames.size() ; iFieldOut++) {
	  for (int iField = 0; iField < grid->fieldsGridWrite.size() ; iField++) {
	    if (outputFieldNames[iFieldOut] == grid->fieldsGridWrite[iField]->label()) {

	      int dimids[3] = {dimidTime, dimidnx, dimidny};
	      int varid;
	      std::string fieldName = grid->fieldsGridWrite[iField]->label();
	      err = nc_def_var(nc_id, fieldName.c_str(), NC_DOUBLE, 3, dimids, &varid);
	      log->check(err == NC_NOERR, "GridWrite::write: Problem defining field: ", fieldName, ": ", nc_strerror(err));

	    }
	  }
	}

	// end definition phase
	err = nc_enddef(nc_id);
	log->check(err == NC_NOERR, "GridWrite::write: Problem endding definition phase: ", nc_strerror(err));

	// write the time
	char* strOut = writeTime.get_char64();

	size_t startTime[2] = {nTimes, 0};
	size_t countTime[2] = {1, strlen(strOut)+1};

	err = nc_put_vara_text(nc_id, varidTime, startTime, countTime, &strOut[0]);
	log->check(err == NC_NOERR, "GridWrite::write: Problem writing time: ", nc_strerror(err));
	delete [] strOut;

	// resolution
	size_t start[1] = {0};
	size_t count[1] = {2};
	err = nc_put_vara_double(nc_id, varidResolution, start, count, &(grid->resolution[0]));
	log->check(err == NC_NOERR, "GridWrite::write: Problem writing resolution: ", nc_strerror(err));

      } else if (err == NC_EEXIST) {

	// open preexisting file
	err = nc_open(filenameOut.c_str(), NC_WRITE, &nc_id);
	log->check(err == NC_NOERR, "GridWrite::write: Problem creating output eulerian file: ", filenameOut, ", :", nc_strerror(err));

	// get current time
	err = nc_inq_dimid(nc_id, "Time", &dimidTime);
	log->check(err == NC_NOERR, "GridWrite::write: Problem defining Time dimension: ", nc_strerror(err));

	err = nc_inq_dimlen(nc_id, dimidTime, &nTimes);
	log->check(err == NC_NOERR, "GridWrite::write: Problem getting Time dimension: ", nc_strerror(err));

	int varidTime;
	err = nc_inq_varid(nc_id, "Time", &varidTime);
	log->check(err == NC_NOERR, "GridWrite::write: Problem getting time varid: ", nc_strerror(err));

	// write the time
	char* strOut = writeTime.get_char64();

	size_t startTime[2] = {nTimes, 0};
	size_t countTime[2] = {1, strlen(strOut)+1};

	err = nc_put_vara_text(nc_id, varidTime, startTime, countTime, &strOut[0]);
	log->check(err == NC_NOERR, "GridWrite::write: Problem writing time: ", nc_strerror(err));
	delete [] strOut;

      } else {
	// failed to create file
	log->abort("GridWrite::write: Failed to create grid output file.");
      }

    } // on master proc

    for (int iFieldOut = 0; iFieldOut < outputFieldNames.size() ; iFieldOut++) {
      for (int iField = 0; iField < grid->fieldsGridWrite.size() ; iField++) {
	if (outputFieldNames[iFieldOut] == grid->fieldsGridWrite[iField]->label()) {

	  this->write_field(grid->fieldsGridWrite[iField], nc_id, nTimes);

	}
      }
    } // field loop

    if (partition->on_master_proc()) {

      int err = nc_close(nc_id);
      log->check(err == NC_NOERR, "GridWrite::write: Problem closing file");

    } // master proc

    // reset output alarm
    outputAlarm.reset();

  } // alarm ringing

}

//------------------------------------------------------------------------------
// class GridOutputStreams
//------------------------------------------------------------------------------

// Constructor for the ParticlesOutputStreams class.
  GridOutputStreams::GridOutputStreams(DEMSI::Partition* partition, DEMSI::Log* log, DEMSI::Configs* configs, DEMSI::Grid* grid, DEMSI::Clock* clock, DEMSI::Time alarmStartTime) {

  // check that have gridded output
  if (configs->exists({"ConfigGroup:griddedOutput","Array:griddedOutputStreams"})) {

    // get gridded output stream names
    std::vector<std::string> outputStreamNames = configs->get_array({"ConfigGroup:griddedOutput","Array:griddedOutputStreams"}, "GriddedOutputStream", "name");

    // iterate over output streams
    for (int iStream=0 ; iStream<outputStreamNames.size() ; iStream++) {

      // output stream configs location
      std::list<std::string> configLocationRoot, configLocation;
      configLocationRoot.push_back("ConfigGroup:griddedOutput");
      configLocationRoot.push_back("Array:griddedOutputStreams");
      configLocationRoot.push_back("GriddedOutputStream:" + outputStreamNames[iStream]);

      // get output stream configs
      std::string writeFilenameTemplate;
      configLocation = configLocationRoot;
      configLocation.push_back("Config:griddedWriteFilenameTemplate");
      configs->get(configLocation, writeFilenameTemplate);

      std::string writeOutputDirectory;
      configLocation = configLocationRoot;
      configLocation.push_back("Config:griddedWriteOutputDirectory");
      configs->get(configLocation, writeOutputDirectory);

      writeFilenameTemplate = writeOutputDirectory + writeFilenameTemplate;

      std::string writeIntervalStr;
      configLocation = configLocationRoot;
      configLocation.push_back("Config:griddedWriteInterval");
      configs->get(configLocation, writeIntervalStr);

      DEMSI::TimeInterval writeInterval(writeIntervalStr, log);

      bool writeClobber;
      configLocation = configLocationRoot;
      configLocation.push_back("Config:griddedWriteClobber");
      configs->get(configLocation, writeClobber);

      // get output field names
      std::list<std::string> outputFieldsLocation = configLocationRoot;
      outputFieldsLocation.push_back("Array:OutputFields");
      std::vector<std::string> outputFieldNames = configs->get_array(outputFieldsLocation, "OutputField");

      // create a new GriddedsWrite object
      DEMSI::GridWrite* gridWrite = new DEMSI::GridWrite(partition, log, grid, writeFilenameTemplate, clock, alarmStartTime, writeInterval, writeClobber, outputFieldNames);
      streams.push_back(gridWrite);

    } // iStream

  } // gridded output exists

}

// Write out all gridded output streams
void GridOutputStreams::write(void) {

  for (int iStream = 0 ; iStream<streams.size() ; iStream++) {
    streams[iStream]->write();
  }

}

} // namespace DEMSI

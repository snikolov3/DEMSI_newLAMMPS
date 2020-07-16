#include "demsi_initial_tessellation_io.h"

#include <netcdf.h>
#include "demsi_file_utils.h"

namespace DEMSI {

//------------------------------------------------------------------------------
// class TessellationWrite
//------------------------------------------------------------------------------

TessellationWrite::TessellationWrite(DEMSI::Partition* partitionIn, DEMSI::Log* logIn, DEMSI::Tessellation* tessellationIn, const std::string streamNameIn, const std::string filenameTemplateIn, DEMSI::Clock* clockIn, DEMSI::Time alarmStartTime, const std::string writeIntervalStr, const bool clobberIn, std::vector<std::string> outputFieldNamesIn) {

  partition = partitionIn;
  log = logIn;
  tessellation = tessellationIn;
  streamName = streamNameIn;
  filenameTemplate = filenameTemplateIn;
  clobber = clobberIn;
  outputFieldNames = outputFieldNamesIn;
  clock = clockIn;

  // check if alarm interval is none
  if (writeIntervalStr == "none") {
    active = false;
  } else {
    active = true;
    DEMSI::TimeInterval writeInterval(writeIntervalStr, log);
    outputAlarm.set(clock, alarmStartTime, writeInterval);
  }

} // TessellationWrite::TessellationWrite

// Write a single Eulerian field to a netcdf file.
void TessellationWrite::write_field(std::pair<double*,std::string> field, const int nc_id, const size_t iTime, int nParticlesInitAllProcs, int* nParticlesInitOnProcs, int* displacements) {

  int err;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int sendCount = tessellation->nParticlesInit;
  double* sendBuffer = new double[sendCount];

  for (int i=0 ; i<sendCount ; i++) {
    sendBuffer[i] = field.first[i];
  } // i

  double* recvBuffer;
  if (partition->on_master_proc()) {
    recvBuffer = new double[nParticlesInitAllProcs];
  }

  err = MPI_Gatherv(sendBuffer, sendCount, MPI_DOUBLE, recvBuffer, nParticlesInitOnProcs, displacements, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "TessellationWrite::write_field: Problem gathering field: ", (std::string) mpiErrBuffer);

  if (partition->on_master_proc()) {

    int varid;

    std::string fieldName = field.second;
    err = nc_inq_varid(nc_id, fieldName.c_str(), &varid);
    log->check(err == NC_NOERR, "TessellationWrite::write_field: Problem getting field varid: ", nc_strerror(err), ", for field: ", fieldName);

    size_t start[2] = {iTime, 0};
    size_t count[2] = {1, (size_t) nParticlesInitAllProcs};

    err = nc_put_vara_double(nc_id, varid, start, count, &recvBuffer[0]);
    log->check(err == NC_NOERR, "TessellationWrite::write_field: Problem writing field: ", nc_strerror(err), ", for field: ", fieldName);

    delete [] recvBuffer;

  } // master proc

  delete [] sendBuffer;

} // TessellationWrite::write_field

// Write the globalIndexInit field to a netcdf file.
void TessellationWrite::write_globalIndexInit_field(const int nc_id, int nParticlesInitAllProcs, int* nParticlesInitOnProcs, int* displacements) {

  int err;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  int sendCount = tessellation->nParticlesInit;
  int* sendBuffer = new int[sendCount];

  for (int i=0 ; i<sendCount ; i++) {
    sendBuffer[i] = tessellation->globalIndexInit[i];
  } // i

  int* recvBuffer;
  if (partition->on_master_proc()) {
    recvBuffer = new int[nParticlesInitAllProcs];
  }
  (*log)(DEMSI::LOG::DEBUG) << nParticlesInitAllProcs << " " << sendCount << std::endl;
  for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
    (*log)(DEMSI::LOG::DEBUG) << iProc << " " << nParticlesInitOnProcs[iProc] << " " << displacements[iProc] << std::endl;
  }
  err = MPI_Gatherv(sendBuffer, sendCount, MPI_INT, recvBuffer, nParticlesInitOnProcs, displacements, MPI_INT, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "TessellationWrite::write_field: Problem gathering field globalIDInit: ", (std::string) mpiErrBuffer);

  if (partition->on_master_proc()) {

    int varid;

    err = nc_inq_varid(nc_id, "globalIndexInit", &varid);
    log->check(err == NC_NOERR, "TessellationWrite::write_field: Problem getting field varid: ", nc_strerror(err), ", for field: globalIndexInit");

    size_t start[1] = {0};
    size_t count[1] = {(size_t) nParticlesInitAllProcs};

    err = nc_put_vara_int(nc_id, varid, start, count, &recvBuffer[0]);
    log->check(err == NC_NOERR, "TessellationWrite::write_field: Problem writing field: ", nc_strerror(err), ", for field: globalIndexInit");

    delete [] recvBuffer;

  } // master proc

  delete [] sendBuffer;

} // TessellationWrite::write_globalIndexInit_field

void TessellationWrite::write(const bool writeNow = false) {

  int nc_id;

  // get the current time from the simulation clock
  Time writeTime = clock->get_time();

  // check if we should write an output
  if ((active and outputAlarm.is_ringing()) or writeNow) {

    int err;
    char mpiErrBuffer[MPI_MAX_ERROR_STRING];
    int mpiErrLen;

    // displacements
    int nParticlesInitOnProcs[partition->nprocs()];
    int displacements[partition->nprocs()];
    int nParticlesInitAllProcs = 0;

    err = MPI_Allgather(&tessellation->nParticlesInit, 1, MPI_INT, nParticlesInitOnProcs, 1, MPI_INT, partition->comm());
    MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
    log->check(err == MPI_SUCCESS, "TessellationWrite::write: MPI_Allgather failed: ", (std::string) mpiErrBuffer);

    displacements[0] = 0;
    for (int iProc = 1 ; iProc < partition->nprocs() ; iProc++) {
      displacements[iProc] = displacements[iProc-1] + nParticlesInitOnProcs[iProc-1];
    } // iProc

    for (int iProc = 0 ; iProc < partition->nprocs() ; iProc++) {
      nParticlesInitAllProcs += nParticlesInitOnProcs[iProc];
    } // iProc

    int dimidTime, dimidStrLen, dimidnParticlesInit, dimidTWO;
    size_t nTimes;

    int newFile;

    if (partition->on_master_proc()) {

      // expand filename template
      std::string filenameOut = expand_filename_template(filenameTemplate, writeTime);

      // try to create new file
      err = nc_create(filenameOut.c_str(), NC_NOCLOBBER, &nc_id);
      if (err == NC_NOERR) {

	newFile = 1;
	nTimes = 0;

	// define dimensions
	err = nc_def_dim(nc_id, "Time", NC_UNLIMITED, &dimidTime);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem defining Time dimension: ", nc_strerror(err));

	err = nc_def_dim(nc_id, "strLen", 64, &dimidStrLen);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem defining strLen dimension: ", nc_strerror(err));

	err = nc_def_dim(nc_id, "nParticlesInit", nParticlesInitAllProcs, &dimidnParticlesInit);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem defining nParticlesInit dimension: ", nc_strerror(err));

	int varidTime;
	int dimidsTime[2] = {dimidTime,dimidStrLen};
	err = nc_def_var(nc_id, "Time", NC_CHAR, 2, dimidsTime, &varidTime);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem defining Time variable: ", nc_strerror(err), ", for file: ", filenameOut);

	// globalIndexInit
	int varidGlobalIndexInit;
	int dimidsGlobalIndexInit[1] = {dimidnParticlesInit};
	err = nc_def_var(nc_id, "globalIndexInit", NC_INT, 1, dimidsGlobalIndexInit, &varidGlobalIndexInit);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem defining field: globalIndexInit: ", nc_strerror(err));

	// output
	for (int iFieldOut = 0; iFieldOut < outputFieldNames.size() ; iFieldOut++) {
	  for (int iField = 0; iField < tessellation->fieldsWrite.size() ; iField++) {
	    if (outputFieldNames[iFieldOut] == tessellation->fieldsWrite[iField].second) {

	      int dimids[2] = {dimidTime, dimidnParticlesInit};
	      int varid;
	      std::string fieldName = tessellation->fieldsWrite[iField].second;
	      err = nc_def_var(nc_id, fieldName.c_str(), NC_DOUBLE, 2, dimids, &varid);
	      log->check(err == NC_NOERR, "TessellationWrite::write: Problem defining field: ", fieldName, ": ", nc_strerror(err));

	    }
	  }
	}

	// end definition phase
	err = nc_enddef(nc_id);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem endding definition phase: ", nc_strerror(err));

	// write the time
	char* strOut = writeTime.get_char64();

	size_t startTime[2] = {nTimes, 0};
	size_t countTime[2] = {1, strlen(strOut)+1};

	err = nc_put_vara_text(nc_id, varidTime, startTime, countTime, &strOut[0]);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem writing time: ", nc_strerror(err));
	delete [] strOut;

      } else if (err == NC_EEXIST) {

	newFile = 0;

	// open preexisting file
	err = nc_open(filenameOut.c_str(), NC_WRITE, &nc_id);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem creating output eulerian file: ", filenameOut, ", :", nc_strerror(err));

	// get current time
	err = nc_inq_dimid(nc_id, "Time", &dimidTime);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem defining Time dimension: ", nc_strerror(err));

	err = nc_inq_dimlen(nc_id, dimidTime, &nTimes);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem getting Time dimension: ", nc_strerror(err));

	int varidTime;
	err = nc_inq_varid(nc_id, "Time", &varidTime);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem getting time varid: ", nc_strerror(err));

	// write the time
	char* strOut = writeTime.get_char64();

	size_t startTime[2] = {nTimes, 0};
	size_t countTime[2] = {1, strlen(strOut)+1};

	err = nc_put_vara_text(nc_id, varidTime, startTime, countTime, &strOut[0]);
	log->check(err == NC_NOERR, "TessellationWrite::write: Problem writing time: ", nc_strerror(err));
	delete [] strOut;

      } else {
	// failed to create file
	log->abort("TessellationWrite::write: Failed to create tessellation output file:", filenameOut, ": ", nc_strerror(err));
      }

    } // on master proc

    // new file or not
    err = MPI_Bcast(&newFile, 1, MPI_INT, partition->master_proc(), partition->comm());
    MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
    log->check(err == MPI_SUCCESS, "TessellationWrite::write: Problem broadcasting newFile: ", (std::string) mpiErrBuffer);

    // write globalIndexInit
    if (newFile) this->write_globalIndexInit_field(nc_id, nParticlesInitAllProcs, nParticlesInitOnProcs, displacements);

    // output fields
    for (int iFieldOut = 0; iFieldOut < outputFieldNames.size() ; iFieldOut++) {
      for (int iField = 0; iField < tessellation->fieldsWrite.size() ; iField++) {
	if (outputFieldNames[iFieldOut] == tessellation->fieldsWrite[iField].second) {

	  this->write_field(tessellation->fieldsWrite[iField], nc_id, nTimes, nParticlesInitAllProcs, nParticlesInitOnProcs, displacements);

	}
      }
    } // field loop

    if (partition->on_master_proc()) {

      int err = nc_close(nc_id);
      log->check(err == NC_NOERR, "TessellationWrite::write: Problem closing file");

    } // master proc

    // reset output alarm
    if (active) outputAlarm.reset();

  } // alarm ringing

} // TessellationWrite::write

// Name of stream
std::string TessellationWrite::name(void) {
  return streamName;
} // TessellationWrite::name

//------------------------------------------------------------------------------
// class TessellationOutputStreams
//------------------------------------------------------------------------------

// Constructor for the ParticlesOutputStreams class.
TessellationOutputStreams::TessellationOutputStreams(DEMSI::Partition* partition, DEMSI::Log* log, DEMSI::Configs* configs, DEMSI::Tessellation* tessellation, DEMSI::Clock* clock, DEMSI::Time alarmStartTime) {

  // check that have output
  if (configs->exists({"ConfigGroup:tessellationOutput","Array:tessellationOutputStreams"})) {

    // get output stream names
    std::vector<std::string> outputStreamNames = configs->get_array({"ConfigGroup:tessellationOutput","Array:tessellationOutputStreams"}, "TessellationOutputStream", "name");

    // iterate over output streams
    for (int iStream=0 ; iStream<outputStreamNames.size() ; iStream++) {

      // stream name
      std::string streamName = outputStreamNames[iStream];

      // output stream configs location
      std::list<std::string> configLocationRoot, configLocation;
      configLocationRoot.push_back("ConfigGroup:tessellationOutput");
      configLocationRoot.push_back("Array:tessellationOutputStreams");
      configLocationRoot.push_back("TessellationOutputStream:" + streamName);

      // get output stream configs
      std::string writeFilenameTemplate;
      configLocation = configLocationRoot;
      configLocation.push_back("Config:tessellationWriteFilenameTemplate");
      configs->get(configLocation, writeFilenameTemplate);

      std::string writeOutputDirectory;
      configLocation = configLocationRoot;
      configLocation.push_back("Config:tessellationWriteOutputDirectory");
      configs->get(configLocation, writeOutputDirectory);

      writeFilenameTemplate = writeOutputDirectory + writeFilenameTemplate;

      std::string writeInterval;
      configLocation = configLocationRoot;
      configLocation.push_back("Config:tessellationWriteInterval");
      configs->get(configLocation, writeInterval);

      bool writeClobber;
      configLocation = configLocationRoot;
      configLocation.push_back("Config:tessellationWriteClobber");
      configs->get(configLocation, writeClobber);

      // get output field names
      std::list<std::string> outputFieldsLocation = configLocationRoot;
      outputFieldsLocation.push_back("Array:OutputFields");
      std::vector<std::string> outputFieldNames = configs->get_array(outputFieldsLocation, "OutputField");

      // create a new TessellationWrite object
      DEMSI::TessellationWrite* tessellationWrite = new DEMSI::TessellationWrite(partition, log, tessellation, streamName, writeFilenameTemplate, clock, alarmStartTime, writeInterval, writeClobber, outputFieldNames);
      streams.push_back(tessellationWrite);

    } // iStream

  } // output exists

} // TessellationOutputStreams::TessellationOutputStreams

// Write out all output streams
void TessellationOutputStreams::write(void) {

  for (int iStream = 0 ; iStream<streams.size() ; iStream++) {
    streams[iStream]->write();
  } // iStream

} // TessellationOutputStreams::write

// Write out particular output stream
void TessellationOutputStreams::write(const std::string streamName, const bool writeNow = false) {

  for (int iStream = 0 ; iStream<streams.size() ; iStream++) {
    if (streams[iStream]->name() == streamName) {
      streams[iStream]->write(writeNow);
    }
  } // iStream

} // TessellationOutputStreams::write

} // namespace DEMSI

#!/bin/bash -ex

source /etc/bashrc

export PATH=/usr/bin:/bin:/sbin:/usr/sbin:$PATH

source /ascldap/users/daibane/env-cee.sh
#needed to avoid conflicts with Dan B.'s custom Python which we use
module unload sems-python
module load sems-openmpi/1.8.7

TEST_DIR=/ascldap/users/daibane/demsi-nightly

if [ ! -d "$TEST_DIR" ]
then
  mkdir -p $TEST_DIR
fi

cd $TEST_DIR

which git
git --version
#grab the latest versions of all but this file
git -C /ascldap/users/daibane/src/DEMSI pull

LOG_FILE=$TEST_DIR/ctest_log.txt
CTEST_FILE=/ascldap/users/daibane/src/DEMSI/dashboards/cee-compute011.sandia.gov/nightly_ctest.cmake

ctest -VV -S $CTEST_FILE &> $LOG_FILE

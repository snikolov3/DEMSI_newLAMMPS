#ifndef DEMSI_COMMUNICATION_H_
#define DEMSI_COMMUNICATION_H_

/*! \file   communication.h
    \brief  Header file for the DEMSI communication classes and functions
*/

#include "demsi_logging.h"
#include "demsi_partition.h"

#include <vector>

namespace DEMSI {

/*! /brief Calculate the modulo of two numbers taking correctly into account negative numbers
    /param i Dividend
    /param n Divisor
*/
int modulo(int i, int n);

/*! \class ShareLists
  \brief This class allows lists to be send around the processors.

  This class allows lists defined on each processor to be sent around the
  processors in a wrapped order of processor ranks. For example, a list
  defined on processor i will be sent first to processor i+1, then i+2
  until the list ends up at processor i-1. Processor n (the last processor)
  send the list to processor 0. All processors do this simultaneously, each
  passing the list to the next processor at the same time. Between transfers
  the lists may be modified. An integer and double list can be send
  simultaneously. Usage:
  DEMSI::ShareLists* shareLists = new DEMSI::ShareLists(&intList,&doubleList,partition,log);
  while (shareLists->iterate()) {
     ...
  }
  delete shareLists;
*/
class ShareLists {
public:

  /*! \brief ShareLists constructor
      \param listIntegerIn Pointer to the integer list to send (or NULL).
      \param listDoubleIn Pointer to the double list to send (or NULL).
      \param partitionIn Pointer to partition object.
      \param logIn Pointer to log object.
   */
  ShareLists(std::vector<int>* listIntegerIn, std::vector<double>* listDoubleIn, DEMSI::Partition* partitionIn, DEMSI::Log* logIn);

  /*! \brief ShareLists default destructor */
  ~ShareLists() = default;

  /*! \brief Move lists to the next processor.
      \return True if lists could be moved. False if moved lists as far as possible.
  */
  bool iterate(void);

  /*! \brief Return processor ID where the current list originally came from.
      \return Processor ID where the current list originally came from.
  */
  int iProc_origin(void);

private:

  /*! Iterator for list sending */
  int iStep;

  /*! Processor ID where the current list originally came from. */
  int iProcTransit;

  /*! Integer list to send around processors */
  std::vector<int>* listInteger;

  /*! Double list to send around processors */
  std::vector<double>* listDouble;

  /*! Flag whether integer list is being sent */
  bool hasIntegerList;

  /*! Flag whether double list is being sent */
  bool hasDoubleList;

  /*! Pointer to partition object */
  DEMSI::Partition* partition;

  /*! Pointer to log object */
  DEMSI::Log* log;

}; // Column class

} // namespace DEMSI

#endif /* DEMSI_COMMUNICATION_H_ */

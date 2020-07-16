#ifndef DEMSI_COLUMN_VARIABLES_H_
#define DEMSI_COLUMN_VARIABLES_H_

/*! \file   demsi_column_variables.h
    \brief  Header file for the DEMSI::ColumnVariables, DEMSI::ColumnDimension, and DEMSI::ColumnDimensions classes
*/

#include <string>
#include <vector>
#include <set>
#include "demsi_configs.h"
#include "demsi_particles.h"
#include "demsi_partition.h"

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \class ColumnDimension
  \brief Class defining a single column dimension.
*/
//------------------------------------------------------------------------------
class ColumnDimension {
public:

  /*! \brief ColumnDimension class constructor
      \param nameIn Name of the dimension.
      \param descriptionIn Description of the dimension.
      \param sizeIn Size of the dimension.
  */
  ColumnDimension(const std::string nameIn, const std::string descriptionIn, const int sizeIn);

  /*! \brief default destructor */
  ~ColumnDimension() = default;

  /*! \brief Return the size of the dimension.
      \return the size of the dimension.
  */
  int size(void) const;

  /*! \brief Return the name of the dimension.
      \return the name of the dimension.
  */
  std::string name(void) const;

  /*! \brief Return the description of the dimension.
      \return the description of the dimension.
  */
  std::string description(void) const;

  /*! \brief less than comparison operator */
  bool operator<(const ColumnDimension &right) const;

  /*! \brief Write the ColumnDimension object to a stream.
      \param os Stream to write to.
      \param columnDimension ColumnDimension object to write to stream.

      The time interval is written to a stream in the format:
      "ColumnDimension: {name = name"
      "                  size = size"
      "                  desc = description}"
  */
  friend std::ostream & operator<<(std::ostream & os, const ColumnDimension & columnDimension);

private:

  /*! Dimension name */
  std::string dimensionName;

  /*! Dimension description */
  std::string dimensionDescription;

  /*! Dimension size */
  int dimensionSize;

};

//------------------------------------------------------------------------------
/*! \class ColumnDimensions
  \brief Class containing all the column dimensions
*/
//------------------------------------------------------------------------------
class ColumnDimensions {
public:

  /*! \brief ColumnDimensions class constructor
      \param configs Pointer to the configs object
      \param logIn Pointer to the log object
      \param useColumn If true we are using the full column physics
  */
  ColumnDimensions(DEMSI::Configs* configs, DEMSI::Log* logIn, const bool useColumn);

  /*! \brief default destructor */
  ~ColumnDimensions() = default;

  /*! \brief add a column dimension to the column dimensions
      \param name The name of the new dimension.
      \param description The description of the new dimension
      \param size The size of the new dimension
  */
  void add(const std::string name, const std::string description, const int size);

  /*! \brief Get a pointer to a column dimension based on its name
      \param name name of the column dimension to return
      \return Pointer to the column dimension
  */
  DEMSI::ColumnDimension* get(const std::string name) const;

  /*! \brief Get a pointer to a column dimension based on its
    position in the column dimensions vector
      \param iDimension index of dimension in the dimensions vector
      \return Pointer to the column dimension
  */
  DEMSI::ColumnDimension* get(const int iDimension) const;

  /*! \brief get the size of a dimension based on its name
      \param name of the column dimension
      \return size of the dimension
  */
  int size(const std::string name) const;

  /*! \brief get the number of column dimensions
      \return the number of column dimensions
  */
  int number(void) const;

  /*! \brief Write all the ColumnDimension objects in the
    columnDimensions object to a stream.
      \param os Stream to write to.
      \param columnDimensions ColumnDimensions object to write contents of to stream.
  */
  friend std::ostream & operator<<(std::ostream & os, const ColumnDimensions & columnDimensions);

private:

  /*! vector of column dimensions */
  std::vector<DEMSI::ColumnDimension*> dimensions;

  /*! Pointer to the log object */
  DEMSI::Log* log;

};

//------------------------------------------------------------------------------
/*! \class ColumnVariable
  \brief Class defining a column variable

  This class defines a column variable that stores column and BGC data for each particle.
  The class is templated to allow multiple types of data to be stored and supports \<double\>
  and \<int\> types. The underlying data is stored as a 1D vector.
*/
//------------------------------------------------------------------------------
template <typename T>
class ColumnVariable {
public:

  /*! \brief ColumnVariable class constructor
      \param particlesIn Pointer to the particles object
      \param columnDimensions Pointer to the columnDimensions object
      \param name Name of the variable
      \param dimensionNames vector of column dimension names to create the
      variable with. nParticles is assumed to be the outermost dimension
      and is not part of this vector.
      \param useVariable Whether to initialize the variable or not.
  */
  ColumnVariable(DEMSI::Particles* particlesIn, DEMSI::ColumnDimensions* columnDimensions, const std::string name, const std::vector<std::string> dimensionNames, const bool useVariable = true);

  /*! \brief default destructor */
  ~ColumnVariable() = default;

  /*! \brief Init the array vector
   */
  void init_array(void);

  /*! \brief Return true if variable is active
      \return True if variable is active
  */
  bool active(void) const;

  /*! \brief Return true if array is allocated
      \return True if array is allocated
  */
  bool allocated(void) const;

  /*! \brief Set the variable to a single value
      \param value Value to set variable to
  */
  void set(const T value);

  /*! \brief Set a variable for one particle to a constant value
      \param iParticle Particle index to set
      \param value Value to set variable to
  */
  void set(const int iParticle, const T value);

  /*! \brief Set the variable to the value of another variable
      \param value Variable whose array will be copied
  */
  void set(const DEMSI::ColumnVariable<T> value);

  /*! \brief Multiplication assignment operator.
      \param factor Factor to multiply variable by.
  */
  void operator*=(const double factor);

  /*! \brief Return the name of the variable
      \return the name of the variable
  */
  std::string name(void) const;

  /*! \brief Return a pointer to the start of the data for a particle.
      \param iParticle index of particle to return data for.
  */
  T* get(const int iParticle);

  /*! \brief Access a single element of the variable with dimensions (nParticles)
      \param iParticle Particle index to access
  */
  T& operator()(const int iParticle);

  /*! \brief Access a single element of the variable with dimensions (nParticles, n1)
      \param iParticle Particle index to access
      \param j index for the second dimension
  */
  T& operator()(const int iParticle, const int j);

  /*! \brief Access a single element of the variable with dimensions (nParticles, n1, n2)
      \param iParticle Particle index to access
      \param j index for the second dimension
      \param k index for the third dimension
  */
  T& operator()(const int iParticle, const int j, const int k);

  /*! \brief Get the dimension pointers of the variable
      \return Vector of ColumnDimension object pointers for this variable
  */
  std::vector<DEMSI::ColumnDimension*> get_dimension_pointers(void) const;

  /*! \brief Get the dimension sizes for the variable
      \return Vector of dimension sizes
  */
  std::vector<int> get_dimension_sizes(void) const;

  /*! \brief get the names of the dimensions for this variable
      \return vector of the names of the dimensions for this variable
  */
  std::vector<std::string> get_dimension_names(void) const;

  /*! \brief Determine if two variable have the same dimensions
      \param other Other variable to compare against
      \return True if variable have same dimensions
  */
  bool same_dimensions(const ColumnVariable<T>* other) const;

  /*! \brief get the variable array size per particle
      \return the variable array size per particle
  */
  int size_per_particle(void) const;

  /*! \brief Set a variable array pointer to the pointer of another array
      \param other The other variable to set the array pointer to
  */
  void set_array_pointer(const ColumnVariable<T>* other);

  /*! \brief Resize the underlying array
      \param newSize The desired new size of the underlying array

      This routine resizes the underlying array, and leaves any data in the array invalid.
  */
  void resize(const int newSize);

  /*! \brief return the number of particles for this variable
      \return The number of particles for this variable
  */
  int num_particles(void) const;

  /*! \brief Perform a global sum of the variable
      \param partition Pointer to the partition object
      \param log Pointer to the log object
  */
  double global_sum(DEMSI::Partition *partition, DEMSI::Log* log);

  /*! \brief Determine if column variable has any NaN values
      \return If column variable has any NaN values
  */
  bool has_nan(void) const;

  /*! \brief Return the minimum value of the variable
      \return The minimum value of the variable
  */
  T min(void) const;

  /*! \brief Return the maximum value of the variable
      \return The maximum value of the variable
  */
  T max(void) const;

  /*! \brief Return the location of the minimum value of the variable
      \return The location of the minimum value of the variable
  */
  int minloc(void) const;

  /*! \brief Return the location of the maximum value of the variable
      \return The location of the maximum value of the variable
  */
  int maxloc(void) const;

  /*! \brief Write variable object metadata to a stream.
      \param os Stream to write to.
      \param columnVariable variable object to write metadata contents of to stream.
  */
  template <typename T2>
  friend std::ostream & operator<<(std::ostream & os, const ColumnVariable<T2> & columnVariable);

private:

  /*! This variable is in use */
  bool useVariable;

  /*! This variable has its array allocated */
  bool arrayAllocated;

  /*! Pointer to the particles object */
  DEMSI::Particles* particles;

  /*! Name of the variable */
  std::string variableName;

  /*! vector of column dimensions for this variable */
  std::vector<DEMSI::ColumnDimension*> variableDimensions;

  /* vector of column dimension sizes */
  std::vector<int> dimSizes;

  /*! Length of array for each particle. */
  int sizePerParticle;

  /*! Data for the variable */
  std::vector<T>* array;

}; // ColumnVariable class

/*! Column field types: double or int */
enum FIELDTYPE {DOUBLE_FIELD, INT_FIELD};

//------------------------------------------------------------------------------
/*! \class ColumnVariables
  \brief Class providing access to all column variables
*/
//------------------------------------------------------------------------------
class ColumnVariables {
public:

  /*! \brief ColumnVariable class constructor
      \param logIn Pointer to log object
      \param partitionIn Pointer to partition object
      \param particlesIn Pointer to particles object
  */
  ColumnVariables(DEMSI::Log* logIn, DEMSI::Partition* partitionIn, DEMSI::Particles* particlesIn);

  /*! \brief default destructor */
  ~ColumnVariables() = default;

  /*! \brief add a pointer to a double column variable to the list of column variables
      \param columnVariable column variable pointer to add to list
  */
  void add(DEMSI::ColumnVariable<double>* columnVariable);

  /*! \brief add a pointer to a double column variable to the list of column variables
      \param columnVariable column variable pointer to add to list
  */
  void add(DEMSI::ColumnVariable<int>* columnVariable);

  /*! \brief Set the inactive variable pointer arrays to valid data */
  void set_arrays_inactive_variables(void);

  /*! \brief Check if a variable exists.
      \param fieldName Name of the field to check if exists
      \return If true named variable exists.
  */
  bool exists(const std::string fieldName) const;

  /*! \brief Return the type of the variable.
      \param fieldName name of field to get type of.
      \return Field type of the named field.
  */
  DEMSI::FIELDTYPE type(const std::string fieldName) const;

  /*! \brief Get the dimension names for a variable of given name
      \param fieldName The field name to get the dimension names of
      \return vector of dimension names for the named variable
  */
  std::vector<std::string> dimension_names(const std::string fieldName) const;

  /*! \brief Get the dimension sizes for a variable of given name
      \param fieldName The field name to get the dimension sizes of
      \return vector of dimension sizes for the named variable
  */
  std::vector<int> dimension_sizes(const std::string fieldName) const;

  /*! \brief Get a double field based on a name
      \param fieldName The field name to get
      \return the double field
  */
  DEMSI::ColumnVariable<double>* double_field(const std::string fieldName);

  /*! \brief Get an int field based on a name
      \param fieldName The field name to get
      \return the int field
  */
  DEMSI::ColumnVariable<int>* int_field(const std::string fieldName);

  /*! \brief Write variable object metadata for all variables to a stream.
      \param os Stream to write to.
      \param columnVariables variable objects to write metadata contents of to stream.
  */
  friend std::ostream & operator<<(std::ostream & os, const ColumnVariables & columnVariables);

  // processor transfer
  /*! \brief Transfer of column variables between processors

    This class transfers DEMSI column particle data between processors after
    LAMMPS has transfered its particle data during the model dynamics. After LAMMPS
    has communicated the particles owned by DEMSI for column physics are out of
    date and some particles will need to be sent to other processors and some got
    from other processors. DEMSI keeps track of the global ids of owned particles to
    facilitate this. When transfering column particle data, DEMSI first determines
    which globalIDs of particles it has gained, which it has lost and which it has
    kept since the last transference. The gained and lost lists are then sent
    around the processors in a circular manner with MPI commands allowing each
    processor to determine if it needs or has any of the particles on the lists
    being sent round. This way processors determine which processors have the
    particles they need to receive/send.
  */
  void processor_transfer(void);

  /*! \brief Initialize variables needed for transfer of column variables between processors */
  void init_processor_transfer(void);

  /*! \brief Initialize transfer column data between processors after remapping.
      \param nParticles Number of new particles
  */
  void init_processor_transfer_remap(const int nParticles);

  /*! \brief Determine if any column variables have NaN values */
  void has_nan(void) const;

  /*! list of pointers to double column variables */
  std::vector<DEMSI::ColumnVariable<double>*> doubleVariables;

  /*! list of pointers to int column variables */
  std::vector<DEMSI::ColumnVariable<int>*> intVariables;

private:

  /*! Pointer to log object */
  DEMSI::Log* log;

  /*! Pointer to partition object */
  DEMSI::Partition* partition;

  /*! Pointer to particles object */
  DEMSI::Particles* particles;

  // processor transfer
  /*! Set of globalIDs owned by this processor before LAMMPS particle transfer. */
  std::set <int> previousOwnedParticles;

  /*! Map from particle globalIDs to indices for before LAMMPS update. */
  std::map<int,int> globalIDsToIndicesPrevious;

};

} // namespace DEMSI

#endif /* COLUMN_VARIABLES_H_ */

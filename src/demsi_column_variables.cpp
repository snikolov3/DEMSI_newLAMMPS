#include "demsi_column_variables.h"

#include <string>
#include <vector>
#include <set>
#include "demsi_configs.h"
#include "demsi_particles.h"
#include <cmath>
#include <limits>

namespace DEMSI {

//------------------------------------------------------------------------------
// ColumnDimension class
//------------------------------------------------------------------------------

// ColumnDimension class constructor
ColumnDimension::ColumnDimension(const std::string nameIn, const std::string descriptionIn, const int sizeIn) {
  dimensionName = nameIn;
  dimensionDescription = descriptionIn;
  dimensionSize = sizeIn;
}

// Return the size of the dimension.
int ColumnDimension::size(void) const {
  return dimensionSize;
}

// Return the name of the dimension.
std::string ColumnDimension::name(void) const {
  return dimensionName;
}

// Return the description of the dimension.
std::string ColumnDimension::description(void) const {
  return dimensionDescription;
}

// less than comparison operator
bool ColumnDimension::operator<(const ColumnDimension &right) const {
  return dimensionName < right.dimensionName;
}

// over load << operator
std::ostream & operator<<(std::ostream & os, const ColumnDimension & columnDimension)
{
  os << "Dimension: {name = " << columnDimension.name() << std::endl;
  os << "            size = " << columnDimension.size() << std::endl;
  os << "            desc = " << columnDimension.description() << "}" << std::endl;
  return os;
}

//------------------------------------------------------------------------------
// ColumnDimensions class
//------------------------------------------------------------------------------

// ColumnDimensions class constructor
ColumnDimensions::ColumnDimensions(DEMSI::Configs* configs, DEMSI::Log* logIn, const bool useColumn) {

  log = logIn;

  if (useColumn) {

    int nCategories;
    configs->get({"ConfigGroup:column_dimensions","Config:nCategories"}, nCategories);
    add("nCategories", "The number of ice thickness categories.", nCategories);

    int nIceLayers;
    configs->get({"ConfigGroup:column_dimensions","Config:nIceLayers"},  nIceLayers);
    add("nIceLayers", "The number of ice layers in the vertical direction.", nIceLayers);

    int nSnowLayers;
    configs->get({"ConfigGroup:column_dimensions","Config:nSnowLayers"}, nSnowLayers);
    add("nSnowLayers", "The number of snow layers in the vertical direction.", nSnowLayers);

    int nAerosols;
    configs->get({"ConfigGroup:column_dimensions","Config:nAerosols"}, nAerosols);
    add("nAerosols", "The number of snow layers in the vertical direction.", nAerosols);

    int nIceLayersP1 = nIceLayers + 1;
    add("nIceLayersP1", "The number of ice layers in the vertical direction plus one.", nIceLayersP1);

    int nCategoriesP1 = nCategories + 1;
    add("nCategoriesP1", "The number of ice thickness categories plus one.", nCategoriesP1);

  } else {

    add("nCategories", "The number of ice thickness categories.", 1);

  } // useColumn

}

// add a column dimension to the column dimensions
void ColumnDimensions::add(const std::string name, const std::string description, const int size) {

  // check not already added
  for (int iDimension = 0 ; iDimension < dimensions.size() ; iDimension++) {
    log->check(dimensions[iDimension]->name() != name, "Dimension already defined: ", name);
  }

  DEMSI::ColumnDimension* dimension = new DEMSI::ColumnDimension(name, description, size);
  dimensions.push_back(dimension);

}

// Get a pointer to a column dimension based on its name
DEMSI::ColumnDimension* ColumnDimensions::get(const std::string name) const {

  for (int iDimension = 0 ; iDimension < dimensions.size() ; iDimension++) {
    if (dimensions[iDimension]->name() == name) {
      return dimensions[iDimension];
    }
  }
  log->abort("ColumnDimensions::get: Could not find dimension: ", name);

}

// Get a pointer to a column dimension based on its position in the column dimensions vector
DEMSI::ColumnDimension* ColumnDimensions::get(const int iDimension) const {

  return dimensions[iDimension];

}

// get the size of a dimension based on its name
int ColumnDimensions::size(const std::string name) const {

  for (int iDimension = 0 ; iDimension < dimensions.size() ; iDimension++) {
    if (dimensions[iDimension]->name() == name) {
      return dimensions[iDimension]->size();
    }
  }
  log->abort("ColumnDimensions::size: Could not find dimension: ", name);

}

// get the number of column dimensions
int ColumnDimensions::number(void) const {
  return dimensions.size();
}

// over load << operator
std::ostream & operator<<(std::ostream & os, const ColumnDimensions & columnDimensions)
{
  os << "ColumnDimensions: {" << std::endl;
  for (int iDimension = 0 ; iDimension < columnDimensions.number() ; iDimension++) {
    DEMSI::ColumnDimension* columnDimension = columnDimensions.get(iDimension);
    os << (*columnDimension);
  }
  os << "}" << std::endl;
  return os;
}

//------------------------------------------------------------------------------
// ColumnVariable class
//------------------------------------------------------------------------------

// ColumnVariable class constructor
template <typename T>
ColumnVariable<T>::ColumnVariable(DEMSI::Particles* particlesIn, DEMSI::ColumnDimensions* columnDimensions, const std::string variableNameIn, const std::vector<std::string> dimensionNames, const bool useVariableIn) {

  useVariable = useVariableIn;

  arrayAllocated = false;

  variableName = variableNameIn;

  particles = particlesIn;

  // get dimensions for this variable
  for (int iDimension = 0 ; iDimension < dimensionNames.size() ; iDimension++) {
    variableDimensions.push_back(columnDimensions->get(dimensionNames[iDimension]));
  }

  // get the size per particle and dim sizes
  sizePerParticle = 1;
  for (int iDimension = 0 ; iDimension < variableDimensions.size() ; iDimension++) {
    sizePerParticle *= variableDimensions[iDimension]->size();
    dimSizes.push_back(variableDimensions[iDimension]->size());
  }

  // initialize data
  array = new std::vector<T>();
  if (useVariable) init_array();

}

// init the array vector
template <typename T>
void ColumnVariable<T>::init_array(void) {
  for (int i = 0 ; i < *(particles->nParticles)*sizePerParticle ; i++) {
    array->push_back((T) 0);
  }
  arrayAllocated = true;
}

// return if variable is active
template <typename T>
bool ColumnVariable<T>::active(void) const {
  return useVariable;
}

// return if variable is allocated
template <typename T>
bool ColumnVariable<T>::allocated(void) const {
  return arrayAllocated;
}

// set a variable to a constant value
template <typename T>
void ColumnVariable<T>::set(const T value) {
  for (int i = 0 ; i < *(particles->nParticles)*sizePerParticle ; i++) {
    (*array)[i] = (T) value;
  }
}

// set a variable for one particle to a constant value
template <typename T>
void ColumnVariable<T>::set(const int iParticle, const T value) {
  for (int i = 0 ; i < sizePerParticle ; i++) {
    (*array)[iParticle * sizePerParticle + i] = (T) value;
  }
}

// set a variable to another variable
template <typename T>
void ColumnVariable<T>::set(const DEMSI::ColumnVariable<T> value) {
  for (int i = 0 ; i < *(particles->nParticles)*sizePerParticle ; i++) {
    (*(this->array))[i] = (*(value.array))[i];
  }
}

// addition assignment
template <typename T>
void ColumnVariable<T>::operator*=(const double factor) {
  for (int i = 0 ; i < *(particles->nParticles)*sizePerParticle ; i++) {
    (*array)[i] = (*array)[i] * factor;
  }
}

// return the variable name
template <typename T>
std::string ColumnVariable<T>::name(void) const {
  return variableName;
}


// Return a pointer to the start of the data for a particle.
template <typename T>
T* ColumnVariable<T>::get(const int iParticle) {

  return &(*array)[iParticle * sizePerParticle];

}

// Access a single element of the variable with dimensions (nParticles)
template <typename T>
T& ColumnVariable<T>::operator()(const int iParticle) {

  return (*array)[iParticle];

}

// Access a single element of the variable with dimensions (nParticles, n1)
template <typename T>
T& ColumnVariable<T>::operator()(const int iParticle, const int j) {

  return (*array)[iParticle * sizePerParticle + j];

}

// Access a single element of the variable with dimensions (nParticles, n1, n2)
template <typename T>
T& ColumnVariable<T>::operator()(const int iParticle, const int j, const int k) {

  return (*array)[iParticle * sizePerParticle + j * dimSizes[1] + k];

}

// Get the dimension pointers of the variable
template <typename T>
std::vector<DEMSI::ColumnDimension*> ColumnVariable<T>::get_dimension_pointers(void) const {

  std::vector<DEMSI::ColumnDimension*> dimensionPointers;

  for (int iDimension = 0 ; iDimension < variableDimensions.size() ; iDimension++) {
    dimensionPointers.push_back(variableDimensions[iDimension]);
  }

  return dimensionPointers;

}

template <typename T>
std::vector<int> ColumnVariable<T>::get_dimension_sizes(void) const {
  return dimSizes;
}

// get the names of the dimensions for this variable
template <typename T>
std::vector<std::string> ColumnVariable<T>::get_dimension_names(void) const {

  std::vector<std::string> dimensionNames;

  for (int iDimension = 0 ; iDimension < variableDimensions.size() ; iDimension++) {
    dimensionNames.push_back(variableDimensions[iDimension]->name());
  }

  return dimensionNames;

}

// see if the dimensions of the variable are the same as another variable
template <typename T>
bool ColumnVariable<T>::same_dimensions(const ColumnVariable<T>* other) const {

  std::vector<DEMSI::ColumnDimension*> thisPointers  = this->get_dimension_pointers();
  std::vector<DEMSI::ColumnDimension*> otherPointers = other->get_dimension_pointers();

  if (thisPointers.size() != otherPointers.size()) {
    return false;
  } else {

    for (int iDimension = 0 ; iDimension < thisPointers.size() ; iDimension++) {
      if (thisPointers[iDimension] != otherPointers[iDimension]) {
	return false;
      }
    } // iDimension

  }

  return true;

}

// get the variable array size per particle
template <typename T>
int ColumnVariable<T>::size_per_particle(void) const {

  return sizePerParticle;

}

// set the array pointer to another variable
template <typename T>
void ColumnVariable<T>::set_array_pointer(const ColumnVariable<T>* other) {
  this->array = other->array;
}

// resize the underlying array
template <typename T>
void ColumnVariable<T>::resize(const int newSize) {
  array->resize(newSize * sizePerParticle);
}

// return the number of particles for this variable
template <typename T>
int ColumnVariable<T>::num_particles(void) const {
  return *(particles->nParticles);
}

// Perform a global sum of the variable
template <typename T>
double ColumnVariable<T>::global_sum(DEMSI::Partition *partition, DEMSI::Log* log) {

  double localSum = 0;
  for (int i = 0 ; i < array->size(); i++) {
    localSum += (double) (*array)[i];
  }

  int err;
  char mpiErrBuffer[MPI_MAX_ERROR_STRING];
  int mpiErrLen;

  double globalSum;
  err = MPI_Reduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem reducing column variable during global sum: ", (std::string) mpiErrBuffer, ", for variable: ", variableName);

  err = MPI_Bcast(&globalSum, 1, MPI_DOUBLE, partition->master_proc(), partition->comm());
  MPI_Error_string(err, mpiErrBuffer, &mpiErrLen);
  log->check(err == MPI_SUCCESS, "Problem broadcasting column variable during global sum: ", (std::string) mpiErrBuffer, ", for variable: ", variableName);

  return globalSum;

}

// return the number of particles for this variable
template <typename T>
bool ColumnVariable<T>::has_nan(void) const {
  bool hasNan = false;
  for (int i = 0 ; i < array->size() ; i++) {
    if (std::isnan((*array)[i])) hasNan = true;
  }
  return hasNan;
}

// Return the minimum value of the variable
template <typename T>
T ColumnVariable<T>::min(void) const {
  T minValue = std::numeric_limits<T>::max();
  for (int i = 0 ; i < array->size() ; i++) {
    minValue = std::min(minValue,(*array)[i]);
  } // i
  return minValue;
}

// Return the maximum value of the variable
template <typename T>
T ColumnVariable<T>::max(void) const {
  T maxValue = std::numeric_limits<T>::min();
  for (int i = 0 ; i < array->size() ; i++) {
    maxValue = std::max(maxValue,(*array)[i]);
  } // i
  return maxValue;
}

// Return the location of the minimum value of the variable
template <typename T>
int ColumnVariable<T>::minloc(void) const {
  T minValue = std::numeric_limits<T>::max();
  int iLoc = -1;
  for (int i = 0 ; i < array->size() ; i++) {
    if ((*array)[i] < minValue) {
      iLoc = i;
      minValue = (*array)[i];
    }
  } // i
  return iLoc;
}

// Return the location of the maximum value of the variable
template <typename T>
int ColumnVariable<T>::maxloc(void) const {
  T maxValue = std::numeric_limits<T>::min();
  int iLoc = -1;
  for (int i = 0 ; i < array->size() ; i++) {
    if ((*array)[i] > maxValue) {
      iLoc = i;
      maxValue = (*array)[i];
    }
  } // i
  return iLoc;
}

// over load << operator
template <typename T>
std::ostream & operator<<(std::ostream & os, const ColumnVariable<T> & columnVariable)
{
  os << "ColumnVariable: {name = " << columnVariable.name() << std::endl;
  os << "                 useVariable       = " << columnVariable.active() << std::endl;
  os << "                 arrayAllocated    = " << columnVariable.allocated() << std::endl;
  os << "                 size per particle = " << columnVariable.size_per_particle() << std::endl;
  std::vector<DEMSI::ColumnDimension*> dimensions = columnVariable.get_dimension_pointers();
  os << "                 nDims = " << dimensions.size() << std::endl;
  os << "dimensions: {" << std::endl;
  for (int iDim = 0 ; iDim < dimensions.size() ; iDim++) {
    os << (*dimensions[iDim]);
  }
  os << "values: {" << std::endl;
  for (int i = 0 ; i < columnVariable.num_particles()*columnVariable.size_per_particle() ; i++) {
    os << i << " " << (*columnVariable.array)[i] << std::endl;
  } // i
  os << "}}" << std::endl;
  return os;
}

// define what types possible for ColumnVariable
template class ColumnVariable<double>;
template class ColumnVariable<int>;
template std::ostream &operator<<(std::ostream &os, ColumnVariable<double> const&);
template std::ostream &operator<<(std::ostream &os, ColumnVariable<int> const&);

//------------------------------------------------------------------------------
// ColumnVariables class
//------------------------------------------------------------------------------

// ColumnVariable class constructor
ColumnVariables::ColumnVariables(DEMSI::Log* logIn, DEMSI::Partition* partitionIn, DEMSI::Particles* particlesIn) {
  log = logIn;
  partition = partitionIn;
  particles = particlesIn;

  init_processor_transfer();
};

// add a pointer to a double column variable to the list of column variables
void ColumnVariables::add(DEMSI::ColumnVariable<double>* columnVariable) {
  doubleVariables.push_back(columnVariable);
}

// add a pointer to a double column variable to the list of column variables
void ColumnVariables::add(DEMSI::ColumnVariable<int>* columnVariable) {
  intVariables.push_back(columnVariable);
}

void ColumnVariables::set_arrays_inactive_variables(void) {

  // iterate over fields and find inactive variables
  for (int iFieldInactive = 0 ; iFieldInactive < doubleVariables.size() ; iFieldInactive++) {
    if (not doubleVariables[iFieldInactive]->active()) {

      bool foundActiveVariable = false;

      // iterate over allocated fields
      for (int iFieldAllocated = 0 ; iFieldAllocated < doubleVariables.size() ; iFieldAllocated++) {
	if (doubleVariables[iFieldAllocated]->allocated() and
	    doubleVariables[iFieldInactive]->same_dimensions(doubleVariables[iFieldAllocated])) {

	  doubleVariables[iFieldInactive]->set_array_pointer(doubleVariables[iFieldAllocated]);
	  foundActiveVariable = true;
	  break;

	} // active variables
      } // iFieldActive

      // if no active variable has the same dimensions init the array
      if (not foundActiveVariable) {
	doubleVariables[iFieldInactive]->init_array();
      }

    } // inactive variables
  } // iFieldInactive

}

// Check if a variable exists
bool ColumnVariables::exists(const std::string fieldName) const {

  for (int iField = 0 ; iField < doubleVariables.size() ; iField++) {
    if (doubleVariables[iField]->name() == fieldName) {
      return true;
    }
  }

  for (int iField = 0 ; iField < intVariables.size() ; iField++) {
    if (intVariables[iField]->name() == fieldName) {
      return true;
    }
  }

  return false;

}

// return the type of the variable
DEMSI::FIELDTYPE ColumnVariables::type(const std::string fieldName) const {

  for (int iField = 0 ; iField < doubleVariables.size() ; iField++) {
    if (doubleVariables[iField]->name() == fieldName) {
      return DEMSI::FIELDTYPE::DOUBLE_FIELD;
    }
  }

  for (int iField = 0 ; iField < intVariables.size() ; iField++) {
    if (intVariables[iField]->name() == fieldName) {
      return DEMSI::FIELDTYPE::INT_FIELD;
    }
  }

  log->abort("ColumnVariables::type: Column variable ", fieldName, "not found");

}

// Get the dimension names for a variable of given name
std::vector<std::string> ColumnVariables::dimension_names(const std::string fieldName) const {

  for (int iField = 0 ; iField < doubleVariables.size() ; iField++) {
    if (doubleVariables[iField]->name() == fieldName) {
      std::vector<std::string> dimNames = doubleVariables[iField]->get_dimension_names();
      return dimNames;
    }
  }

  for (int iField = 0 ; iField < intVariables.size() ; iField++) {
    if (intVariables[iField]->name() == fieldName) {
      std::vector<std::string> dimNames = doubleVariables[iField]->get_dimension_names();
      return dimNames;
    }
  }

  log->abort("ColumnVariables::dimension_names: Column variable ", fieldName, "not found");

}

// Get the dimension names for a variable of given name
std::vector<int> ColumnVariables::dimension_sizes(const std::string fieldName) const {

  for (int iField = 0 ; iField < doubleVariables.size() ; iField++) {
    if (doubleVariables[iField]->name() == fieldName) {
      std::vector<int> dimSizes = doubleVariables[iField]->get_dimension_sizes();
      return dimSizes;
    }
  }

  for (int iField = 0 ; iField < intVariables.size() ; iField++) {
    if (intVariables[iField]->name() == fieldName) {
      std::vector<int> dimSizes = doubleVariables[iField]->get_dimension_sizes();
      return dimSizes;
    }
  }

  log->abort("ColumnVariables::dimension_sizes: Column variable ", fieldName, "not found");

}


// Get a double field based on a name
DEMSI::ColumnVariable<double>* ColumnVariables::double_field(const std::string fieldName) {

  for (int iField = 0 ; iField < doubleVariables.size() ; iField++) {
    if (doubleVariables[iField]->name() == fieldName) {
      return doubleVariables[iField];
    }
  }

  log->abort("ColumnVariables::double_field: Column variable ", fieldName, "not found");

}

// Get an int field based on a name
DEMSI::ColumnVariable<int>* ColumnVariables::int_field(const std::string fieldName) {

  for (int iField = 0 ; iField < intVariables.size() ; iField++) {
    if (intVariables[iField]->name() == fieldName) {
      return intVariables[iField];
    }
  }

  log->abort("ColumnVariables::int_field: Column variable ", fieldName, "not found");

}

// Check if a variable exists
void ColumnVariables::has_nan(void) const {

  for (int iField = 0 ; iField < doubleVariables.size() ; iField++) {
    if (doubleVariables[iField]->has_nan()) {
      (*log)() << "Column variable: " << doubleVariables[iField]->name() << " has NaNs" << std::endl;
    }
  }

  for (int iField = 0 ; iField < intVariables.size() ; iField++) {
    if (intVariables[iField]->has_nan()) {
      (*log)() << "Column variable: " << intVariables[iField]->name() << " has NaNs" << std::endl;
    }
  }

}

// Write variable object metadata for all variables to a stream.
std::ostream & operator<<(std::ostream & os, const ColumnVariables & columnVariables) {

  os << "ColumnVariables: {" << std::endl;
  os << "Double: {" << std::endl;
  for (int iVar = 0 ; iVar < columnVariables.doubleVariables.size() ; iVar++) {
    os << (*columnVariables.doubleVariables[iVar]);
  }
  os << "} Int: {" << std::endl;
  for (int iVar = 0 ; iVar < columnVariables.intVariables.size() ; iVar++) {
    os << (*columnVariables.intVariables[iVar]);
  }
  os << "}}" << std::endl;
}

} // namespace DEMSI

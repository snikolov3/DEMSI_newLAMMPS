#ifndef DEMSI_CONFIGS_H_
#define DEMSI_CONFIGS_H_

#include <iostream>
#include <string>
#include <list>
#include <algorithm>
#include <exception>
#include <vector>

#include "demsi_tinyxml2.h"
#include "demsi_logging.h"

namespace DEMSI {

enum CONFIGTYPE {INTEGER, DOUBLE, STRING, BOOLEAN};

/*! \class Configs
    \brief Class defining model configuration parameters.

    The Configs class allows an xml configuration file to be read in and
    parameters stored in that file parsed and stored. These parameters are
    then available for use within the code.
*/
class Configs {

public:

  /*! \brief Configs constructor
      \param configFilename Input xml config filename.
      \param logIn Pointer to the Log object to use for logging.
  */
  Configs(const std::string configFilename, DEMSI::Log* logIn);

  /*! Default destructor.*/
  ~Configs() = default;

  /*! \brief Determine if a config is in the configs file.
      \param elementLocation Location of the configuration in the xml config
      file. Each entry in the list moves the position one level down in the
      xml file. Each entry is a string in the form "Tag:Name" where Tag is
      the tag name of the xml entry and Name is the "name" attribute value.
      \return If the config exists in the configs file.

      This should be used VERY sparingly, as defaults in code are discouraged.
  */
  bool exists(std::list<std::string> elementLocation);

  /*! \brief Get a std::string configuration attribute.
      \param elementLocation Location of the configuration in the xml config
      file. Each entry in the list moves the position one level down in the
      xml file. Each entry is a string in the form "Tag:Name" where Tag is
      the tag name of the xml entry and Name is the "name" attribute value.
      \param attributeName The desired attribute name from the desired
      xml file element.
      \param value The std::string configuration attribute value to get.
  */
  void get_attribute(std::list<std::string> elementLocation, const std::string attributeName, std::string &value);

  /*! \brief Get a integer configuration attribute.
      \param elementLocation Location of the configuration in the xml config
      file. Each entry in the list moves the position one level down in the
      xml file. Each entry is a string in the form "Tag:Name" where Tag is
      the tag name of the xml entry and Name is the "name" attribute value.
      \param attributeName The desired attribute name from the desired
      xml file element.
      \param value The integer configuration attribute value to get.
  */
  void get_attribute(std::list<std::string> elementLocation, const std::string attributeName, int &value);

  /*! \brief Get a double configuration attribute.
      \param elementLocation Location of the configuration in the xml config
      file. Each entry in the list moves the position one level down in the
      xml file. Each entry is a string in the form "Tag:Name" where Tag is
      the tag name of the xml entry and Name is the "name" attribute value.
      \param attributeName The desired attribute name from the desired
      xml file element.
      \param value The double configuration attribute value to get.
  */
  void get_attribute(std::list<std::string> elementLocation, const std::string attributeName, double &value);

  /*! \brief Get a boolean configuration attribute.
      \param elementLocation Location of the configuration in the xml config
      file. Each entry in the list moves the position one level down in the
      xml file. Each entry is a string in the form "Tag:Name" where Tag is
      the tag name of the xml entry and Name is the "name" attribute value.
      \param attributeName The desired attribute name from the desired
      xml file element.
      \param value The boolean configuration attribute value to get.
  */
  void get_attribute(std::list<std::string> elementLocation, const std::string attributeName, bool &value);

  /*! \brief Return an array of std::string from the xml config file. By default
      this is the contents of the elements with tag type arrayType, but can be an
      attribute of those elements if attributeName is specified.
      \param elementLocation The location of the array in the xml config file.
      \param arrayType The element type to extract from the array.
      \param attributeName The attribute to extract from the array elements. If
      this is not specified the element contents are returned.
      \return The desired array of strings.
  */
  std::vector<std::string> get_array(std::list<std::string> elementLocation, std::string arrayType, std::string attributeName="");

  /*! \brief Return an array of double from the xml config file. By default
      this is the contents of the elements with tag type arrayType, but can be an
      attribute of those elements if attributeName is specified.
      \param elementLocation The location of the array in the xml config file.
      \param arrayType The element type to extract from the array.
      \param attributeName The attribute to extract from the array elements. If
      this is not specified the element contents are returned.
      \return The desired array of doubles.
  */
  std::vector<double> get_double_array(std::list<std::string> elementLocation, std::string arrayType, std::string attributeName="");

  /*! \brief Get the type of a configuration parameter.
      \param elementLocation Location of the configuration in the xml config file.
      \return The type of the configuration parameter.
  */
  DEMSI::CONFIGTYPE get_config_type(std::list<std::string> elementLocation);

  /*! \brief Get an integer configuration value.
      \param elementLocation Location of the configuration in the xml config file.
      \param config The integer configuration parameter value to get.
  */
  void get(std::list<std::string> elementLocation, int &config);

  /*! \brief Get a double configuration value.
      \param elementLocation Location of the configuration in the xml config file.
      \param config The double configuration parameter value to get.
  */
  void get(std::list<std::string> elementLocation, double &config);

  /*! \brief Get a string configuration value.
      \param elementLocation Location of the configuration in the xml config file.
      \param config The string configuration parameter value to get.
  */
  void get(std::list<std::string> elementLocation, std::string &config);

  /*! \brief Get a boolean configuration value.
      \param elementLocation Location of the configuration in the xml config file.
      \param config The boolean configuration parameter value to get.
  */
  void get(std::list<std::string> elementLocation, bool &config);

  /*! \brief Get an integer configuration value.
      \param elementLocation Location of the configuration in the xml config file.
      \param config The integer configuration parameter value to get.
      \param hasDefault Return default (0) if config not present in configs file.
  */
  void get(std::list<std::string> elementLocation, int &config, bool hasDefault);

  /*! \brief Get a double configuration value.
      \param elementLocation Location of the configuration in the xml config file.
      \param config The double configuration parameter value to get.
      \param hasDefault Return default (0.0) if config not present in configs file.
  */
  void get(std::list<std::string> elementLocation, double &config, bool hasDefault);

  /*! \brief Get a string configuration value.
      \param elementLocation Location of the configuration in the xml config file.
      \param config The string configuration parameter value to get.
      \param hasDefault Return default ("") if config not present in configs file.
  */
  void get(std::list<std::string> elementLocation, std::string &config, bool hasDefault);

  /*! \brief Get a boolean configuration value.
      \param elementLocation Location of the configuration in the xml config file.
      \param config The boolean configuration parameter value to get.
      \param hasDefault Return default (false) if config not present in configs file.
  */
  void get(std::list<std::string> elementLocation, bool &config, bool hasDefault);

  /*! \brief Return a boolean config if exists, false otherwise
      \param elementLocation Location of the configuration in the xml config
      file.
      \return Config boolean if exists or false otherwise.
  */
  bool is_true(std::list<std::string> elementLocation);

private:

  // xml configs file object
  tinyxml2::XMLDocument xmlDoc;

  // Pointer to logging object
  DEMSI::Log* log;

  /*! \brief Get a pointer to an element of the xml config file.
      \param elementLocation List of strings defining the location of the
      desired element in the xml config file. Each string defines a
      parent element of the desired element in the format "Tag:Name"
      where Tag is the tag name of the xml entry and Name is the "name"
      attribute value of the xml file.
      \return Pointer to the desired element.
  */
  tinyxml2::XMLElement* get_element(std::list<std::string> elementLocation);

};

} // namespace DEMSI

#endif

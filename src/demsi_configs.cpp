#include "demsi_logging.h"
#include "demsi_configs.h"

#include <vector>

namespace DEMSI {

// recursive function to check configs
void check_children_options(tinyxml2::XMLElement* parent, DEMSI::Log* log) {

  tinyxml2::XMLElement* child = parent->FirstChildElement();
  while (child != nullptr) {

    if ((std::string) child->Value() == "Config") {

      std::string configName  = (std::string) child->Attribute("name");
      std::string configType  = (std::string) child->Attribute("type");
      std::string configValue = (std::string) child->Attribute("value");

      // check types
      if (configType == "Double") {
	try {
	  double value = std::stod(configValue);
	} catch(std::exception& e) {
	  log->abort("Non-double input '", configValue, "' for double config '", configName, "'.");
	}
      } else if (configType == "Integer") {
	try {
	  int value = std::stoi(configValue);
	} catch(std::exception& e) {
	  log->abort("Non-integer input '", configValue, "' for integer config '", configName, "'.");
	}
      } else if (configType == "Boolean") {
	log->check(configValue == "true" or configValue == "false", "Unknown boolean value '", configValue, "' for config '", configName, "'.");
      } else if (configType == "String") {

      } else {
	log->abort("Unknown config type '", configType, "' for config '", configName, "'.");
      } // configType

      // check options
      tinyxml2::XMLElement* options = child->FirstChildElement();
      while (options != nullptr) {

	if ((std::string) options->Value() == "Array" and
	    (std::string) options->Attribute("name") == "Options") {

	  std::vector<std::string> optionValues;

	  tinyxml2::XMLElement* option = options->FirstChildElement();
	  while (option != nullptr) {

	    if ((std::string) option->Value() == "Option") {

	      optionValues.push_back(option->GetText());

	    } // check type option

	    option = option->NextSiblingElement();
	  } // option

	  log->check(std::find(optionValues.begin(), optionValues.end(), configValue) != optionValues.end(),
		     "Option '", configValue,"' for config '", configName, "' is not an allowed option.");

	} // check is options array

	options = options->NextSiblingElement();
      } // options

    } // child is config

    // recursive
    check_children_options(child, log);

    child = child->NextSiblingElement();
  } // level3

} // check_level

// class constructor
Configs::Configs(const std::string configFilename, DEMSI::Log* logIn) {

  log = logIn;

  // load the configs file
  tinyxml2::XMLError xmlResult = xmlDoc.LoadFile(configFilename.c_str());
  log->check(xmlResult == tinyxml2::XML_SUCCESS, "Problem opening config file: ", configFilename);

  // check options
  tinyxml2::XMLNode* root = xmlDoc.FirstChild();
  log->check(root != nullptr, "Failed to parse doc first child");

  // check specified config option against allowed options
  check_children_options(root->ToElement(), log);

};

// Turn an element location list into a single string for diagnostics
std::string element_location_str(std::list<std::string> elementLocation) {

  std::string elementLocationStr = "{";

  std::list<std::string>::const_iterator itr;
  for (itr = elementLocation.begin(); itr != elementLocation.end(); ++itr) {
    elementLocationStr = elementLocationStr + *itr + ", ";
  }

  elementLocationStr = elementLocationStr + "}";
  return elementLocationStr;

}

// Get a pointer to an element of the xml config file.
tinyxml2::XMLElement* Configs::get_element(std::list<std::string> elementLocationIn) {

  std::list<std::string> elementLocation = elementLocationIn;

  std::string elementType;
  std::string elementName;
  {
    std::istringstream strStream(elementLocation.front());
    std::getline(strStream, elementType, ':');
    std::getline(strStream, elementName);
    elementLocation.pop_front();
  }

  tinyxml2::XMLNode* root = xmlDoc.FirstChild();
  log->check(root != nullptr, "Failed to parse doc first child: ", element_location_str(elementLocationIn));

  tinyxml2::XMLElement* element = root->FirstChildElement();
  log->check(element != nullptr, "Failed to parse root first child element: ", element_location_str(elementLocationIn));

  while (element) {

    if ((std::string) element->Value() == elementType and element->Attribute("name") == elementName) {
      if (elementLocation.size() > 0) {
	element = element->FirstChildElement();
	std::istringstream strStream(elementLocation.front());
	std::getline(strStream, elementType, ':');
	std::getline(strStream, elementName);
	elementLocation.pop_front();
      } else {
	return element;
      }
    } else {
      element = element->NextSiblingElement();
    }

  }

  return nullptr;

}

// check if a config exists
bool Configs::exists(std::list<std::string> elementLocation) {

  tinyxml2::XMLElement* element = get_element(elementLocation);
  if (element == nullptr) {
    return false;
  } else {
    return true;
  }

}

// Get a std::string configuration attribute.
void Configs::get_attribute(std::list<std::string> elementLocation, const std::string attributeName, std::string &value) {

  tinyxml2::XMLElement* element = get_element(elementLocation);
  log->check(element != nullptr, "Cant find config element: ", element_location_str(elementLocation));
  value = element->Attribute(attributeName.c_str());

}

// Get a integer configuration attribute.
void Configs::get_attribute(std::list<std::string> elementLocation, const std::string attributeName, int &value) {

  tinyxml2::XMLElement* element = get_element(elementLocation);
  log->check(element != nullptr, "Cant find config element: ", element_location_str(elementLocation));
  tinyxml2::XMLError xmlResult = element->QueryIntAttribute(attributeName.c_str(), &value);
  log->check(xmlResult == tinyxml2::XML_SUCCESS, "Failed to parse integer value for config: ", element_location_str(elementLocation));
  return;

}

// Get a double configuration attribute.
void Configs::get_attribute(std::list<std::string> elementLocation, const std::string attributeName, double &value) {

  tinyxml2::XMLElement* element = get_element(elementLocation);
  log->check(element != nullptr, "Cant find config element: ", element_location_str(elementLocation));
  tinyxml2::XMLError xmlResult = element->QueryDoubleAttribute(attributeName.c_str(), &value);
  log->check(xmlResult == tinyxml2::XML_SUCCESS, "Failed to parse double value for config: ", element_location_str(elementLocation));
  return;

}

// Get a boolean configuration attribute.
void Configs::get_attribute(std::list<std::string> elementLocation, const std::string attributeName, bool &value) {

  tinyxml2::XMLElement* element = get_element(elementLocation);
  log->check(element != nullptr, "Cant find config element: ", element_location_str(elementLocation));
  tinyxml2::XMLError xmlResult = element->QueryBoolAttribute(attributeName.c_str(), &value);
  log->check(xmlResult == tinyxml2::XML_SUCCESS, "Failed to parse bool value for config: ", element_location_str(elementLocation));
  return;

}

// Return an array of std::string from the xml config file.
std::vector<std::string> Configs::get_array(std::list<std::string> elementLocation, const std::string arrayType, const std::string attributeName) {

  std::vector<std::string> array;

  tinyxml2::XMLElement* element = get_element(elementLocation);
  log->check(element != nullptr, "Cant find config element: ", element_location_str(elementLocation));
  element = element->FirstChildElement();
  while (element) {
    if ((std::string) element->Value() == arrayType) {
      if (attributeName == "") {
	array.push_back(element->GetText());
      } else {
	array.push_back(element->Attribute(attributeName.c_str()));
      }
    }
    element = element->NextSiblingElement();
  }

  return array;

}

// Return an array of double from the xml config file.
std::vector<double> Configs::get_double_array(std::list<std::string> elementLocation, const std::string arrayType, const std::string attributeName) {

  std::vector<double> array;

  tinyxml2::XMLElement* element = get_element(elementLocation);
  log->check(element != nullptr, "Cant find config element: ", element_location_str(elementLocation));
  element = element->FirstChildElement();
  while (element) {
    if ((std::string) element->Value() == arrayType) {
      if (attributeName == "") {
	array.push_back(std::atof(element->GetText()));
      } else {
	array.push_back(std::atof(element->Attribute(attributeName.c_str())));
      }
    }
    element = element->NextSiblingElement();
  }

  return array;

}

// get the config type
DEMSI::CONFIGTYPE Configs::get_config_type(std::list<std::string> elementLocation) {

  std::string configType;
  get_attribute(elementLocation, (std::string) "type", configType);

  if (configType == "Integer") {
    return DEMSI::CONFIGTYPE::INTEGER;
  } else if (configType == "Double") {
    return DEMSI::CONFIGTYPE::DOUBLE;
  } else if (configType == "String") {
    return DEMSI::CONFIGTYPE::STRING;
  } else if (configType == "Boolean") {
    return DEMSI::CONFIGTYPE::BOOLEAN;
  } else {
    log->abort("Unknown config type for config: ", element_location_str(elementLocation));
  }

}

// Get an integer configuration value.
void Configs::get(std::list<std::string> elementLocation, int &config) {
  log->check(get_config_type(elementLocation) == DEMSI::CONFIGTYPE::INTEGER, "Wrong config type request for config: ", element_location_str(elementLocation));
  get_attribute(elementLocation, "value", config);
}

// Get an double configuration value.
void Configs::get(std::list<std::string> elementLocation, double &config) {
  log->check(get_config_type(elementLocation) == DEMSI::CONFIGTYPE::DOUBLE, "Wrong config type request for config: ", element_location_str(elementLocation));
  get_attribute(elementLocation, "value", config);
}

// Get an string configuration value.
void Configs::get(std::list<std::string> elementLocation, std::string &config) {
  log->check(get_config_type(elementLocation) == DEMSI::CONFIGTYPE::STRING, "Wrong config type request for config: ", element_location_str(elementLocation));
  get_attribute(elementLocation, "value", config);
}

// Get an boolean configuration value.
void Configs::get(std::list<std::string> elementLocation, bool &config) {
  log->check(get_config_type(elementLocation) == DEMSI::CONFIGTYPE::BOOLEAN, "Wrong config type request for config: ", element_location_str(elementLocation));
  get_attribute(elementLocation, "value", config);
}

// Get an integer configuration value.
void Configs::get(std::list<std::string> elementLocation, int &config, bool hasDefault) {
  config = 0;
  if (not hasDefault or (hasDefault and exists(elementLocation))) {
    log->check(get_config_type(elementLocation) == DEMSI::CONFIGTYPE::INTEGER, "Wrong config type request for config: ", element_location_str(elementLocation));
    get_attribute(elementLocation, "value", config);
  }
}

// Get an double configuration value.
void Configs::get(std::list<std::string> elementLocation, double &config, bool hasDefault) {
  config = 0.0;
  if (not hasDefault or (hasDefault and exists(elementLocation))) {
    log->check(get_config_type(elementLocation) == DEMSI::CONFIGTYPE::DOUBLE, "Wrong config type request for config: ", element_location_str(elementLocation));
    get_attribute(elementLocation, "value", config);
  }
}

// Get an string configuration value.
void Configs::get(std::list<std::string> elementLocation, std::string &config, bool hasDefault) {
  config = "";
  if (not hasDefault or (hasDefault and exists(elementLocation))) {
    log->check(get_config_type(elementLocation) == DEMSI::CONFIGTYPE::STRING, "Wrong config type request for config: ", element_location_str(elementLocation));
    get_attribute(elementLocation, "value", config);
  }
}

// Get an boolean configuration value.
void Configs::get(std::list<std::string> elementLocation, bool &config, bool hasDefault) {
  config = false;
  if (not hasDefault or (hasDefault and exists(elementLocation))) {
    log->check(get_config_type(elementLocation) == DEMSI::CONFIGTYPE::BOOLEAN, "Wrong config type request for config: ", element_location_str(elementLocation));
    get_attribute(elementLocation, "value", config);
  }
}

// Return a boolean config if exists, false otherwise
bool Configs::is_true(std::list<std::string> elementLocation) {
  if (this->exists(elementLocation)) {
    bool isTrue;
    this->get(elementLocation, isTrue);
    return isTrue;
  } else {
    return false;
  }
}

} // namespace DEMSI

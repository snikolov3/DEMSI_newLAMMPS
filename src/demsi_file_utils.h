#ifndef DEMSI_FILE_UTILS_H_
#define DEMSI_FILE_UTILS_H_

#include <string>
#include "demsi_time.h"

/*! \file   demsi_file_utils.h
    \brief  Header file for collection of file utility functions.
*/

namespace DEMSI {

/*! /brief Replace a token in a string with another string.
    /param str String to modify.
    /param from Toekn to replace.
    /param to String to replace token with.
*/
void replace_time_token(std::string& str, const std::string& from, const std::string& to);

/*! /brief Expand the time tokens in a filename template.
    /param filenameTemplate Filename template to expand.
    /param fileTime Time to expand the filename template with.
*/
std::string expand_filename_template(const std::string filenameTemplate, const Time fileTime);

} // namespace DEMSI

#endif /* FILE_UTILS_H_ */

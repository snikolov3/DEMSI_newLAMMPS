#include "demsi_file_utils.h"

namespace DEMSI {

// Replace a single token in a filename template with a substitute string.
void replace_time_token(std::string& str, const std::string& from, const std::string& to) {
  size_t startPos = str.find(from);
  if (startPos != std::string::npos) {
    str.replace(startPos, from.length(), to);
  }
}

// Expand a filename template with a specific read time.
std::string expand_filename_template(const std::string filenameTemplate, const Time fileTime) {

  std::string filename = filenameTemplate;

  // get components of time
  int years, months, days, hours, minutes, seconds, milliseconds;
  fileTime.get_components(years, months, days, hours, minutes, seconds, milliseconds);
  months = months + 1;
  days = days + 1;

  // replace tokens in template
  // year
  char strYear[5];
  sprintf(strYear, "%4.4d", years);
  replace_time_token(filename, (std::string) "$Y", (std::string) strYear);

  // months
  char strMonth[3];
  sprintf(strMonth, "%2.2d", months);
  replace_time_token(filename, (std::string) "$M", (std::string) strMonth);

  // days
  char strDay[3];
  sprintf(strDay, "%2.2d", days);
  replace_time_token(filename, (std::string) "$D", (std::string) strDay);

  // hours
  char strHour[3];
  sprintf(strHour, "%2.2d", hours);
  replace_time_token(filename, (std::string) "$h", (std::string) strHour);

  // minutes
  char strMinute[3];
  sprintf(strMinute, "%2.2d", minutes);
  replace_time_token(filename, (std::string) "$m", (std::string) strMinute);

  // seconds
  char strSecond[3];
  sprintf(strSecond, "%2.2d", seconds);
  replace_time_token(filename, (std::string) "$s", (std::string) strSecond);

  // milliseconds
  char strMillisecond[4];
  sprintf(strMillisecond, "%3.3d", milliseconds);
  replace_time_token(filename, (std::string) "$X", (std::string) strMillisecond);

  return filename;

}

} // namespace DEMSI

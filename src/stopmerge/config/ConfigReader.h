// Copyright 2016
// University of Freiburg - Chair of Algorithms and Data Structures
// Author: Patrick Brosi

#ifndef STOPMERGE_CONFIG_CONFIGREADER_H_
#define STOPMERGE_CONFIG_CONFIGREADER_H_

#include "stopmerge/config/StopMergeConfig.h"

namespace stopmerge {
namespace config {

class ConfigReader {
 public:
  ConfigReader();
  void read(StopMergeConfig* cfg, int argc, char** argv) const;
  void help(const char* bin) const;
};

}  // namespace config
}  // namespace stopmerge

#endif  // STOPMERGE_CONFIG_CONFIGREADER_H_

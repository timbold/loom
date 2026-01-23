// Copyright 2016
// University of Freiburg - Chair of Algorithms and Data Structures
// Author: Patrick Brosi

#ifndef STOPMERGE_CONFIG_STOPMERGECONFIG_H_
#define STOPMERGE_CONFIG_STOPMERGECONFIG_H_

#include <string>

namespace stopmerge {
namespace config {

struct StopMergeConfig {
  std::string mergeStops = "never";               // never|auto|always
  std::string parallelCorridors = "off";          // off|auto|always

  double mergeStopsSnapDistM = 40.0;
  double mergeStopsRadiusM = 40.0;
  double mergeStopsChainageM = 40.0;
  double mergeStopsNodeZoneM = 15.0;
  double mergeStopsCorridorAngleDeg = 25.0;
  size_t mergeStopsMaxClusterSize = 6;

  double parallelPairMaxDistM = 25.0;
  double parallelPairMaxAngleDeg = 20.0;
  double parallelPairMinOverlapRatio = 0.35;
  bool parallelPairRequireLines = true;
  double parallelPairMinGroupLengthM = 150.0;
  double parallelPairMaxMedianDistM = 30.0;

  double hubRadiusM = 80.0;
  size_t mergeStopsMaxLocalDensity = 8;
  bool hubAutoRequireBaseName = true;
  std::string mergeStopsHubFallback = "auto";  // off|auto|always

  std::string mergeStopsDebugCsv;
  std::string mergeStopsDebugRejectsCsv;
  std::string parallelPairDebugCsv;
  std::string stopMergeMapJson;
  std::string mergeStopsOverrides;
};

}  // namespace config
}  // namespace stopmerge

#endif  // STOPMERGE_CONFIG_STOPMERGECONFIG_H_

// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2GRAPH_CONFIG_GTFS2GEOCONFIG_H_
#define GTFS2GRAPH_CONFIG_GTFS2GEOCONFIG_H_

#include <cstddef>
#include <set>
#include <string>
#include "ad/cppgtfs/gtfs/flat/Route.h"

namespace gtfs2graph {
namespace config {

enum class CloseCircularRoutesMode { Never, Auto, Always };

struct Config {
  std::string inputFeedPath;

  double pruneThreshold = 0.0;

  CloseCircularRoutesMode closeCircularRoutesMode =
      CloseCircularRoutesMode::Never;
  double circularMaxCloseDistanceM = 5000.0;
  double circularDefaultSpeedKmh = 20.0;
  size_t circularMinStops = 4;
  size_t circularMinUniqueStops = 3;
  double circularProximityRatio = 0.15;

  std::set<ad::cppgtfs::gtfs::flat::Route::TYPE> useMots;
};

}  // namespace config
}  // namespace gtfs2graph

#endif  // GTFS2GRAPH_CONFIG_GTFS2TOPOCONFIG_H_

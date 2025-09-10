// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_RENDERGRAPH_LANDMARK_H_
#define SHARED_RENDERGRAPH_LANDMARK_H_

#include "util/geo/Geo.h"
#include <string>

namespace shared {
namespace rendergraph {

struct Landmark {
  std::string iconPath;
  std::string label;
  std::string color = "#474747";
  util::geo::DPoint coord;
  double size = 200;
  std::string cssClass = "landmark";
};

} // namespace rendergraph
} // namespace shared

#endif // SHARED_RENDERGRAPH_LANDMARK_H_

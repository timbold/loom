// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_
#define TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_

#include "shared/rendergraph/Landmark.h"
#include <string>
#include <vector>

namespace transitmapper {
namespace config {

using shared::rendergraph::Landmark;

struct Config {
  double lineWidth = 20;
  double lineSpacing = 10;

  double lineLabelSize = 40;
  // Maximum allowed bend angle in radians for line label candidates.
  double lineLabelBendAngle = 0.3490658503988659; // ~20 degrees
  // Maximum allowed ratio of polyline length to straight-line distance.
  double lineLabelLengthRatio = 1.1;
  double stationLabelSize = 60;
  double stationLineOverlapPenalty = 15;
  // Maximum font size for station labels in SVG output; -1 for no limit.
  double fontSvgMax = 11;
  // Gap between consecutive route label boxes.
  double routeLabelBoxGap = 10;
  // Gap between the terminus station label and the first route label box.
  double routeLabelTerminusGap = 80;

  std::string renderMethod = "svg";

  std::string mvtPath = ".";

  bool writeStats = false;

  double outputResolution = 0.1;
  double inputSmoothing = 1;
  double innerGeometryPrecision = 3;

  double outputPadding = -1;
  double paddingTop = -1;
  double paddingRight = -1;
  double paddingBottom = -1;
  double paddingLeft = -1;

  double ratio = -1;
  double tlRatio = -1;

  double outlineWidth = 1;
  std::string outlineColor;

  bool renderStations = true;
  bool renderNodeFronts = false;
  bool renderNodeCircles = false;
  bool renderEdges = true;
  bool renderLabels = false;
  bool renderRouteLabels = false;
  bool highlightTerminals = false;
  bool dontLabelDeg2 = false;
  bool fromDot = false;

  bool randomColors = false;

  bool renderNodeConnections = true;
  bool tightStations = false;

  std::vector<size_t> mvtZooms;

  bool renderDirMarkers = false;
  bool renderMarkersTail = false;
  bool renderBiDirMarker = false;
  size_t crowdedLineThresh = 3;
  double sharpTurnAngle = 0.7853981633974483; // 45 degrees in radians
  std::string worldFilePath;

  std::vector<Landmark> landmarks;

  bool renderMe = false;
  Landmark meLandmark;
};

} // namespace config
} // namespace transitmapper

#endif // TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_

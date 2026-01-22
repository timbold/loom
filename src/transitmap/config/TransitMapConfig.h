// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_
#define TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_

#include <string>
#include <vector>

namespace transitmapper {
namespace config {

struct Config {
  double lineWidth = 20;
  double lineSpacing = 10;

  double lineLabelSize = 40;
  double stationLabelSize = 60;

  std::string renderMethod = "svg";

  bool writeStats = false;

  double outputResolution = 0.1;
  double inputSmoothing = 1;
  double innerGeometryPrecision = 3;

  double outputPadding = -1;

  double outlineWidth = 1;
  std::string outlineColor;

  bool renderStations = true;
  bool renderNodeFronts = false;
  bool renderNodeCircles = false;
  bool renderEdges = true;
  bool renderLabels = false;
  bool dontLabelDeg2 = false;
  bool fromDot = false;

  bool randomColors = false;

  bool renderNodeConnections = true;
  bool tightStations = false;

  bool renderDirMarkers = false;
  std::string worldFilePath;

  std::string mbtilesPath;
  std::string paper = "A4L";
  int canvasWidth = 1000;
  std::string canvasUnit = "px";
  std::vector<int> zoomLevels;
  int maxTiles = 512;
  double oversample = 1.25;
  double backgroundPadPct = 0.03;
  double backgroundOpacity = 1.0;
};

}  // namespace config
}  // namespace transitmapper

#endif  // TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_

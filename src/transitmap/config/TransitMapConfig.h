// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_
#define TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_

#include "shared/rendergraph/Landmark.h"
#include "util/geo/Geo.h"
#include "util/log/Log.h"
#include <string>
#include <vector>

namespace transitmapper {
namespace config {

using shared::rendergraph::Landmark;

enum class TerminusLabelAnchor {
  StationLabel,
  StopFootprint,
  Node,
};

struct Config {
  double lineWidth = 20;
  double lineSpacing = 10;

  double lineLabelSize = 40;
  // Maximum allowed bend angle in radians for line label candidates.
  double lineLabelBendAngle = 0.3490658503988659; // ~20 degrees
  // Maximum allowed ratio of polyline length to straight-line distance.
  double lineLabelLengthRatio = 1.1;
  double stationLabelSize = 60;
  // Text size for the optional "YOU ARE HERE" label.
  double meLabelSize = 80;
  // Size of the star marker for --me.
  double meStarSize = 150;
  double stationLineOverlapPenalty = 15;
  // Count distinct transit lines when penalizing station-line overlaps.
  bool stationLineOverlapPerLine = false;
  // Penalize candidates when far label ends crowd nearby edges.
  double stationLabelFarCrowdRadius = 0;
  double stationLabelFarCrowdPenalty = 25;
  double sidePenaltyWeight = 2.5;
  // Penalty for placing station labels on the opposite side of a connecting
  // edge.
  double sameSidePenalty = 100;
  // Scale factor applied to sameSidePenalty in the crowding relief pass.
  double crowdingSameSideScale = 0.5;
  // Additional passes to reposition station labels after initial placement.
  int repositionLabel = 0;
  // Scale factor for the station crowding penalty.
  double clusterPenScale = 1.0;
  // Penalty (positive) or bonus (negative) for labels outside the map bounds.
  double outsidePenalty = -5.0;
  std::vector<double> orientationPenalties = {0, 3, 6, 4, 1, 5, 6, 2};
  // Maximum font size for station labels in SVG output; -1 for no limit.
  double fontSvgMax = 11;
  // Gap between consecutive route label boxes.
  double routeLabelBoxGap = 10;
  // Gap between the terminus station label and the first route label box.
  double routeLabelTerminusGap = 80;
  // Control the geometry used to anchor terminus route label stacks.
  TerminusLabelAnchor terminusLabelAnchor = TerminusLabelAnchor::StationLabel;
  // Maximum number of lateral shifts to try when avoiding collisions.
  int terminusLabelMaxLateralShift = 2;
  // Arrange route labels in multiple columns at termini when enabled.
  bool compactTerminusLabel = false;
  // Stack route labels above edges into multiple rows when enabled.
  bool compactRouteLabel = false;

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

  bool renderDirMarkers = false;
  size_t dirMarkerSpacing = 1;
  bool renderMarkersTail = false;
  bool tailIgnoreSharpAngle = false;
  bool renderBiDirMarker = false;
  size_t crowdedLineThresh = 3;
  double sharpTurnAngle = 0.7853981633974483; // 45 degrees in radians
  std::string bgMapPath;
  // Background map coordinates are interpreted as latitude/longitude by
  // default and converted to Web Mercator. Set when the coordinates are
  // already in Web Mercator projection.
  bool bgMapWebmerc = false;
  // Extend output bounds with background map geometry when enabled.
  bool extendWithBgMap = false;
  // Opacity for background map geometry.
  double bgMapOpacity = 1.0;
  std::string worldFilePath;

  // Ensure output covers at least a specific geographic bounding box.
  bool geoLock = false;
  util::geo::Box<double> geoLockBox;

  // Landmark coordinates are interpreted as latitude/longitude by default
  // and converted to Web Mercator. Set when coordinates are already in Web
  // Mercator projection.
  bool landmarksWebmerc = false;

  std::vector<Landmark> landmarks;

  // Render landmarks even when overlapping existing geometry.
  bool renderOverlappingLandmarks = true;
  // Radius (in steps) for searching free positions for landmark icons.
  int landmarkSearchRadius = 10;
  // Maximum iterations for force-directed landmark displacement.
  int displacementIterations = 100;
  // Cooling factor for landmark displacement step size.
  double displacementCooling = 0.9;

  std::string meStation;
  std::string meStationFill = "#f00";
  std::string meStationBorder;
  bool meStationWithBg = false;
  std::string meStationBgFill = "#f5f5f5";
  std::string meStationBgStroke = "#d0d0d0";
  std::string meStationTextColor = "#3a3a3a";

  bool renderMe = false;
  bool renderMeLabel = false;
  Landmark meLandmark;

  int logLevel = INFO;
};

} // namespace config
} // namespace transitmapper

#endif // TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_

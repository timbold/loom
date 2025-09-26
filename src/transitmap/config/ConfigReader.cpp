// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
#include <getopt.h>
#include <limits.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#ifdef _WIN32
#include <windows.h>
#endif

#include "shared/rendergraph/Landmark.h"
#include "transitmap/_config.h"
#include "transitmap/config/ConfigReader.h"
#include "transitmap/util/String.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/log/Log.h"

using shared::rendergraph::Landmark;
using transitmapper::config::Config;
using transitmapper::config::ConfigReader;
using transitmapper::config::TerminusLabelAnchor;

namespace {

// Return the directory part of a path. If no directory is present, '.' is
// returned. Both POSIX '/' and Windows '\\' separators are handled.
std::string dirName(const std::string &path) {
  size_t pos = path.find_last_of("/\\");
  if (pos == std::string::npos)
    return ".";
  return path.substr(0, pos);
}

// Join two paths using '/' as separator if necessary.
std::string joinPath(const std::string &base, const std::string &rel) {
  if (base.empty() || base == ".")
    return rel;
  char last = base.back();
  if (last == '/' || last == '\\')
    return base + rel;
  return base + "/" + rel;
}

// Return the canonical/absolute version of a path. If resolving fails,
// the original path is returned.
std::string canonicalPath(const std::string &path) {
#ifdef _WIN32
  char buf[_MAX_PATH];
  if (_fullpath(buf, path.c_str(), _MAX_PATH)) {
    return std::string(buf);
  }
#else
  char buf[PATH_MAX];
  if (realpath(path.c_str(), buf)) {
    return std::string(buf);
  }
#endif
  return path;
}

// Keep track of landmark files that have already been processed so the
// same file is not parsed twice when provided multiple times.
std::unordered_set<std::string> processedLandmarkFiles;

// Explicit option codes for getopt_long so we don't clash with ASCII values.
// GEO_LOCK used to be 58 which conflicted with ';' when interpreted as a
// character. 260 is safely outside the signed char range and unique in this
// file.
constexpr int OPT_GEO_LOCK = 260;
// Assign a unique code outside the ASCII range for long options without a
// short equivalent to avoid clashes with character options such as 'D'.
constexpr int OPT_REPOSITION_LABEL = 261;
constexpr int OPT_TERMINUS_LABEL_ANCHOR = 262;
constexpr int OPT_STATION_LABEL_FAR_CROWD_RADIUS = 263;
constexpr int OPT_STATION_LABEL_FAR_CROWD_PENALTY = 264;
constexpr int OPT_STATION_LINE_OVERLAP_PER_LINE = 265;
constexpr int OPT_TERMINUS_LABEL_MAX_SHIFT = 266;
constexpr int OPT_ME_WITH_BG = 267;
constexpr int OPT_ME_BG_FILL = 268;
constexpr int OPT_ME_BG_STROKE = 269;
constexpr int OPT_ME_LABEL_COLOR = 270;
constexpr int OPT_TERMINUS_HIGHLIGHT_FILL = 271;
constexpr int OPT_TERMINUS_HIGHLIGHT_STROKE = 272;
constexpr int OPT_STATION_LABEL_ANGLE_STEPS = 273;
constexpr int OPT_STATION_LABEL_ANGLE_STEP_DEG = 274;
constexpr int OPT_TERMINUS_ANGLE_PENALTY = 275;
constexpr int OPT_ME_STAR = 276;
constexpr int OPT_NO_SINGLE_ROUTE_LABELS = 278;
bool toBool(const std::string &v) {
  std::string s = util::toLower(v);
  return s == "1" || s == "true" || s == "yes" || s == "on";
}

TerminusLabelAnchor parseTerminusLabelAnchor(const std::string &arg) {
  std::string value = util::toLower(arg);
  if (value == "station-label" || value == "station") {
    return TerminusLabelAnchor::StationLabel;
  }
  if (value == "stop-footprint" || value == "footprint" || value == "stop") {
    return TerminusLabelAnchor::StopFootprint;
  }
  if (value == "node" || value == "point") {
    return TerminusLabelAnchor::Node;
  }

  LOG(WARN) << "Unknown terminus label anchor '" << arg
            << "', using station-label";
  return TerminusLabelAnchor::StationLabel;
}

void applyOption(Config *cfg, int c, const std::string &arg,
                 const std::string &baseDir = "") {
  switch (c) {
  case 55:
    cfg->logLevel = atoi(arg.c_str());
    break;
  case 2:
    cfg->lineWidth = atof(arg.c_str());
    break;
  case 3:
    cfg->lineSpacing = atof(arg.c_str());
    break;
  case 4:
    cfg->outlineWidth = atof(arg.c_str());
    break;
  case 5:
    cfg->lineLabelSize = atof(arg.c_str());
    break;
  case 35:
    cfg->lineLabelBendAngle = atof(arg.c_str());
    break;
  case 36:
    cfg->lineLabelLengthRatio = atof(arg.c_str());
    break;
  case 6:
    cfg->stationLabelSize = atof(arg.c_str());
    break;
  case OPT_STATION_LABEL_ANGLE_STEPS: {
    int steps = atoi(arg.c_str());
    cfg->stationLabelAngleSteps = steps > 0 ? static_cast<size_t>(steps) : 0;
    break;
  }
  case OPT_STATION_LABEL_ANGLE_STEP_DEG:
    cfg->stationLabelAngleStepDeg = atof(arg.c_str());
    break;
  case 40:
    cfg->meLabelSize = atof(arg.c_str());
    cfg->meLandmark.fontSize = cfg->meLabelSize;
    cfg->meLabelSizeExplicit = true;
    break;
  case 38:
    cfg->fontSvgMax = atof(arg.c_str());
    break;
  case 37:
    cfg->stationLineOverlapPenalty = atof(arg.c_str());
    break;
  case OPT_STATION_LINE_OVERLAP_PER_LINE:
    cfg->stationLineOverlapPerLine = arg.empty() ? true : toBool(arg);
    break;
  case OPT_STATION_LABEL_FAR_CROWD_RADIUS:
    cfg->stationLabelFarCrowdRadius = atof(arg.c_str());
    break;
  case OPT_STATION_LABEL_FAR_CROWD_PENALTY:
    cfg->stationLabelFarCrowdPenalty = atof(arg.c_str());
    break;
  case 61:
    cfg->sidePenaltyWeight = atof(arg.c_str());
    break;
  case 67:
    cfg->sameSidePenalty = atof(arg.c_str());
    break;
  case OPT_REPOSITION_LABEL:
    cfg->repositionLabel = atoi(arg.c_str());
    break;
  case 69:
    cfg->crowdingSameSideScale = atof(arg.c_str());
    break;
  case 62: {
    auto parts = util::split(arg, ',');
    if (parts.size() == 8) {
      cfg->orientationPenalties.clear();
      for (const auto &p : parts) {
        cfg->orientationPenalties.push_back(atof(p.c_str()));
      }
    }
    break;
  }
  case OPT_TERMINUS_ANGLE_PENALTY:
    cfg->terminusAnglePenalty = atof(arg.c_str());
    break;
  case 65:
    cfg->clusterPenScale = atof(arg.c_str());
    break;
  case 66:
    cfg->outsidePenalty = atof(arg.c_str());
    break;
  case 32:
    cfg->routeLabelBoxGap = atof(arg.c_str());
    break;
  case 34:
    cfg->routeLabelTerminusGap = atof(arg.c_str());
    break;
  case OPT_TERMINUS_LABEL_ANCHOR:
    cfg->terminusLabelAnchor = parseTerminusLabelAnchor(arg);
    break;
  case OPT_TERMINUS_LABEL_MAX_SHIFT:
    cfg->terminusLabelMaxLateralShift =
        std::max(0, atoi(arg.c_str()));
    break;
  case 33:
    cfg->highlightTerminals = arg.empty() ? true : toBool(arg);
    break;
  case OPT_TERMINUS_HIGHLIGHT_FILL:
    cfg->terminusHighlightFill = arg;
    break;
  case OPT_TERMINUS_HIGHLIGHT_STROKE:
    cfg->terminusHighlightStroke = arg;
    break;
  case 48:
    cfg->compactTerminusLabel = arg.empty() ? true : toBool(arg);
    break;
  case 49:
    cfg->compactRouteLabel = arg.empty() ? true : toBool(arg);
    break;
  case OPT_NO_SINGLE_ROUTE_LABELS:
    cfg->renderSingleRouteLabel = arg.empty() ? false : !toBool(arg);
    break;
  case 7:
    cfg->renderStations = arg.empty() ? false : !toBool(arg);
    break;
  case 'l':
    cfg->renderLabels = arg.empty() ? true : toBool(arg);
    break;
  case 'r':
    cfg->renderRouteLabels = arg.empty() ? true : toBool(arg);
    break;
  case 9:
    cfg->tightStations = arg.empty() ? true : toBool(arg);
    break;
  case 10:
    cfg->renderDirMarkers = arg.empty() ? true : toBool(arg);
    break;
  case 50:
    cfg->dirMarkerSpacing = atoi(arg.c_str());
    break;
  case 20:
    cfg->renderMarkersTail = arg.empty() ? true : toBool(arg);
    break;
  case 47:
    cfg->tailIgnoreSharpAngle = arg.empty() ? true : toBool(arg);
    break;
  case 11:
    cfg->renderNodeConnections = arg.empty() ? false : !toBool(arg);
    break;
  case 12:
    cfg->outputResolution = atof(arg.c_str());
    break;
  case 13:
    cfg->outputPadding = atof(arg.c_str());
    break;
  case 23:
    cfg->paddingTop = atof(arg.c_str());
    break;
  case 24:
    cfg->paddingRight = atof(arg.c_str());
    break;
  case 25:
    cfg->paddingBottom = atof(arg.c_str());
    break;
  case 26:
    cfg->paddingLeft = atof(arg.c_str());
    break;
  case 14:
    cfg->inputSmoothing = atof(arg.c_str());
    break;
  case 27:
    cfg->ratio = atof(arg.c_str());
    break;
  case 31:
    cfg->tlRatio = atof(arg.c_str());
    break;
  case 15:
    cfg->renderNodeFronts = arg.empty() ? true : toBool(arg);
    break;
  case 28:
    cfg->crowdedLineThresh = atoi(arg.c_str());
    break;
  case 29: {
    double ang = atof(arg.c_str());
    if (ang > M_PI) {
      ang = ang * M_PI / 180.0;
    }
    if (ang <= 0 || ang > M_PI) {
      std::cerr << "Error: sharp-turn-angle " << ang << " is out of range (0-π)"
                << std::endl;
      exit(1);
    }
    cfg->sharpTurnAngle = ang;
    break;
  }
  case 30:
    cfg->renderBiDirMarker = arg.empty() ? true : toBool(arg);
    break;
  case 52:
    cfg->bgMapPath = arg;
    break;
  case 53:
    cfg->bgMapWebmerc = arg.empty() ? true : toBool(arg);
    break;
  case 60:
    cfg->bgMapOpacity = atof(arg.c_str());
    break;
  case 57:
    cfg->extendWithBgMap = arg.empty() ? true : toBool(arg);
    break;
  case OPT_GEO_LOCK:
    cfg->geoLock = arg.empty() ? true : toBool(arg);
    if (cfg->geoLock && cfg->geoLockBox.getLowerLeft().getX() >
                            cfg->geoLockBox.getUpperRight().getX()) {
      util::geo::DPoint ll(106.7105, 47.8521);
      util::geo::DPoint ur(107.0209, 47.9504);
      ll = util::geo::latLngToWebMerc(ll);
      ur = util::geo::latLngToWebMerc(ur);
      cfg->geoLockBox = util::geo::Box<double>(ll, ur);
    }
    break;
  case 59: {
    auto parts = util::split(arg, ',');
    if (parts.size() == 4) {
      double south = atof(parts[0].c_str());
      double west = atof(parts[1].c_str());
      double north = atof(parts[2].c_str());
      double east = atof(parts[3].c_str());
      util::geo::DPoint ll(west, south);
      util::geo::DPoint ur(east, north);
      ll = util::geo::latLngToWebMerc(ll);
      ur = util::geo::latLngToWebMerc(ur);
      cfg->geoLockBox = util::geo::Box<double>(ll, ur);
      cfg->geoLock = true;
    }
    break;
  }
  case 54:
    cfg->landmarksWebmerc = arg.empty() ? true : toBool(arg);
    break;
  case 18:
    cfg->randomColors = arg.empty() ? true : toBool(arg);
    break;
  case 19:
    cfg->writeStats = arg.empty() ? true : toBool(arg);
    break;
  case 21: {
    auto parts = util::split(arg, ',');
    if (parts.size() >= 3) {
      Landmark l;
      double lat = atof(parts[1].c_str());
      double lon = atof(parts[2].c_str());
      util::geo::DPoint p(lon, lat);
      if (!cfg->landmarksWebmerc) {
        p = util::geo::latLngToWebMerc(p);
      }
      l.coord = p;
      std::string first = parts[0];
      bool isWord = first.rfind("word:", 0) == 0;
      if (isWord) {
        l.label = util::trim(first.substr(5));
        if (parts.size() >= 4)
          l.fontSize = atof(parts[3].c_str());
        if (parts.size() >= 5)
          l.color = parts[4];
        if (parts.size() >= 6)
          l.opacity = atof(parts[5].c_str());
      } else {
        l.iconPath = first;
        if (!baseDir.empty() && !l.iconPath.empty() &&
            l.iconPath.find(':') == std::string::npos && l.iconPath[0] != '/' &&
            l.iconPath[0] != '\\') {
          l.iconPath = joinPath(baseDir, l.iconPath);
        }
        if (parts.size() >= 4)
          l.size = atof(parts[3].c_str());
      }
      cfg->landmarks.push_back(l);
    }
    break;
  }
  case 22: {
    std::string canon = canonicalPath(arg);
    if (processedLandmarkFiles.find(canon) != processedLandmarkFiles.end()) {
      LOG(WARN) << "Landmarks file '" << arg
                << "' already processed, skipping.";
      break;
    }
    std::ifstream in(arg.c_str());
    if (in.good()) {
      processedLandmarkFiles.insert(canon);
      std::string l;
      std::string dir = dirName(arg);
      while (std::getline(in, l)) {
        // Landmark parsing (case 21) handles coordinate conversion.
        applyOption(cfg, 21, util::trim(l), dir);
      }
    }
    break;
  }
  case 51:
    cfg->renderOverlappingLandmarks = arg.empty() ? true : toBool(arg);
    break;
  case 56:
    cfg->landmarkSearchRadius = atoi(arg.c_str());
    break;
  case 63:
    cfg->displacementIterations = atoi(arg.c_str());
    break;
  case 64:
    cfg->displacementCooling = atof(arg.c_str());
    break;
  case 41:
    cfg->meStarSize = atof(arg.c_str());
    cfg->meStarSizeExplicit = true;
    break;
  case 42:
    cfg->renderMeLabel = arg.empty() ? true : toBool(arg);
    if (cfg->renderMeLabel) {
      cfg->meLandmark.label = "YOU ARE HERE";
      cfg->meLandmark.fontSize = cfg->meLabelSize;
    }
    break;
  case OPT_ME_STAR: {
    cfg->forceMeStar = arg.empty() ? true : toBool(arg);
    if (cfg->forceMeStar) {
      cfg->meLandmark.color = cfg->meStationFill;
    }
    break;
  }
  case 39: {
    auto parts = util::split(arg, ',');
    if (parts.size() == 2) {
      cfg->renderMe = true;
      double lat = atof(parts[0].c_str());
      double lon = atof(parts[1].c_str());
      util::geo::DPoint p(lon, lat);
      if (!cfg->landmarksWebmerc) {
        p = util::geo::latLngToWebMerc(p);
      }
      cfg->meLandmark.coord = p;
      cfg->meLandmark.size = cfg->meStarSize;
      cfg->meLandmark.label = "";
    }
    break;
  }
  case 43: {
    std::string trimmed = util::trim(arg);
    std::string sanitized = util::sanitizeStationLabel(trimmed);
    cfg->meStation = sanitized;
    cfg->meStationId = sanitized;
    cfg->meStationLabel = trimmed;
    break;
  }
  case 44:
    cfg->meStationFill = arg;
    cfg->meLandmark.color = cfg->meStationFill;
    break;
  case 45:
    cfg->meStationBorder = arg;
    break;
  case OPT_ME_WITH_BG:
    cfg->meStationWithBg = arg.empty() ? true : toBool(arg);
    cfg->highlightMeStationLabel = cfg->meStationWithBg;
    break;
  case OPT_ME_BG_FILL:
    cfg->meStationBgFill = arg;
    break;
  case OPT_ME_BG_STROKE:
    cfg->meStationBgStroke = arg;
    break;
  case OPT_ME_LABEL_COLOR:
    cfg->meStationTextColor = arg;
    break;
  case 'D':
    cfg->fromDot = arg.empty() ? true : toBool(arg);
    break;
  case 16:
    cfg->dontLabelDeg2 = arg.empty() ? true : toBool(arg);
    break;
  default:
    break;
  }
}

} // namespace

static const char *YEAR = &__DATE__[7];
static const char *COPY =
    "University of Freiburg - Chair of Algorithms and Data Structures";
static const char *AUTHORS = "Patrick Brosi <brosi@informatik.uni-freiburg.de>";

// _____________________________________________________________________________
ConfigReader::ConfigReader() {}

// _____________________________________________________________________________
void ConfigReader::help(const char *bin) const {
  std::cout
      << std::setfill(' ') << std::left << "transitmap (part of LOOM) "
      << VERSION_FULL << "\n(built " << __DATE__ << " " << __TIME__ << ")"
      << "\n\n(C) 2017-" << YEAR << " " << COPY << "\n"
      << "Authors: " << AUTHORS << "\n\n"
      << "Usage: " << bin << " < linegraph.json\n\n"
      << "Allowed options:\n\n"
      << "General:\n"
      << std::setw(37) << "  -v [ --version ]"
      << "print version\n"
      << std::setw(37) << "  -h [ --help ]"
      << "show this help message\n"
      << std::setw(37) << "  --config arg"
      << "read options from config file\n"
      << std::setw(37) << "  --log-level arg (=2)"
      << "log verbosity 0-4\n"
      << std::setw(37) << "  --line-width arg (=20)"
      << "width of a single transit line\n"
      << std::setw(37) << "  --line-spacing arg (=10)"
      << "spacing between transit lines\n"
      << std::setw(37) << "  --outline-width arg (=1)"
      << "width of line outlines\n"
      << std::setw(37) << "  --render-dir-markers"
      << "render line direction markers\n"
      << std::setw(37) << "  --render-markers-tail"
      << "add tail to direction markers\n"
      << std::setw(37) << "  --dir-marker-spacing arg (=1)"
      << "edges between forced direction markers\n"
      << std::setw(37) << "  --bi-dir-marker"
      << "render markers for bidirectional edges\n"
      << std::setw(37) << "  --crowded-line-thresh arg (=3)"
      << "lines on edge to trigger direction marker\n"
      << std::setw(37) << "  --sharp-turn-angle arg (=0.785398)"
      << "turn angle in radians (0-PI) to trigger direction marker; "
      << "values >PI are treated as degrees\n"
      << std::setw(37) << "  --tail-ignore-sharp-angle"
      << "ignore sharp angle check for marker tail\n"
      << std::setw(37) << "  -l [ --labels ]"
      << "render labels\n"
      << std::setw(37) << "  -r [ --route-labels ]"
      << "render route names at line termini\n"
      << std::setw(37) << "  --line-label-textsize arg (=40)"
      << "textsize for line labels\n"
      << std::setw(37) << "  --line-label-bend-angle arg (=0.349066)"
      << "max bend angle in radians for line label candidates\n"
      << std::setw(37) << "  --line-label-length-ratio arg (=1.1)"
      << "max length/straight distance ratio for line label candidates\n"
      << std::setw(37) << "  --station-label-textsize arg (=60)"
      << "textsize for station labels\n"
      << std::setw(37) << "  --station-label-angle-steps arg (=24)"
      << "number of station label orientations to sample\n"
      << std::setw(37)
      << "  --station-label-angle-step-deg arg (=15)"
      << "degrees between successive station label orientations\n"
      << std::setw(37) << "  --me-label-textsize arg (=80)"
      << "textsize for 'me' label\n"
      << std::setw(37) << "  --font-svg-max arg (=11)"
      << "max font size for station labels in SVG, -1 for no limit\n"
      << std::setw(37) << "  --station-line-overlap-penalty arg (=15)"
      << "penalty multiplier for station-line overlaps\n"
      << std::setw(37) << "  --station-line-overlap-per-line"
      << "count unique lines instead of edges for station overlaps\n"
      << std::setw(37)
      << "  --station-label-far-crowd-radius arg (=0)"
      << "radius from far label end to penalize nearby edges\n"
      << std::setw(37)
      << "  --station-label-far-crowd-penalty arg (=25)"
      << "penalty when far label end crowds nearby edges\n"
      << std::setw(37) << "  --side-penalty-weight arg (=2.5)"
      << "weight for station label side preference penalties\n"
      << std::setw(37) << "  --same-side-penalty arg (=100)"
      << "penalty for station labels on opposite sides\n"
      << std::setw(37) << "  --crowding-same-side-scale arg (=0.5)"
      << "scale factor for same-side penalty in crowding pass\n"
      << std::setw(37) << "  --reposition-label arg (=0)"
      << "re-run label placement to reduce crowding n times\n"
      << std::setw(37)
      << "  --orientation-penalties arg (=0,3,6,4,1,5,6,2)"
      << "penalties for 8 label orientations\n"
      << std::setw(37) << "  --terminus-angle-penalty arg (=3)"
      << "penalty for non-axis-aligned terminus station labels\n"
      << std::setw(37) << "  --cluster-pen-scale arg (=1)"
      << "scale factor for station crowding penalty\n"
      << std::setw(37) << "  --outside-penalty arg (=-5)"
      << "penalty or bonus for labels outside map bounds\n"
      << std::setw(37) << "  --route-label-gap arg (=20)"
      << "gap between route label boxes\n"
      << std::setw(37) << "  --route-label-terminus-gap arg (=100)"
      << "gap between terminus station label and route labels\n"
      << std::setw(37) << "  --terminus-label-anchor arg (=station-label)"
      << "anchor geometry for terminus route labels (station-label|stop-footprint|node)\n"
      << std::setw(37) << "  --terminus-label-max-shift arg (=2)"
      << "max lateral shifts tested when avoiding terminus collisions\n"
      << std::setw(37) << "  --highlight-terminal"
      << "highlight terminus stations\n"
      << std::setw(37) << "  --terminus-highlight-fill arg (=black)"
      << "fill color when highlighting terminus stations\n"
      << std::setw(37) << "  --terminus-highlight-stroke arg (=#BAB6B6)"
      << "stroke color when highlighting terminus stations\n"
      << std::setw(37) << "  --compact-terminal-label"
      << "arrange terminus route labels in multiple columns\n"
      << std::setw(37) << "  --compact-route-label"
      << "stack route labels above edges in multiple rows\n"
      << std::setw(37) << "  --no-single-route-labels"
      << "omit route labels when only one line terminates\n"
      << std::setw(37) << "  --no-deg2-labels"
      << "no labels for deg-2 stations\n"
      << "Misc:\n"
      << std::setw(37) << "  -D [ --from-dot ]"
      << "input is in dot format\n"
      << std::setw(37) << "  --padding arg (=-1)"
      << "padding, -1 for auto\n"
      << std::setw(37) << "  --padding-top arg (=-1)"
      << "top padding, -1 for auto\n"
      << std::setw(37) << "  --padding-right arg (=-1)"
      << "right padding, -1 for auto\n"
      << std::setw(37) << "  --padding-bottom arg (=-1)"
      << "bottom padding, -1 for auto\n"
      << std::setw(37) << "  --padding-left arg (=-1)"
      << "left padding, -1 for auto\n"
      << std::setw(37) << "  --smoothing arg (=1)"
      << "input line smoothing\n"
      << std::setw(37) << "  --ratio arg (=-1)"
      << "output width/height ratio\n"
      << std::setw(37) << "  --tl-ratio arg"
      << "top-left anchored width/height ratio\n"
      << std::setw(37) << "  --random-colors"
      << "fill missing colors with random colors\n"
      << std::setw(37) << "  --tight-stations"
      << "don't expand node fronts for stations\n"
      << std::setw(37) << "  --no-render-stations"
      << "don't render stations\n"
      << std::setw(37) << "  --no-render-node-connections"
      << "don't render inner node connections\n"
      << std::setw(37) << "  --render-node-fronts"
      << "render node fronts\n"
      << std::setw(37) << "  --bg-map arg"
      << "GeoJSON file with background geometry (lat/lon, WGS84)\n"
      << std::setw(37) << "  --bg-map-webmerc"
      << "background GeoJSON already in Web Mercator\n"
      << std::setw(37) << "  --bg-map-opacity arg (=1)"
      << "opacity for background map geometry\n"
      << std::setw(37) << "  --extend-with-bgmap"
      << "expand output bounds using background map geometry\n"
      << std::setw(37) << "  --geo-lock"
      << "ensure output covers default bbox\n"
      << std::setw(37) << "  --geo-lock-bbox arg"
      << "lock output to bbox south,west,north,east\n"
      << std::setw(37) << "  --landmark arg"
      << "add landmark word:text,lat,lon[,fontSize[,color[,opacity]]] or "
         "iconPath,lat,lon[,size] (size optional)\n"
      << std::setw(37) << "  --landmarks arg"
      << "read landmarks from file, one word:text,lat,lon[,fontSize[,color"
         "[,opacity]]] or iconPath,lat,lon[,size] per line\n"
      << std::setw(37) << "  --force-landmarks"
      << "render landmarks even if they overlap existing geometry (default)\n"
      << std::setw(37) << "  --landmark-search-radius arg (=10)"
      << "search radius for moving overlapping landmarks\n"
      << std::setw(37) << "  --displacement-iterations arg (=100)"
      << "max iterations for landmark displacement\n"
      << std::setw(37) << "  --displacement-cooling arg (=0.9)"
      << "cooling factor for landmark displacement\n"
      << std::setw(37) << "  --landmarks-webmerc"
      << "landmark and --me coordinates already in Web Mercator\n"
      << std::setw(37) << "  --me arg"
      << "mark current location lat,lon with star\n"
      << std::setw(37) << "  --me-size arg (=150)"
      << "size of 'me' star\n"
      << std::setw(37) << "  --me-label"
      << "add 'YOU ARE HERE' text\n"
      << std::setw(37) << "  --me-star[=<bool>]"
      << "force rendering of 'me' star (default false, fill from --me-station-fill)\n"
      << std::setw(37) << "  --me-station arg"
      << "mark current location by station label\n"
      << std::setw(37) << "  --me-with-bg"
      << "render station badge with background\n"
      << std::setw(37) << "  --me-bg-fill arg (=#f5f5f5)"
      << "badge fill color for --me background\n"
      << std::setw(37) << "  --me-bg-stroke arg (=#d0d0d0)"
      << "badge stroke color for --me background\n"
      << std::setw(37) << "  --me-label-color arg (=#3a3a3a)"
      << "text color for --me badge\n"
      << std::setw(37) << "  --me-station-fill arg (=#f00)"
      << "fill color for 'me' marker\n"
      << std::setw(37) << "  --me-station-border arg"
      << "border color for 'me' marker\n"
      << std::setw(37) << "  --print-stats"
      << "write stats to stdout\n";
}

// _____________________________________________________________________________
void ConfigReader::read(Config *cfg, int argc, char **argv) const {

  std::unordered_map<std::string, int> optMap = {
      {"line-width", 2},
      {"line-spacing", 3},
      {"outline-width", 4},
      {"log-level", 55},
      {"from-dot", 'D'},
      {"no-deg2-labels", 16},
      {"line-label-textsize", 5},
      {"line-label-bend-angle", 35},
      {"line-label-length-ratio", 36},
      {"station-label-textsize", 6},
      {"station-label-angle-steps", OPT_STATION_LABEL_ANGLE_STEPS},
      {"station-label-angle-step-deg", OPT_STATION_LABEL_ANGLE_STEP_DEG},
      {"me-label-textsize", 40},
      {"font-svg-max", 38},
      {"station-line-overlap-penalty", 37},
      {"station-line-overlap-per-line", OPT_STATION_LINE_OVERLAP_PER_LINE},
      {"station-label-far-crowd-radius", OPT_STATION_LABEL_FAR_CROWD_RADIUS},
      {"station-label-far-crowd-penalty", OPT_STATION_LABEL_FAR_CROWD_PENALTY},
      {"side-penalty-weight", 61},
      {"same-side-penalty", 67},
      {"reposition-label", OPT_REPOSITION_LABEL},
      {"crowding-same-side-scale", 69},
      {"orientation-penalties", 62},
      {"terminus-angle-penalty", OPT_TERMINUS_ANGLE_PENALTY},
      {"cluster-pen-scale", 65},
      {"outside-penalty", 66},
      {"route-label-gap", 32},
      {"route-label-terminus-gap", 34},
      {"terminus-label-anchor", OPT_TERMINUS_LABEL_ANCHOR},
      {"terminus-label-max-shift", OPT_TERMINUS_LABEL_MAX_SHIFT},
      {"highlight-terminal", 33},
      {"terminus-highlight-fill", OPT_TERMINUS_HIGHLIGHT_FILL},
      {"terminus-highlight-stroke", OPT_TERMINUS_HIGHLIGHT_STROKE},
      {"compact-terminal-label", 48},
      {"compact-route-label", 49},
      {"no-single-route-labels", OPT_NO_SINGLE_ROUTE_LABELS},
      {"no-render-stations", 7},
      {"labels", 'l'},
      {"route-labels", 'r'},
      {"tight-stations", 9},
      {"render-dir-markers", 10},
      {"render-markers-tail", 20},
      {"dir-marker-spacing", 50},
      {"tail-ignore-sharp-angle", 47},
      {"no-render-node-connections", 11},
      {"resolution", 12},
      {"padding", 13},
      {"padding-top", 23},
      {"padding-right", 24},
      {"padding-bottom", 25},
      {"padding-left", 26},
      {"smoothing", 14},
      {"ratio", 27},
      {"tl-ratio", 31},
      {"render-node-fronts", 15},
      {"crowded-line-thresh", 28},
      {"sharp-turn-angle", 29},
      {"bi-dir-marker", 30},
      {"random-colors", 18},
      {"print-stats", 19},
      {"landmark", 21},
      {"landmarks", 22},
      {"force-landmarks", 51},
      {"landmark-search-radius", 56},
      {"displacement-iterations", 63},
      {"displacement-cooling", 64},
      {"landmarks-webmerc", 54},
      {"me-size", 41},
      {"me-label", 42},
      {"me-star", OPT_ME_STAR},
      {"me", 39},
      {"me-station", 43},
      {"me-with-bg", OPT_ME_WITH_BG},
      {"me-bg-fill", OPT_ME_BG_FILL},
      {"me-bg-stroke", OPT_ME_BG_STROKE},
      {"me-label-color", OPT_ME_LABEL_COLOR},
      {"me-station-fill", 44},
      {"me-station-border", 45},
      {"bg-map", 52},
      {"bg-map-webmerc", 53},
      {"bg-map-opacity", 60},
      {"extend-with-bgmap", 57},
      {"geo-lock", OPT_GEO_LOCK},
      {"geo-lock-bbox", 59}};

  auto parseIni = [&](const std::string &path) {
    std::ifstream in(path.c_str());
    if (!in.good())
      return;
    std::string line;
    while (std::getline(in, line)) {
      auto pos = line.find('#');
      if (pos != std::string::npos)
        line = line.substr(0, pos);
      line = util::trim(line);
      if (line.empty())
        continue;
      pos = line.find('=');
      std::string key = util::trim(line.substr(0, pos));
      std::string val =
          pos == std::string::npos ? "" : util::trim(line.substr(pos + 1));
      auto it = optMap.find(key);
      if (it != optMap.end()) {
        applyOption(cfg, it->second, val);
      }
    }
  };

#ifdef _WIN32
  const char *homeEnv = std::getenv("USERPROFILE");
#else
  const char *homeEnv = std::getenv("HOME");
#endif
  if (homeEnv) {
    parseIni(std::string(homeEnv) + "/.loom.ini");
  }

#ifdef _WIN32
  char buf[MAX_PATH];
  DWORD len = GetModuleFileNameA(NULL, buf, MAX_PATH);
  std::string exe(buf, len);
#else
  char buf[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
  if (len != -1) {
    buf[len] = '\0';
  }
  std::string exe(buf);
#endif
  std::string binDir = dirName(exe);
  parseIni(joinPath(binDir, ".loom.ini"));

  // Check command line for an explicit configuration file before parsing
  // other options so that later flags override config values.
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    std::string path;
    if (arg.rfind("--config=", 0) == 0) {
      path = arg.substr(9);
    } else if (arg == "--config" && i + 1 < argc) {
      path = argv[++i];
    }
    if (!path.empty()) {
      parseIni(path);
    }
  }

  struct option ops[] = {
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"config", required_argument, 0, 46},
      {"line-width", required_argument, 0, 2},
      {"line-spacing", required_argument, 0, 3},
      {"outline-width", required_argument, 0, 4},
      {"log-level", required_argument, 0, 55},
      {"from-dot", no_argument, 0, 'D'},
      {"no-deg2-labels", no_argument, 0, 16},
      {"line-label-textsize", required_argument, 0, 5},
      {"line-label-bend-angle", required_argument, 0, 35},
      {"line-label-length-ratio", required_argument, 0, 36},
      {"station-label-textsize", required_argument, 0, 6},
      {"station-label-angle-steps", required_argument, 0,
       OPT_STATION_LABEL_ANGLE_STEPS},
      {"station-label-angle-step-deg", required_argument, 0,
       OPT_STATION_LABEL_ANGLE_STEP_DEG},
      {"me-label-textsize", required_argument, 0, 40},
      {"font-svg-max", required_argument, 0, 38},
      {"station-line-overlap-penalty", required_argument, 0, 37},
      {"station-line-overlap-per-line", no_argument, 0,
       OPT_STATION_LINE_OVERLAP_PER_LINE},
      {"station-label-far-crowd-radius", required_argument, 0,
       OPT_STATION_LABEL_FAR_CROWD_RADIUS},
      {"station-label-far-crowd-penalty", required_argument, 0,
       OPT_STATION_LABEL_FAR_CROWD_PENALTY},
      {"side-penalty-weight", required_argument, 0, 61},
      {"same-side-penalty", required_argument, 0, 67},
      {"crowding-same-side-scale", required_argument, 0, 69},
      {"reposition-label", required_argument, 0, OPT_REPOSITION_LABEL},
      {"orientation-penalties", required_argument, 0, 62},
      {"terminus-angle-penalty", required_argument, 0,
       OPT_TERMINUS_ANGLE_PENALTY},
      {"cluster-pen-scale", required_argument, 0, 65},
      {"outside-penalty", required_argument, 0, 66},
      {"route-label-gap", required_argument, 0, 32},
      {"route-label-terminus-gap", required_argument, 0, 34},
      {"terminus-label-anchor", required_argument, 0, OPT_TERMINUS_LABEL_ANCHOR},
      {"terminus-label-max-shift", required_argument, 0,
       OPT_TERMINUS_LABEL_MAX_SHIFT},
      {"highlight-terminal", no_argument, 0, 33},
      {"terminus-highlight-fill", required_argument, 0,
       OPT_TERMINUS_HIGHLIGHT_FILL},
      {"terminus-highlight-stroke", required_argument, 0,
       OPT_TERMINUS_HIGHLIGHT_STROKE},
      {"compact-terminal-label", no_argument, 0, 48},
      {"compact-route-label", no_argument, 0, 49},
      {"no-single-route-labels", no_argument, 0, OPT_NO_SINGLE_ROUTE_LABELS},
      {"no-render-stations", no_argument, 0, 7},
      {"labels", no_argument, 0, 'l'},
      {"route-labels", no_argument, 0, 'r'},
      {"tight-stations", no_argument, 0, 9},
      {"render-dir-markers", no_argument, 0, 10},
      {"render-markers-tail", no_argument, 0, 20},
      {"dir-marker-spacing", required_argument, 0, 50},
      {"tail-ignore-sharp-angle", no_argument, 0, 47},
      {"no-render-node-connections", no_argument, 0, 11},
      {"resolution", required_argument, 0, 12},
      {"padding", required_argument, 0, 13},
      {"padding-top", required_argument, 0, 23},
      {"padding-right", required_argument, 0, 24},
      {"padding-bottom", required_argument, 0, 25},
      {"padding-left", required_argument, 0, 26},
      {"smoothing", required_argument, 0, 14},
      {"ratio", required_argument, 0, 27},
      {"tl-ratio", required_argument, 0, 31},
      {"render-node-fronts", no_argument, 0, 15},
      {"crowded-line-thresh", required_argument, 0, 28},
      {"sharp-turn-angle", required_argument, 0, 29},
      {"bi-dir-marker", no_argument, 0, 30},
      {"random-colors", no_argument, 0, 18},
      {"print-stats", no_argument, 0, 19},
      {"landmark", required_argument, 0, 21},
      {"landmarks", required_argument, 0, 22},
      {"force-landmarks", no_argument, 0, 51},
      {"landmark-search-radius", required_argument, 0, 56},
      {"displacement-iterations", required_argument, 0, 63},
      {"displacement-cooling", required_argument, 0, 64},
      {"landmarks-webmerc", no_argument, 0, 54},
      {"me-size", required_argument, 0, 41},
      {"me-label", no_argument, 0, 42},
      {"me-star", optional_argument, 0, OPT_ME_STAR},
      {"me", required_argument, 0, 39},
      {"me-station", required_argument, 0, 43},
      {"me-with-bg", optional_argument, 0, OPT_ME_WITH_BG},
      {"me-bg-fill", required_argument, 0, OPT_ME_BG_FILL},
      {"me-bg-stroke", required_argument, 0, OPT_ME_BG_STROKE},
      {"me-label-color", required_argument, 0, OPT_ME_LABEL_COLOR},
      {"me-station-fill", required_argument, 0, 44},
      {"me-station-border", required_argument, 0, 45},
      {"bg-map", required_argument, 0, 52},
      {"bg-map-webmerc", no_argument, 0, 53},
      {"bg-map-opacity", required_argument, 0, 60},
      {"extend-with-bgmap", no_argument, 0, 57},
      {"geo-lock", optional_argument, 0, OPT_GEO_LOCK},
      {"geo-lock-bbox", required_argument, 0, 59},
      {0, 0, 0, 0}};
  int c;
  while ((c = getopt_long(argc, argv, ":hvlrD", ops, 0)) != -1) {
    if (c == 'h') {
      help(argv[0]);
      exit(0);
    } else if (c == 'v') {
      std::cout << "transitmap - (LOOM " << VERSION_FULL << ")" << std::endl;
      exit(0);
    } else if (c == ':') {
      std::cerr << argv[optind - 1] << " requires an argument" << std::endl;
      exit(1);
    } else if (c == '?') {
      std::cerr << argv[optind - 1] << " option unknown" << std::endl;
      exit(1);
    } else {
      std::string arg = optarg ? std::string(optarg) : "";
      if (c != 46) {
        applyOption(cfg, c, arg);
      }
    }
  }

  if (cfg->lineWidth < 0) {
    std::cerr << "Error: line width " << cfg->lineWidth << " is negative!"
              << std::endl;
    exit(1);
  }

  if (cfg->outlineWidth < 0) {
    std::cerr << "Error: outline width " << cfg->outlineWidth << " is negative!"
              << std::endl;
    exit(1);
  }

  if (cfg->outputPadding < 0) {
    cfg->outputPadding = (cfg->lineWidth + cfg->lineSpacing);
  }

  bool topUnset = cfg->paddingTop < 0;
  bool rightUnset = cfg->paddingRight < 0;
  bool bottomUnset = cfg->paddingBottom < 0;
  bool leftUnset = cfg->paddingLeft < 0;

  if (topUnset && rightUnset && bottomUnset && leftUnset) {
    cfg->paddingTop = cfg->outputPadding;
    cfg->paddingRight = cfg->outputPadding;
    cfg->paddingBottom = cfg->outputPadding;
    cfg->paddingLeft = cfg->outputPadding;
  } else {
    if (topUnset)
      cfg->paddingTop = 0;
    if (rightUnset)
      cfg->paddingRight = 0;
    if (bottomUnset)
      cfg->paddingBottom = 0;
    if (leftUnset)
      cfg->paddingLeft = 0;
  }

  if (cfg->ratio != -1 && cfg->ratio <= 0) {
    std::cerr << "Error: ratio " << cfg->ratio << " is not positive!"
              << std::endl;
    exit(1);
  }
  if (cfg->tlRatio != -1 && cfg->tlRatio <= 0) {
    std::cerr << "Error: tl-ratio " << cfg->tlRatio << " is not positive!"
              << std::endl;
    exit(1);
  }
  if (cfg->stationLabelAngleSteps == 0) {
    std::cerr << "Error: station-label-angle-steps must be positive!"
              << std::endl;
    exit(1);
  }
  if (cfg->stationLabelAngleSteps % 4 != 0) {
    std::cerr << "Error: station-label-angle-steps "
              << cfg->stationLabelAngleSteps
              << " must be divisible by 4 to support terminus labels!"
              << std::endl;
    exit(1);
  }
  if (cfg->stationLabelAngleStepDeg <= 0) {
    std::cerr << "Error: station-label-angle-step-deg "
              << cfg->stationLabelAngleStepDeg << " is not positive!"
              << std::endl;
    exit(1);
  }
}

// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <float.h>
#include <getopt.h>

#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#include "transitmap/_config.h"
#include "transitmap/config/ConfigReader.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/log/Log.h"

using std::exception;
using transitmapper::config::ConfigReader;
using transitmapper::config::Landmark;

static const char *YEAR = &__DATE__[7];
static const char *COPY =
    "University of Freiburg - Chair of Algorithms and Data Structures";
static const char *AUTHORS = "Patrick Brosi <brosi@informatik.uni-freiburg.de>";

// _____________________________________________________________________________
ConfigReader::ConfigReader() {}

// _____________________________________________________________________________
void ConfigReader::help(const char *bin) const {
  std::cout << std::setfill(' ') << std::left << "transitmap (part of LOOM) "
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
            << std::setw(37) << "  --render-engine arg (=svg)"
#ifdef PROTOBUF_FOUND
            << "Render engine, either 'svg' or 'mvt'\n"
#else
            << "Render engine, only 'svg' supported\n"
#endif
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
            << std::setw(37) << "  --bi-dir-marker"
            << "render markers for bidirectional edges\n"
            << std::setw(37) << "  --crowded-line-thresh arg (=3)"
            << "lines on edge to trigger direction marker\n"
            << std::setw(37) << "  --sharp-turn-angle arg (=0.785398)"
            << "turn angle in radians to trigger direction marker\n"
            << std::setw(37) << "  -l [ --labels ]"
            << "render labels\n"
            << std::setw(37) << "  -r [ --route-labels ]"
            << "render route names at line termini\n"
            << std::setw(37) << "  --line-label-textsize arg (=40)"
            << "textsize for line labels\n"
            << std::setw(37) << "  --station-label-textsize arg (=60)"
            << "textsize for station labels\n"
            << std::setw(37) << "  --route-label-gap arg (=5)"
            << "gap between station and route labels\n"
            << std::setw(37) << "  --highlight-terminal"
            << "highlight terminus stations\n"
            << std::setw(37) << "  --no-deg2-labels"
            << "no labels for deg-2 stations\n"
#ifdef PROTOBUF_FOUND
            << std::setw(37) << "  -z [ --zoom ] (=14)"
            << "zoom level to write for MVT tiles, comma separated or range\n"
            << std::setw(37) << "  --mvt-path (=.)"
            << "path for MVT tiles\n\n"
#endif
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
            << std::setw(37) << "  --landmark arg"
            << "add landmark lat,lon or iconPath,lat,lon,size\n"
            << std::setw(37) << "  --landmarks arg"
            << "read landmarks from file, one iconPath,lat,lon,size per line\n"
            << std::setw(37) << "  --print-stats"
            << "write stats to stdout\n";
}

// _____________________________________________________________________________
void ConfigReader::read(Config *cfg, int argc, char **argv) const {
  struct option ops[] = {{"version", no_argument, 0, 'v'},
                         {"help", no_argument, 0, 'h'},
                         {"render-engine", required_argument, 0, 1},
                         {"line-width", required_argument, 0, 2},
                         {"line-spacing", required_argument, 0, 3},
                         {"outline-width", required_argument, 0, 4},
                         {"from-dot", no_argument, 0, 'D'},
                         {"no-deg2-labels", no_argument, 0, 16},
                         {"line-label-textsize", required_argument, 0, 5},
                         {"station-label-textsize", required_argument, 0, 6},
                         {"route-label-gap", required_argument, 0, 32},
                         {"highlight-terminal", no_argument, 0, 33},
                         {"no-render-stations", no_argument, 0, 7},
                         {"labels", no_argument, 0, 'l'},
                         {"route-labels", no_argument, 0, 'r'},
                         {"tight-stations", no_argument, 0, 9},
                         {"render-dir-markers", no_argument, 0, 10},
                         {"render-markers-tail", no_argument, 0, 20},
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
                         {"zoom", required_argument, 0, 'z'},
                         {"mvt-path", required_argument, 0, 17},
                         {"random-colors", no_argument, 0, 18},
                         {"print-stats", no_argument, 0, 19},
                         {"landmark", required_argument, 0, 21},
                         {"landmarks", required_argument, 0, 22},
                         {0, 0, 0, 0}};

  std::string zoom;

  int c;
  while ((c = getopt_long(argc, argv, ":hvlrDz:", ops, 0)) != -1) {
    switch (c) {
    case 'h':
      help(argv[0]);
      exit(0);
    case 'v':
      std::cout << "transitmap - (LOOM " << VERSION_FULL << ")" << std::endl;
      exit(0);
    case 1:
      cfg->renderMethod = optarg;
      break;
    case 2:
      cfg->lineWidth = atof(optarg);
      break;
    case 3:
      cfg->lineSpacing = atof(optarg);
      break;
    case 4:
      cfg->outlineWidth = atof(optarg);
      break;
    case 5:
      cfg->lineLabelSize = atof(optarg);
      break;
    case 6:
      cfg->stationLabelSize = atof(optarg);
      break;
    case 32:
      cfg->routeLabelGap = atof(optarg);
      break;
    case 33:
      cfg->highlightTerminals = true;
      break;
    case 7:
      cfg->renderStations = false;
      break;
    case 'l':
      cfg->renderLabels = true;
      break;
    case 'r':
      cfg->renderRouteLabels = true;
      break;
    case 9:
      cfg->tightStations = true;
      break;
    case 10:
      cfg->renderDirMarkers = true;
      break;
    case 20:
      cfg->renderMarkersTail = true;
      break;
    case 11:
      cfg->renderNodeConnections = false;
      break;
    case 12:
      cfg->outputResolution = atof(optarg);
      break;
    case 13:
      cfg->outputPadding = atof(optarg);
      break;
    case 23:
      cfg->paddingTop = atof(optarg);
      break;
    case 24:
      cfg->paddingRight = atof(optarg);
      break;
    case 25:
      cfg->paddingBottom = atof(optarg);
      break;
    case 26:
      cfg->paddingLeft = atof(optarg);
      break;
    case 14:
      cfg->inputSmoothing = atof(optarg);
      break;
    case 27:
      cfg->ratio = atof(optarg);
      break;
    case 31:
      cfg->tlRatio = atof(optarg);
      if (cfg->paddingRight < 0) cfg->paddingRight = 500;
      if (cfg->paddingBottom < 0) cfg->paddingBottom = 500;
      if (cfg->paddingTop < 0) cfg->paddingTop = 2000;
      if (cfg->paddingLeft < 0) cfg->paddingLeft = 500;
      break;
    case 15:
      cfg->renderNodeFronts = true;
      break;
    case 28:
      cfg->crowdedLineThresh = atoi(optarg);
      break;
    case 29:
      cfg->sharpTurnAngle = atof(optarg);
      break;
    case 30:
      cfg->renderBiDirMarker = true;
      break;
    case 16:
      cfg->dontLabelDeg2 = true;
      break;
    case 17:
      cfg->mvtPath = optarg;
      break;
    case 18:
      cfg->randomColors = true;
      break;
    case 19:
      cfg->writeStats = true;
      break;
    case 21: {
      auto parts = util::split(optarg, ',');
      Landmark lm;

      if (parts.size() == 2) {
        double lat = atof(parts[0].c_str());
        double lon = atof(parts[1].c_str());
        lm.coord = util::geo::latLngToWebMerc<double>(lat, lon);
        lm.size = 200;
      } else if (parts.size() == 4) {
        lm.iconPath = parts[0];
        double lat = atof(parts[1].c_str());
        double lon = atof(parts[2].c_str());
        lm.size = atof(parts[3].c_str());
        lm.coord = util::geo::latLngToWebMerc<double>(lat, lon);
      } else {
        std::cerr << "Error while parsing landmark " << optarg << std::endl;
        exit(1);
      }

      cfg->landmarks.push_back(lm);
      break;
    }
    case 22: {
      std::ifstream infile(optarg);
      if (!infile.good()) {
        std::cerr << "Could not open landmarks file " << optarg << std::endl;
        exit(1);
      }
      std::string line;
      while (std::getline(infile, line)) {
        if (line.empty())
          continue;
        auto parts = util::split(line, ',');
        Landmark lm;
        if (parts.size() == 2) {
          double lat = atof(parts[0].c_str());
          double lon = atof(parts[1].c_str());
          lm.coord = util::geo::latLngToWebMerc<double>(lat, lon);
          lm.size = 200;
        } else if (parts.size() == 4) {
          lm.iconPath = parts[0];
          double lat = atof(parts[1].c_str());
          double lon = atof(parts[2].c_str());
          lm.size = atof(parts[3].c_str());
          lm.coord = util::geo::latLngToWebMerc<double>(lat, lon);
        } else {
          std::cerr << "Error while parsing landmark " << line << std::endl;
          exit(1);
        }
        cfg->landmarks.push_back(lm);
      }
      break;
    }
    case 'D':
      cfg->fromDot = true;
      break;
    case 'z':
      zoom = optarg;
      break;
    case ':':
      std::cerr << argv[optind - 1];
      std::cerr << " requires an argument" << std::endl;
      exit(1);
    case '?':
      std::cerr << argv[optind - 1];
      std::cerr << " option unknown" << std::endl;
      exit(1);
      break;
    default:
      std::cerr << "Error while parsing arguments" << std::endl;
      exit(1);
      break;
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

  for (auto range : util::split(zoom, ',')) {
    util::replaceAll(range, " ", "");
    util::replaceAll(range, "=", "");
    auto parts = util::split(range, '-');
    if (parts.size() > 2) {
      std::cerr << "Error while parsing zoom range" << zoom << std::endl;
      exit(1);
    }

    int from = atoi(parts.front().c_str());
    int to = atoi(parts.back().c_str());

    if (from > to) {
      int a = from;
      from = to;
      to = a;
    }

    if (from < 0 || from > 25 || to < 0 || to > 25) {
      std::cerr << "Error while parsing zoom range" << zoom << std::endl;
      exit(1);
    }

    for (int z = from; z <= to; z++) {
      cfg->mvtZooms.push_back(z);
    }
  }

  if (cfg->mvtZooms.size() == 0)
    cfg->mvtZooms.push_back(14);

  if (cfg->outputPadding < 0) {
    cfg->outputPadding = (cfg->lineWidth + cfg->lineSpacing);
  }
  if (cfg->paddingTop < 0)
    cfg->paddingTop = cfg->outputPadding;
  if (cfg->paddingRight < 0)
    cfg->paddingRight = cfg->outputPadding;
  if (cfg->paddingBottom < 0)
    cfg->paddingBottom = cfg->outputPadding;
  if (cfg->paddingLeft < 0)
    cfg->paddingLeft = cfg->outputPadding;

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
}

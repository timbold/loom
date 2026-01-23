// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <float.h>
#include <getopt.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "transitmap/_config.h"
#include "transitmap/config/ConfigReader.h"
#include "util/log/Log.h"
#include "util/String.h"

using std::exception;
using transitmapper::config::ConfigReader;

static const char* YEAR = &__DATE__[7];
static const char* COPY =
    "University of Freiburg - Chair of Algorithms and Data Structures";
static const char* AUTHORS = "Patrick Brosi <brosi@informatik.uni-freiburg.de>";

namespace {
std::vector<int> parseZoomLevels(const std::string& spec) {
  std::set<int> levels;
  std::stringstream ss(spec);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token = util::trim(token);
    if (token.empty()) continue;
    size_t dash = token.find('-');
    if (dash == std::string::npos) {
      levels.insert(atoi(token.c_str()));
      continue;
    }
    std::string lhs = util::trim(token.substr(0, dash));
    std::string rhs = util::trim(token.substr(dash + 1));
    if (lhs.empty() || rhs.empty()) continue;
    int a = atoi(lhs.c_str());
    int b = atoi(rhs.c_str());
    if (a > b) std::swap(a, b);
    for (int z = a; z <= b; ++z) levels.insert(z);
  }
  return std::vector<int>(levels.begin(), levels.end());
}

bool isPaperValid(const std::string& paper) {
  return paper == "A4" || paper == "A4L" || paper == "A3" || paper == "A3L";
}

bool isCanvasUnitValid(const std::string& unit) {
  return unit == "px" || unit == "mm";
}

transitmapper::config::Config::MergedStopView parseMergedStopView(
    const std::string& raw) {
  std::string s = util::trim(raw);
  if (s == "none") return transitmapper::config::Config::MergedStopView::NONE;
  if (s == "double_ring")
    return transitmapper::config::Config::MergedStopView::DOUBLE_RING;
  if (s == "merge_count")
    return transitmapper::config::Config::MergedStopView::MERGE_COUNT;
  if (s == "simple_cross")
    return transitmapper::config::Config::MergedStopView::SIMPLE_CROSS;

  std::cerr << "Error: invalid --merged-stop-view '" << raw << "' "
            << "(expected none|double_ring|merge_count|simple_cross)"
            << std::endl;
  exit(1);
}
}  // namespace

// _____________________________________________________________________________
ConfigReader::ConfigReader() {}

// _____________________________________________________________________________
void ConfigReader::help(const char* bin) const {
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
            << "Render engine, only 'svg' supported\n"
            << std::setw(37) << "  --line-width arg (=20)"
            << "width of a single transit line\n"
            << std::setw(37) << "  --line-spacing arg (=10)"
            << "spacing between transit lines\n"
            << std::setw(37) << "  --outline-width arg (=1)"
            << "width of line outlines\n"
            << std::setw(37) << "  --render-dir-markers"
            << "render line direction markers\n"
            << std::setw(37) << "  -l [ --labels ]"
            << "render labels\n"
            << std::setw(37) << "  --line-label-textsize arg (=40)"
            << "textsize for line labels\n"
            << std::setw(37) << "  --station-label-textsize arg (=60)"
            << "textsize for station labels\n"
            << std::setw(37) << "  --no-deg2-labels"
            << "no labels for deg-2 stations\n"
            << std::setw(37) << "  --show-merged-stop-members arg (=tooltip)"
            << " off|tooltip|inline|multiline\n"
            << std::setw(37) << "  --merged-stop-view arg (=double_ring)"
            << " none|double_ring|merge_count|simple_cross\n"
            << std::setw(37) << "  --debug-stop-labels arg"
            << "members (alias for multiline)\n"
            << "Misc:\n"
            << std::setw(37) << "  -D [ --from-dot ]"
            << "input is in dot format\n"
            << std::setw(37) << "  --padding arg (=-1)"
            << "padding, -1 for auto\n"
            << std::setw(37) << "  --padding-top arg (=-1)"
            << "top padding, -1 to use --padding\n"
            << std::setw(37) << "  --padding-right arg (=-1)"
            << "right padding, -1 to use --padding\n"
            << std::setw(37) << "  --padding-bottom arg (=-1)"
            << "bottom padding, -1 to use --padding\n"
            << std::setw(37) << "  --padding-left arg (=-1)"
            << "left padding, -1 to use --padding\n"
            << std::setw(37) << "  --smoothing arg (=1)"
            << "input line smoothing\n"
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
            << std::setw(37) << "  --print-stats"
            << "write stats to stdout\n";
  std::cout << "Background:\n"
            << std::setw(37) << "  --mbtiles arg"
            << "path to MBTiles raster background\n"
            << std::setw(37) << "  --paper arg (=A4L)"
            << "paper ratio: A4|A4L|A3|A3L\n"
            << std::setw(37) << "  --canvas-width arg (=1000)"
            << "canvas width in canvas units\n"
            << std::setw(37) << "  --canvas-unit arg (=px)"
            << "canvas units: px|mm\n"
            << std::setw(37) << "  --zoom-levels arg"
            << "override zooms (range or list, e.g. 8-12 or 10,12)\n"
            << std::setw(37) << "  --max-tiles arg (=512)"
            << "max tile count for background mosaic\n"
            << std::setw(37) << "  --oversample arg (=1.25)"
            << "oversample factor for auto zoom\n"
            << std::setw(37) << "  --background-pad-pct arg (=0.03)"
            << "extra padding for background bbox\n"
            << std::setw(37) << "  --background-opacity arg (=1.0)"
            << "background image opacity\n";
}

// _____________________________________________________________________________
void ConfigReader::read(Config* cfg, int argc, char** argv) const {
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
                         {"no-render-stations", no_argument, 0, 7},
                         {"labels", no_argument, 0, 'l'},
                         {"show-merged-stop-members", required_argument, 0, 17},
                         {"merged-stop-view", required_argument, 0, 34},
                         {"debug-stop-labels", required_argument, 0, 33},
                         {"tight-stations", no_argument, 0, 9},
                         {"render-dir-markers", no_argument, 0, 10},
                         {"no-render-node-connections", no_argument, 0, 11},
                         {"resolution", required_argument, 0, 12},
                         {"padding", required_argument, 0, 13},
                         {"smoothing", required_argument, 0, 14},
                         {"render-node-fronts", no_argument, 0, 15},
                         {"random-colors", no_argument, 0, 18},
                         {"print-stats", no_argument, 0, 19},
                         {"padding-top", required_argument, 0, 29},
                         {"padding-right", required_argument, 0, 30},
                         {"padding-bottom", required_argument, 0, 31},
                         {"padding-left", required_argument, 0, 32},
                         {"mbtiles", required_argument, 0, 20},
                         {"paper", required_argument, 0, 21},
                         {"canvas-width", required_argument, 0, 22},
                         {"canvas-unit", required_argument, 0, 23},
                         {"zoom-levels", required_argument, 0, 24},
                         {"max-tiles", required_argument, 0, 25},
                         {"oversample", required_argument, 0, 26},
                         {"background-pad-pct", required_argument, 0, 27},
                         {"background-opacity", required_argument, 0, 28},
                         {0, 0, 0, 0}};

  int c;
  while ((c = getopt_long(argc, argv, ":hvlD", ops, 0)) != -1) {
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
      case 17:
        cfg->showMergedStopMembers = optarg;
        break;
      case 33:
        if (std::string(optarg) == "members") {
          cfg->showMergedStopMembers = "multiline";
        }
        break;
      case 34:
        cfg->mergedStopView = parseMergedStopView(optarg);
        break;
      case 7:
        cfg->renderStations = false;
        break;
      case 'l':
        cfg->renderLabels = true;
        break;
      case 9:
        cfg->tightStations = true;
        break;
      case 10:
        cfg->renderDirMarkers = true;
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
      case 29:
        cfg->outputPaddingTop = atof(optarg);
        break;
      case 30:
        cfg->outputPaddingRight = atof(optarg);
        break;
      case 31:
        cfg->outputPaddingBottom = atof(optarg);
        break;
      case 32:
        cfg->outputPaddingLeft = atof(optarg);
        break;
      case 14:
        cfg->inputSmoothing = atof(optarg);
        break;
      case 15:
        cfg->renderNodeFronts = true;
        break;
      case 16:
        cfg->dontLabelDeg2 = true;
        break;
      case 18:
        cfg->randomColors = true;
        break;
      case 19:
        cfg->writeStats = true;
        break;
      case 20:
        cfg->mbtilesPath = optarg;
        break;
      case 21:
        cfg->paper = optarg;
        break;
      case 22:
        cfg->canvasWidth = atoi(optarg);
        break;
      case 23:
        cfg->canvasUnit = optarg;
        break;
      case 24:
        cfg->zoomLevels = parseZoomLevels(optarg);
        break;
      case 25:
        cfg->maxTiles = atoi(optarg);
        break;
      case 26:
        cfg->oversample = atof(optarg);
        break;
      case 27:
        cfg->backgroundPadPct = atof(optarg);
        break;
      case 28:
        cfg->backgroundOpacity = atof(optarg);
        break;
      case 'D':
        cfg->fromDot = true;
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

  if (cfg->renderMethod != "svg") {
    std::cerr << "Error: render engine " << cfg->renderMethod
              << " is not supported" << std::endl;
    exit(1);
  }

  if (!isPaperValid(cfg->paper)) {
    std::cerr << "Error: paper " << cfg->paper << " is invalid" << std::endl;
    exit(1);
  }

  if (!isCanvasUnitValid(cfg->canvasUnit)) {
    std::cerr << "Error: canvas unit " << cfg->canvasUnit << " is invalid"
              << std::endl;
    exit(1);
  }

  if (cfg->canvasWidth <= 0) {
    std::cerr << "Error: canvas width " << cfg->canvasWidth << " is invalid"
              << std::endl;
    exit(1);
  }

  if (cfg->maxTiles <= 0) {
    std::cerr << "Error: max tiles " << cfg->maxTiles << " is invalid"
              << std::endl;
    exit(1);
  }

  if (cfg->oversample <= 0) {
    std::cerr << "Error: oversample " << cfg->oversample << " is invalid"
              << std::endl;
    exit(1);
  }

  if (cfg->backgroundPadPct < 0) {
    std::cerr << "Error: background pad pct " << cfg->backgroundPadPct
              << " is invalid" << std::endl;
    exit(1);
  }

  if (cfg->backgroundOpacity < 0 || cfg->backgroundOpacity > 1.0) {
    std::cerr << "Error: background opacity " << cfg->backgroundOpacity
              << " is invalid" << std::endl;
    exit(1);
  }

  if (cfg->outputPadding < 0) {
    cfg->outputPadding = (cfg->lineWidth + cfg->lineSpacing);
  }

  if (cfg->outputPaddingTop < 0) cfg->outputPaddingTop = cfg->outputPadding;
  if (cfg->outputPaddingRight < 0) cfg->outputPaddingRight = cfg->outputPadding;
  if (cfg->outputPaddingBottom < 0)
    cfg->outputPaddingBottom = cfg->outputPadding;
  if (cfg->outputPaddingLeft < 0) cfg->outputPaddingLeft = cfg->outputPadding;
}

// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <float.h>
#include <getopt.h>

#include <algorithm>
#include <cctype>
#include <exception>
#include <iostream>
#include <string>

#include "ad/cppgtfs/gtfs/flat/Route.h"
#include "gtfs2graph/_config.h"
#include "gtfs2graph/config/ConfigReader.h"
#include "util/String.h"
#include "util/log/Log.h"

using gtfs2graph::config::CloseCircularRoutesMode;
using gtfs2graph::config::ConfigReader;

using std::exception;
using std::string;
using std::vector;

static const char* YEAR = &__DATE__[7];
static const char* COPY =
    "University of Freiburg - Chair of Algorithms and Data Structures";
static const char* AUTHORS = "Patrick Brosi <brosi@informatik.uni-freiburg.de>";

// _____________________________________________________________________________
ConfigReader::ConfigReader() {}

// _____________________________________________________________________________
void ConfigReader::help(const char* bin) const {
  std::cout
      << std::setfill(' ') << std::left << "gtfs2graph (part of LOOM) "
      << VERSION_FULL << "\n(built " << __DATE__ << " " << __TIME__ << ")"
      << "\n\n(C) 2017-" << YEAR << " " << COPY << "\n"
      << "Authors: " << AUTHORS << "\n\n"
      << "Usage: " << bin << " <GTFS FEED>\n\n"
      << "Allowed options:\n\n"
      << "General:\n"
      << std::setw(36) << "  -v [ --version ]"
      << "print version\n"
      << std::setw(36) << "  -h [ --help ]"
      << "show this help message\n"
      << std::setw(36) << "  -m [ --mots ] arg (=all)"
      << "MOTs to calculate shapes for, comma sep.,\n"
      << std::setw(36) << " "
      << "  either as string "
         "{all, tram | streetcar,\n"
      << std::setw(36) << " "
      << "  subway | metro, rail | train, bus,\n"
      << std::setw(36) << " "
      << "  ferry | boat | ship, cablecar, gondola,\n"
      << std::setw(36) << " "
      << "  funicular, coach, mono-rail | monorail,\n"
      << std::setw(36) << " "
      << "  trolley | trolleybus | trolley-bus} or\n"
      << std::setw(36) << " "
      << "  as GTFS mot codes\n"
      << std::setw(36) << " "
      << "  funicular, coach} or as GTFS mot codes\n"
      << std::setw(36) << "  -p [ --prune-threshold ] arg (=0)"
      << "Threshold for pruning of seldomly occuring\n"
      << std::setw(36) << " " << "  lines, between 0 and 1\n";
  std::cout << "\nCircular closure:\n"
            << std::setw(36)
            << "  --close-circular-routes arg (=never)"
            << "Close circular trips: never, auto, always\n"
            << std::setw(36)
            << "  --circular-max-close-distance-m arg (=5000)"
            << "Max closure distance in meters, 0 to disable\n"
            << std::setw(36)
            << "  --circular-default-speed-kmh arg (=20)"
            << "Fallback speed for closure time\n"
            << std::setw(36) << "  --circular-min-stops arg (=4)"
            << "Min stops for AUTO\n"
            << std::setw(36) << "  --circular-min-unique-stops arg (=3)"
            << "Min unique stops for AUTO\n"
            << std::setw(36)
            << "  --circular-proximity-ratio arg (=0.15)"
            << "AUTO requires dist <= ratio * bbox diagonal\n";
}

// _____________________________________________________________________________
void ConfigReader::read(Config* cfg, int argc, char** argv) const {
  std::string motStr = "all";
  double pruneThreshold = 0;

  enum {
    OPT_CLOSE_CIRCULAR_ROUTES = 1000,
    OPT_CIRCULAR_MAX_CLOSE_DISTANCE,
    OPT_CIRCULAR_DEFAULT_SPEED_KMH,
    OPT_CIRCULAR_MIN_STOPS,
    OPT_CIRCULAR_MIN_UNIQUE_STOPS,
    OPT_CIRCULAR_PROXIMITY_RATIO,
  };

  struct option ops[] = {{"version", no_argument, 0, 'v'},
                         {"help", no_argument, 0, 'h'},
                         {"mots", required_argument, 0, 'm'},
                         {"prune-threshold", required_argument, 0, 'p'},
                         {"close-circular-routes", required_argument, 0,
                          OPT_CLOSE_CIRCULAR_ROUTES},
                         {"circular-max-close-distance-m", required_argument, 0,
                          OPT_CIRCULAR_MAX_CLOSE_DISTANCE},
                         {"circular-default-speed-kmh", required_argument, 0,
                          OPT_CIRCULAR_DEFAULT_SPEED_KMH},
                         {"circular-min-stops", required_argument, 0,
                          OPT_CIRCULAR_MIN_STOPS},
                         {"circular-min-unique-stops", required_argument, 0,
                          OPT_CIRCULAR_MIN_UNIQUE_STOPS},
                         {"circular-proximity-ratio", required_argument, 0,
                          OPT_CIRCULAR_PROXIMITY_RATIO},
                         {0, 0, 0, 0}};

  int c;
  while ((c = getopt_long(argc, argv, ":hvim:p:", ops, 0)) != -1) {
    switch (c) {
      case 'h':
        help(argv[0]);
        exit(0);
      case 'v':
        std::cout << "gtfs2graph - (LOOM " << VERSION_FULL << ")" << std::endl;
        exit(0);
      case 'm':
        motStr = optarg;
        break;
      case 'p':
        pruneThreshold = atof(optarg);
        break;
      case OPT_CLOSE_CIRCULAR_ROUTES: {
        std::string mode = optarg;
        std::transform(mode.begin(), mode.end(), mode.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        if (mode == "never") {
          cfg->closeCircularRoutesMode = CloseCircularRoutesMode::Never;
        } else if (mode == "auto") {
          cfg->closeCircularRoutesMode = CloseCircularRoutesMode::Auto;
        } else if (mode == "always") {
          cfg->closeCircularRoutesMode = CloseCircularRoutesMode::Always;
        } else {
          std::cerr << "Invalid value for --close-circular-routes: " << optarg
                    << std::endl;
          exit(1);
        }
        break;
      }
      case OPT_CIRCULAR_MAX_CLOSE_DISTANCE:
        cfg->circularMaxCloseDistanceM = atof(optarg);
        break;
      case OPT_CIRCULAR_DEFAULT_SPEED_KMH:
        cfg->circularDefaultSpeedKmh = atof(optarg);
        break;
      case OPT_CIRCULAR_MIN_STOPS:
        cfg->circularMinStops = static_cast<size_t>(atoi(optarg));
        break;
      case OPT_CIRCULAR_MIN_UNIQUE_STOPS:
        cfg->circularMinUniqueStops = static_cast<size_t>(atoi(optarg));
        break;
      case OPT_CIRCULAR_PROXIMITY_RATIO:
        cfg->circularProximityRatio = atof(optarg);
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

  if (optind == argc) {
    std::cerr << "No input GTFS feed specified." << std::endl;
    exit(1);
  }

  cfg->inputFeedPath = argv[optind];
  cfg->pruneThreshold = pruneThreshold;

  for (auto sMotStr : util::split(motStr, ',')) {
    for (auto mot :
         ad::cppgtfs::gtfs::flat::Route::getTypesFromString(sMotStr)) {
      cfg->useMots.insert(mot);
    }
  }
}

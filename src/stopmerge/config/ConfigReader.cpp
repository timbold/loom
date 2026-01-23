// Copyright 2016
// University of Freiburg - Chair of Algorithms and Data Structures
// Author: Patrick Brosi

#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <string>

#include "stopmerge/_config.h"
#include "stopmerge/config/ConfigReader.h"

using stopmerge::config::ConfigReader;
using stopmerge::config::StopMergeConfig;

static const char* YEAR = &__DATE__[7];
static const char* COPY =
    "University of Freiburg - Chair of Algorithms and Data Structures";
static const char* AUTHORS = "Patrick Brosi <brosi@informatik.uni-freiburg.de>";

ConfigReader::ConfigReader() {}

void ConfigReader::help(const char* bin) const {
  std::cout << std::setfill(' ') << std::left << "stopmerge (part of LOOM) "
            << VERSION_FULL << "\n(built " << __DATE__ << " " << __TIME__ << ")"
            << "\n\n(C) 2017-" << YEAR << " " << COPY << "\n"
            << "Authors: " << AUTHORS << "\n\n"
            << "Usage: " << bin << " < linegraph.json\n\n"
            << "Allowed options:\n\n"
            << "Core merge toggles:\n"
            << std::setw(46) << "  --merge-stops arg (=never)"
            << "never|auto|always\n"
            << std::setw(46) << "  --merge-stops-parallel-corridors arg (=off)"
            << "off|auto|always\n\n"
            << "Stop-to-edge snapping & clustering:\n"
            << std::setw(46) << "  --merge-stops-snap-dist-m arg (=40)"
            << "snap distance in meters\n"
            << std::setw(46) << "  --merge-stops-radius-m arg (=40)"
            << "merge radius in meters\n"
            << std::setw(46) << "  --merge-stops-chainage-m arg (=40)"
            << "axis chainage threshold in meters\n"
            << std::setw(46) << "  --merge-stops-node-zone-m arg (=15)"
            << "node zone threshold in meters\n"
            << std::setw(46) << "  --merge-stops-corridor-angle-deg arg (=25)"
            << "corridor collinearity in degrees\n"
            << std::setw(46) << "  --merge-stops-max-cluster-size arg (=6)"
            << "max cluster size in auto mode\n\n"
            << "Parallel corridor refinement:\n"
            << std::setw(46) << "  --parallel-pair-max-dist-m arg (=25)"
            << "max min-distance between groups\n"
            << std::setw(46) << "  --parallel-pair-max-angle-deg arg (=20)"
            << "max parallel angle difference\n"
            << std::setw(46) << "  --parallel-pair-min-overlap-ratio arg (=0.35)"
            << "min overlap ratio along axis\n"
            << std::setw(46) << "  --parallel-pair-require-lines arg (=true)"
            << "require shared line ids\n"
            << std::setw(46) << "  --parallel-pair-min-group-length-m arg (=150)"
            << "min corridor length\n"
            << std::setw(46) << "  --parallel-pair-max-median-dist-m arg (=30)"
            << "max median distance safeguard\n\n"
            << "Hub protection:\n"
            << std::setw(46) << "  --hub-radius-m arg (=80)"
            << "hub radius for density guard\n"
            << std::setw(46) << "  --merge-stops-max-local-density arg (=8)"
            << "max station density in hub radius\n"
            << std::setw(46) << "  --hub-auto-require-base-name arg (=true)"
            << "require shared base name at hubs\n\n"
            << "Artifacts / debug:\n"
            << std::setw(46) << "  --merge-stops-debug-csv arg"
            << "write merge decisions CSV\n"
            << std::setw(46) << "  --parallel-pair-debug-csv arg"
            << "write pairing decisions CSV\n"
            << std::setw(46) << "  --stop-merge-map-json arg"
            << "write merge mapping JSON\n\n"
            << "Manual overrides:\n"
            << std::setw(46) << "  --merge-stops-overrides arg"
            << "override JSON file\n\n"
            << "General:\n"
            << std::setw(46) << "  -v [ --version ]"
            << "print version\n"
            << std::setw(46) << "  -h [ --help ]"
            << "show this help message\n";
}

void ConfigReader::read(StopMergeConfig* cfg, int argc, char** argv) const {
  struct option ops[] = {
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"merge-stops", required_argument, 0, 1},
      {"merge-stops-snap-dist-m", required_argument, 0, 2},
      {"merge-stops-radius-m", required_argument, 0, 3},
      {"merge-stops-chainage-m", required_argument, 0, 4},
      {"merge-stops-node-zone-m", required_argument, 0, 5},
      {"merge-stops-corridor-angle-deg", required_argument, 0, 6},
      {"merge-stops-max-cluster-size", required_argument, 0, 7},
      {"merge-stops-parallel-corridors", required_argument, 0, 8},
      {"parallel-pair-max-dist-m", required_argument, 0, 9},
      {"parallel-pair-max-angle-deg", required_argument, 0, 10},
      {"parallel-pair-min-overlap-ratio", required_argument, 0, 11},
      {"parallel-pair-require-lines", required_argument, 0, 12},
      {"parallel-pair-min-group-length-m", required_argument, 0, 13},
      {"parallel-pair-max-median-dist-m", required_argument, 0, 14},
      {"hub-radius-m", required_argument, 0, 15},
      {"merge-stops-max-local-density", required_argument, 0, 16},
      {"hub-auto-require-base-name", required_argument, 0, 17},
      {"merge-stops-debug-csv", required_argument, 0, 18},
      {"parallel-pair-debug-csv", required_argument, 0, 19},
      {"stop-merge-map-json", required_argument, 0, 20},
      {"merge-stops-overrides", required_argument, 0, 21},
      {0, 0, 0, 0}};

  int c;
  while ((c = getopt_long(argc, argv, ":hv", ops, 0)) != -1) {
    switch (c) {
      case 'h':
        help(argv[0]);
        exit(0);
      case 'v':
        std::cout << "stopmerge - (LOOM " << VERSION_FULL << ")" << std::endl;
        exit(0);
      case 1:
        cfg->mergeStops = optarg;
        break;
      case 2:
        cfg->mergeStopsSnapDistM = atof(optarg);
        break;
      case 3:
        cfg->mergeStopsRadiusM = atof(optarg);
        break;
      case 4:
        cfg->mergeStopsChainageM = atof(optarg);
        break;
      case 5:
        cfg->mergeStopsNodeZoneM = atof(optarg);
        break;
      case 6:
        cfg->mergeStopsCorridorAngleDeg = atof(optarg);
        break;
      case 7:
        cfg->mergeStopsMaxClusterSize = static_cast<size_t>(atoi(optarg));
        break;
      case 8:
        cfg->parallelCorridors = optarg;
        break;
      case 9:
        cfg->parallelPairMaxDistM = atof(optarg);
        break;
      case 10:
        cfg->parallelPairMaxAngleDeg = atof(optarg);
        break;
      case 11:
        cfg->parallelPairMinOverlapRatio = atof(optarg);
        break;
      case 12:
        cfg->parallelPairRequireLines =
            (std::string(optarg) == "true" || std::string(optarg) == "1");
        break;
      case 13:
        cfg->parallelPairMinGroupLengthM = atof(optarg);
        break;
      case 14:
        cfg->parallelPairMaxMedianDistM = atof(optarg);
        break;
      case 15:
        cfg->hubRadiusM = atof(optarg);
        break;
      case 16:
        cfg->mergeStopsMaxLocalDensity = static_cast<size_t>(atoi(optarg));
        break;
      case 17:
        cfg->hubAutoRequireBaseName =
            (std::string(optarg) == "true" || std::string(optarg) == "1");
        break;
      case 18:
        cfg->mergeStopsDebugCsv = optarg;
        break;
      case 19:
        cfg->parallelPairDebugCsv = optarg;
        break;
      case 20:
        cfg->stopMergeMapJson = optarg;
        break;
      case 21:
        cfg->mergeStopsOverrides = optarg;
        break;
      case ':':
        std::cerr << argv[optind - 1] << " requires an argument" << std::endl;
        exit(1);
      case '?':
        std::cerr << argv[optind - 1] << " option unknown" << std::endl;
        exit(1);
      default:
        std::cerr << "Error while parsing arguments" << std::endl;
        exit(1);
    }
  }
}

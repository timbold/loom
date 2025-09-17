#include "transitmap/tests/ConfigParseTest.h"

#include "transitmap/config/ConfigReader.h"
#include "util/Misc.h"
#include <cmath>

using transitmapper::config::Config;
using transitmapper::config::ConfigReader;

void ConfigParseTest::run() {
  Config cfg;
  const char* argv[] = {
      "prog",                         "--station-label-far-crowd-radius",
      "12.5",                        "--station-label-far-crowd-penalty",
      "33",                          "--side-penalty-weight", "4.5",
      "--cluster-pen-scale",          "2.0",
      "--outside-penalty",           "7.5",
      "--orientation-penalties",     "1,2,3,4,5,6,7,8",
      "--displacement-iterations",   "5",
      "--displacement-cooling",      "0.5",
      "--same-side-penalty",         "42",
      "--crowding-same-side-scale",  "0.25",
      "--reposition-label",          "2"};
  ConfigReader reader;
  reader.read(&cfg, 23, const_cast<char**>(argv));
  TEST(std::abs(cfg.stationLabelFarCrowdRadius - 12.5) < 1e-9);
  TEST(std::abs(cfg.stationLabelFarCrowdPenalty - 33.0) < 1e-9);
  TEST(std::abs(cfg.sidePenaltyWeight - 4.5) < 1e-9);
  TEST(cfg.orientationPenalties.size(), ==, 8);
  TEST(cfg.orientationPenalties[0], ==, 1);
  TEST(cfg.orientationPenalties[7], ==, 8);
  TEST(cfg.displacementIterations, ==, 5);
  TEST(std::abs(cfg.displacementCooling - 0.5) < 1e-9);
  TEST(std::abs(cfg.clusterPenScale - 2.0) < 1e-9);
  TEST(std::abs(cfg.outsidePenalty - 7.5) < 1e-9);
  TEST(std::abs(cfg.sameSidePenalty - 42) < 1e-9);
  TEST(std::abs(cfg.crowdingSameSideScale - 0.25) < 1e-9);
  TEST(cfg.repositionLabel, ==, 2);
}

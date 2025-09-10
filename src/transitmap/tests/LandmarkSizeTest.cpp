#include <string>

#include "transitmap/tests/LandmarkSizeTest.h"
#include "transitmap/config/TransitMapConfig.h"
#include "shared/rendergraph/Landmark.h"
#include "util/geo/Geo.h"
#include "util/Misc.h"

using transitmapper::config::Config;
using shared::rendergraph::Landmark;

// Forward declaration of the function under test
std::pair<double, double> getLandmarkSizePx(const Landmark &lm,
                                            const Config *cfg);

void LandmarkSizeTest::run() {
  Config cfg;
  cfg.outputResolution = 1.0;
  cfg.stationLabelSize = 60.0;

  Landmark lm;
  lm.coord = util::geo::DPoint(0, 0);
  lm.size = 10.0;

  // ASCII label
  lm.label = "Hello";
  auto dims = getLandmarkSizePx(lm, &cfg);
  double expectedH = lm.size * cfg.outputResolution;
  double expectedW = util::toWStr(lm.label).size() * (expectedH * 0.6);
  TEST(dims.first, ==, expectedW);
  TEST(dims.second, ==, expectedH);

  // Multi-byte UTF-8 label
  lm.label = "\xE4\xBD\xA0\xE5\xA5\xBD";  // "你好"
  dims = getLandmarkSizePx(lm, &cfg);
  expectedW = util::toWStr(lm.label).size() * (expectedH * 0.6);
  TEST(dims.first, ==, expectedW);
  TEST(dims.second, ==, expectedH);
}


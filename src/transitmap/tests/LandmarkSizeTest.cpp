#include <string>

#include "shared/rendergraph/Landmark.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/tests/LandmarkSizeTest.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"

using shared::rendergraph::Landmark;
using transitmapper::config::Config;

// Forward declaration of the function under test
std::pair<double, double> getLandmarkSizePx(const Landmark &lm,
                                            const Config *cfg);

void LandmarkSizeTest::run() {
  Config cfg;
  cfg.outputResolution = 1.0;
  cfg.stationLabelSize = 60.0;

  Landmark lm;
  lm.coord = util::geo::DPoint(0, 0);
  lm.fontSize = 10.0;

  // ASCII label
  lm.label = "Hello";
  auto dims = getLandmarkSizePx(lm, &cfg);
  double expectedH = lm.fontSize;
  double expectedW = util::toWStr(lm.label).size() * (expectedH * 0.6);
  TEST(dims.first, ==, expectedW);
  TEST(dims.second, ==, expectedH);

  // Multi-byte UTF-8 label
  lm.label = "\xE4\xBD\xA0\xE5\xA5\xBD"; // "你好"
  dims = getLandmarkSizePx(lm, &cfg);
  expectedW = util::toWStr(lm.label).size() * (expectedH * 0.6);
  TEST(dims.first, ==, expectedW);
  TEST(dims.second, ==, expectedH);

  // Width-capping when the computed width exceeds the maximum allowed width
  lm.label = "Hello";  // 5 characters
  lm.fontSize = 200.0; // yields width greater than maxWidth
  dims = getLandmarkSizePx(lm, &cfg);
  double maxWidth = cfg.stationLabelSize * cfg.outputResolution * 0.6 * 10.0;
  double w = util::toWStr(lm.label).size() * (lm.fontSize * 0.6);
  double f = maxWidth / w;
  expectedW = maxWidth;
  expectedH = lm.fontSize * f;
  TEST(dims.first, ==, expectedW);
  TEST(dims.second, ==, expectedH);
}

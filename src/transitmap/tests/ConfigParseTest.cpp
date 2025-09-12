#include "transitmap/tests/ConfigParseTest.h"

#include "transitmap/config/ConfigReader.h"
#include "util/Misc.h"
#include <cmath>

using transitmapper::config::Config;
using transitmapper::config::ConfigReader;

void ConfigParseTest::run() {
  Config cfg;
  const char* argv[] = {"prog", "--side-penalty-weight", "4.5",
                        "--orientation-penalties", "1,2,3,4,5,6,7,8"};
  ConfigReader reader;
  reader.read(&cfg, 5, const_cast<char**>(argv));
  TEST(std::abs(cfg.sidePenaltyWeight - 4.5) < 1e-9);
  TEST(cfg.orientationPenalties.size(), ==, 8);
  TEST(cfg.orientationPenalties[0], ==, 1);
  TEST(cfg.orientationPenalties[7], ==, 8);
}

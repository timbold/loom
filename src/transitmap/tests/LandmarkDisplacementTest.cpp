#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include "shared/rendergraph/Landmark.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/output/SvgRenderer.h"
#include "transitmap/tests/LandmarkDisplacementTest.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"

using shared::rendergraph::Landmark;
using shared::rendergraph::RenderGraph;
using transitmapper::config::Config;
using transitmapper::output::SvgRenderer;

void LandmarkDisplacementTest::run() {
  Config cfg;
  cfg.outputResolution = 1.0;
  cfg.renderOverlappingLandmarks = false;
  cfg.displacementIterations = 50;
  cfg.displacementCooling = 0.8;

  std::string iconPath = "test_icon.svg";
  std::ofstream icon(iconPath);
  icon << "<svg width=\"10\" height=\"10\"></svg>";
  icon.close();

  Landmark l1;
  l1.coord = util::geo::DPoint(0, 0);
  l1.iconPath = iconPath;
  l1.size = 10.0;
  Landmark l2 = l1;

  RenderGraph g;
  g.addLandmark(l1);
  g.addLandmark(l2);

  std::ostringstream out;
  SvgRenderer s(&out, &cfg);
  s.print(g);

  std::string svg = out.str();
  size_t pos1 = svg.find("<use");
  size_t pos2 = svg.find("<use", pos1 + 1);
  TEST(pos1 != std::string::npos);
  TEST(pos2 != std::string::npos);

  auto getAttr = [&](size_t pos, const std::string &name) {
    size_t a = svg.find(name + "=\"", pos);
    TEST(a != std::string::npos);
    size_t start = a + name.size() + 2;
    size_t end = svg.find("\"", start);
    return std::stod(svg.substr(start, end - start));
  };

  double x1 = getAttr(pos1, "x");
  double y1 = getAttr(pos1, "y");
  double w1 = getAttr(pos1, "width");
  double h1 = getAttr(pos1, "height");
  double x2 = getAttr(pos2, "x");
  double y2 = getAttr(pos2, "y");
  double w2 = getAttr(pos2, "width");
  double h2 = getAttr(pos2, "height");

  bool overlap = !(x1 + w1 <= x2 || x2 + w2 <= x1 ||
                   y1 + h1 <= y2 || y2 + h2 <= y1);
  TEST(overlap, ==, false);

  std::remove(iconPath.c_str());
}


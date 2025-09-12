#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>

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

  // Now test that displacement is capped by the search radius.
  Config cfg2;
  cfg2.outputResolution = 1.0;
  cfg2.renderOverlappingLandmarks = false;
  cfg2.displacementIterations = 50;
  cfg2.displacementCooling = 0.8;
  cfg2.landmarkSearchRadius = 1;  // allow only one step of movement

  // recreate icon file
  std::ofstream icon2(iconPath);
  icon2 << "<svg width=\"10\" height=\"10\"></svg>";
  icon2.close();

  Landmark l3;
  l3.coord = util::geo::DPoint(0, 0);
  l3.iconPath = iconPath;
  l3.size = 10.0;
  Landmark l4 = l3;

  RenderGraph g2;
  g2.addLandmark(l3);
  g2.addLandmark(l4);

  std::ostringstream out2;
  SvgRenderer s2(&out2, &cfg2);
  s2.print(g2);

  std::string svg2 = out2.str();
  size_t p1 = svg2.find("<use");
  size_t p2 = svg2.find("<use", p1 + 1);
  TEST(p1 != std::string::npos);
  TEST(p2 != std::string::npos);

  auto getAttr2 = [&](size_t pos, const std::string &name) {
    size_t a = svg2.find(name + "=\"", pos);
    TEST(a != std::string::npos);
    size_t start = a + name.size() + 2;
    size_t end = svg2.find("\"", start);
    return std::stod(svg2.substr(start, end - start));
  };

  double ax1 = getAttr2(p1, "x");
  double ay1 = getAttr2(p1, "y");
  double aw1 = getAttr2(p1, "width");
  double ah1 = getAttr2(p1, "height");
  double ax2 = getAttr2(p2, "x");
  double ay2 = getAttr2(p2, "y");
  double aw2 = getAttr2(p2, "width");
  double ah2 = getAttr2(p2, "height");

  double c1x = ax1 + aw1 / 2.0;
  double c1y = ay1 + ah1 / 2.0;
  double c2x = ax2 + aw2 / 2.0;
  double c2y = ay2 + ah2 / 2.0;
  double disp = std::sqrt((c2x - c1x) * (c2x - c1x) +
                          (c2y - c1y) * (c2y - c1y));
  double maxDisp = cfg2.landmarkSearchRadius *
                   (std::max(aw1, ah1) / 2.0);
  TEST(disp <= maxDisp + 1e-6);

  // With such a small radius the landmarks still overlap.
  bool overlap2 = !(ax1 + aw1 <= ax2 || ax2 + aw2 <= ax1 ||
                    ay1 + ah1 <= ay2 || ay2 + ah2 <= ay1);
  TEST(overlap2, ==, true);

  std::remove(iconPath.c_str());
}


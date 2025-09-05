#include <cmath>
#include <sstream>

#define private public
#include "transitmap/output/SvgRenderer.h"
#undef private

#include "transitmap/config/ConfigReader.h"
#include "transitmap/tests/LandmarkProjectionTest.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"

using transitmapper::config::Config;
using transitmapper::config::ConfigReader;
using transitmapper::output::SvgRenderer;
using shared::rendergraph::RenderGraph;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNode;
using shared::linegraph::LineNodePL;
using shared::linegraph::NodeFront;
using util::geo::DPoint;
using util::geo::PolyLine;

namespace {
LineEdge* makeEdge(RenderGraph& g, const Line* line, double length) {
  LineNode* a = g.addNd(LineNodePL(DPoint(0, 0)));
  LineNode* b = g.addNd(LineNodePL(DPoint(length, 0)));
  PolyLine<double> pl;
  pl << DPoint(0, 0) << DPoint(length, 0);
  LineEdgePL epl(pl);
  LineEdge* e = g.addEdg(a, b, epl);
  e->pl().addLine(line, b);
  NodeFront nfA(a, e);
  PolyLine<double> nfAP;
  nfAP << DPoint(0, 0) << DPoint(0, 1);
  nfA.setInitialGeom(nfAP);
  a->pl().addFront(nfA);
  NodeFront nfB(b, e);
  PolyLine<double> nfBP;
  nfBP << DPoint(length, 0) << DPoint(length, 1);
  nfB.setInitialGeom(nfBP);
  b->pl().addFront(nfB);
  return e;
}
}  // namespace

void LandmarkProjectionTest::run() {
  Config cfg;
  const char* argv[] = {"prog", "--landmark", "word:Test,0.0001,0.0001,1"};
  ConfigReader reader;
  reader.read(&cfg, 3, const_cast<char**>(argv));

  RenderGraph g;
  Line line("L1", "L1", "#000");
  makeEdge(g, &line, 10.0);
  for (const auto& lm : cfg.landmarks) {
    g.addLandmark(lm);
  }

  std::ostringstream out;
  SvgRenderer s(&out, &cfg);
  s.print(g);

  std::string svg = out.str();
  size_t pos = svg.find("<text");
  TEST(pos != std::string::npos, ==, true);
  size_t xPos = svg.find("x=\"", pos);
  size_t xEnd = svg.find("\"", xPos + 3);
  double x = std::stod(svg.substr(xPos + 3, xEnd - (xPos + 3)));
  size_t yPos = svg.find("y=\"", pos);
  size_t yEnd = svg.find("\"", yPos + 3);
  double y = std::stod(svg.substr(yPos + 3, yEnd - (yPos + 3)));

  // Compute expected coordinates in the SVG after Web Mercator conversion.
  util::geo::DPoint merc =
      util::geo::latLngToWebMerc(DPoint(0.0001, 0.0001));
  util::geo::Box<double> box = g.getBBox();
  box = util::geo::pad(box, g.getMaxLineNum() *
                                (cfg.lineWidth + cfg.lineSpacing));
  double xOff = box.getLowerLeft().getX();
  double yOff = box.getLowerLeft().getY();
  double svgH = (box.getUpperRight().getY() - yOff) * cfg.outputResolution;
  double expX = (merc.getX() - xOff) * cfg.outputResolution;
  double expY = svgH - (merc.getY() - yOff) * cfg.outputResolution;

  TEST(std::abs(x - expX) < 1e-6);
  TEST(std::abs(y - expY) < 1e-6);
}


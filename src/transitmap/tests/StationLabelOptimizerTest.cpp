#include "transitmap/tests/StationLabelOptimizerTest.h"

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "shared/linegraph/Line.h"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/label/Labeller.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"

using shared::linegraph::Line;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNode;
using shared::linegraph::LineNodePL;
using shared::rendergraph::RenderGraph;
using transitmapper::config::Config;
using transitmapper::label::Labeller;

namespace {

double computeSide(const LineNode* a, const LineNode* b, size_t deg,
                   double stepDeg) {
  double edgeVecX = b->pl().getGeom()->getX() - a->pl().getGeom()->getX();
  double edgeVecY = b->pl().getGeom()->getY() - a->pl().getGeom()->getY();
  double ang = deg * stepDeg * M_PI / 180.0;
  double vecX = std::cos(ang);
  double vecY = std::sin(ang);
  return edgeVecX * vecY - edgeVecY * vecX;
}

}  // namespace

void StationLabelOptimizerTest::run() {
  Config cfg;
  cfg.outputResolution = 1.0;
  cfg.sameSidePenalty = 500.0;
  cfg.sidePenaltyWeight = 10.0;

  RenderGraph g;
  auto* n1 = g.addNd(LineNodePL(util::geo::DPoint(0.0, 0.0)));
  auto* n2 = g.addNd(LineNodePL(util::geo::DPoint(60.0, 0.0)));

  std::vector<std::unique_ptr<Line>> lines;
  lines.emplace_back(new Line("L", "L", "#000"));

  util::geo::PolyLine<double> edgeGeom;
  edgeGeom << *n1->pl().getGeom() << *n2->pl().getGeom();
  auto* edge = g.addEdg(n1, n2, LineEdgePL(edgeGeom));
  edge->pl().addLine(lines.back().get(), n2);

  auto stop1 = shared::linegraph::Station("A", "A", *n1->pl().getGeom());
  stop1.pos = util::geo::DPoint(n1->pl().getGeom()->getX(),
                                 n1->pl().getGeom()->getY() + 10.0);
  n1->pl().addStop(stop1);

  auto stop2 = shared::linegraph::Station("B", "B", *n2->pl().getGeom());
  stop2.pos = util::geo::DPoint(n2->pl().getGeom()->getX(),
                                 n2->pl().getGeom()->getY() - 10.0);
  n2->pl().addStop(stop2);

  Labeller labeller(&cfg);
  labeller.label(g, false);

  size_t deg1 = n1->pl().stops()[0].labelDeg;
  size_t deg2 = n2->pl().stops()[0].labelDeg;

  TEST(deg1 != std::numeric_limits<size_t>::max());
  TEST(deg2 != std::numeric_limits<size_t>::max());

  size_t angleSteps = cfg.stationLabelAngleSteps;
  double angleStepDeg = cfg.stationLabelAngleStepDeg;
  size_t deg1Mod = deg1 % angleSteps;
  size_t deg2Mod = deg2 % angleSteps;

  double side1 = computeSide(n1, n2, deg1Mod, angleStepDeg);
  double side2 = computeSide(n2, n1, deg2Mod, angleStepDeg);
  TEST(side1 * side2 >= 0.0);

  TEST(deg1Mod == deg2Mod);
  TEST(deg2Mod != 3 * angleSteps / 4);
}

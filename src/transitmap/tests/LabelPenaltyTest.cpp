#include <memory>
#include <vector>

#include "transitmap/tests/LabelPenaltyTest.h"

#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/label/Labeller.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"
#include "util/geo/Line.h"

using transitmapper::config::Config;
using transitmapper::label::Labeller;
using transitmapper::label::StationLabel;
using transitmapper::label::Overlaps;
using util::geo::PolyLine;
using util::geo::MultiLine;
using util::geo::DPoint;
using shared::linegraph::Line;
using shared::linegraph::LineEdgePL;
using shared::rendergraph::RenderGraph;
using shared::linegraph::Station;

namespace transitmapper::label {

class LabellerOverlapTestAccess {
 public:
  static Overlaps getOverlaps(Labeller& labeller,
                              const util::geo::MultiLine<double>& band,
                              const shared::linegraph::LineNode* node,
                              const shared::rendergraph::RenderGraph& g,
                              double radius) {
    return labeller.getOverlaps(band, node, g, radius);
  }
};

}  // namespace transitmapper::label

namespace {

const shared::linegraph::LineNode* buildOverlapScenario(
    RenderGraph* g, std::vector<std::unique_ptr<Line>>* lines,
    MultiLine<double>* band) {
  lines->clear();
  band->clear();

  auto* nodeA = g->addNd(shared::linegraph::LineNodePL(DPoint(0.0, 0.0)));
  auto* nodeB = g->addNd(shared::linegraph::LineNodePL(DPoint(100.0, 0.0)));
  auto* nodeC = g->addNd(shared::linegraph::LineNodePL(DPoint(200.0, 0.0)));

  PolyLine<double> geomAB;
  geomAB << *nodeA->pl().getGeom() << *nodeB->pl().getGeom();
  auto* edgeAB = g->addEdg(nodeA, nodeB, LineEdgePL(geomAB));

  PolyLine<double> geomBC;
  geomBC << *nodeB->pl().getGeom() << *nodeC->pl().getGeom();
  auto* edgeBC = g->addEdg(nodeB, nodeC, LineEdgePL(geomBC));

  lines->emplace_back(new Line("L1", "Line 1", "#000"));
  lines->emplace_back(new Line("L2", "Line 2", "#111"));
  lines->emplace_back(new Line("L3", "Line 3", "#222"));
  lines->emplace_back(new Line("L4", "Line 4", "#333"));

  edgeAB->pl().addLine(lines->at(0).get(), nodeB);
  edgeAB->pl().addLine(lines->at(1).get(), nodeB);
  edgeAB->pl().addLine(lines->at(2).get(), nodeB);

  edgeBC->pl().addLine(lines->at(0).get(), nodeC);
  edgeBC->pl().addLine(lines->at(3).get(), nodeC);

  util::geo::Line<double> base;
  base.push_back(*nodeA->pl().getGeom());
  base.push_back(*nodeC->pl().getGeom());
  band->push_back(base);
  band->push_back(base);

  return nodeA;
}

}  // namespace

void LabelPenaltyTest::run() {
  Config cfgA;
  Config cfgB;
  cfgA.sidePenaltyWeight = 1.0;
  cfgB.sidePenaltyWeight = 5.0;
  cfgA.orientationPenalties = {0, 100, 0, 0, 0, 0, 0, 0};
  cfgB.orientationPenalties = cfgA.orientationPenalties;

  PolyLine<double> geom;
  MultiLine<double> band;
  std::vector<const shared::linegraph::Line*> lines;
  Overlaps ov{0, 0, 0, 0, 0};
  Station st("id", "name", util::geo::DPoint());

  double sideA = 2 * cfgA.sidePenaltyWeight;
  double sideB = 2 * cfgB.sidePenaltyWeight;

  StationLabel lblA(geom, band, 10, false, lines, 0, 0, ov, sideA,
                    cfgA.stationLineOverlapPenalty, 0, 0, false,
                    cfgA.clusterPenScale, cfgA.outsidePenalty,
                    &cfgA.orientationPenalties, st);
  StationLabel lblB(geom, band, 10, false, lines, 0, 0, ov, sideB,
                    cfgB.stationLineOverlapPenalty, 0, 0, false,
                    cfgB.clusterPenScale, cfgB.outsidePenalty,
                    &cfgB.orientationPenalties, st);
  TEST(lblB.getPen() > lblA.getPen());

  StationLabel lblC(geom, band, 10, false, lines, 1, 0, ov, sideA,
                    cfgA.stationLineOverlapPenalty, 0, 0, false,
                    cfgA.clusterPenScale, cfgA.outsidePenalty,
                    &cfgA.orientationPenalties, st);
  TEST(lblC.getPen() > lblA.getPen());

  StationLabel lblFar(geom, band, 10, false, lines, 0, 0, ov, sideA,
                      cfgA.stationLineOverlapPenalty, 0,
                      cfgA.stationLabelFarCrowdPenalty, false,
                      cfgA.clusterPenScale, cfgA.outsidePenalty,
                      &cfgA.orientationPenalties, st);
  TEST(lblFar.getPen() > lblA.getPen());

  // crowding penalty scale
  Config cfgC;
  Config cfgD;
  cfgC.clusterPenScale = 1.0;
  cfgD.clusterPenScale = 5.0;
  StationLabel lblD(geom, band, 10, false, lines, 0, 0, ov, 0,
                    cfgC.stationLineOverlapPenalty, 1.0, 0, false,
                    cfgC.clusterPenScale, cfgC.outsidePenalty,
                    &cfgC.orientationPenalties, st);
  StationLabel lblE(geom, band, 10, false, lines, 0, 0, ov, 0,
                    cfgD.stationLineOverlapPenalty, 1.0, 0, false,
                    cfgD.clusterPenScale, cfgD.outsidePenalty,
                    &cfgD.orientationPenalties, st);
  TEST(lblE.getPen() > lblD.getPen());

  // outside penalty bonus/penalty
  Config cfgOutPen;
  Config cfgOutBon;
  cfgOutPen.outsidePenalty = 5.0;
  cfgOutBon.outsidePenalty = -5.0;
  StationLabel lblOutPen(geom, band, 10, false, lines, 0, 0, ov, 0,
                         cfgOutPen.stationLineOverlapPenalty, 0, 0, true,
                         cfgOutPen.clusterPenScale, cfgOutPen.outsidePenalty,
                         &cfgOutPen.orientationPenalties, st);
  StationLabel lblOutBon(geom, band, 10, false, lines, 0, 0, ov, 0,
                         cfgOutBon.stationLineOverlapPenalty, 0, 0, true,
                         cfgOutBon.clusterPenScale, cfgOutBon.outsidePenalty,
                         &cfgOutBon.orientationPenalties, st);
  TEST(lblOutPen.getPen() > lblOutBon.getPen());

  // Station line overlap counting per edge (default).
  {
    Config cfg;
    cfg.stationLineOverlapPerLine = false;
    RenderGraph g;
    std::vector<std::unique_ptr<Line>> lines;
    MultiLine<double> overlapBand;
    auto* anchor = buildOverlapScenario(&g, &lines, &overlapBand);
    Labeller labeller(&cfg);
    auto overlaps = transitmapper::label::LabellerOverlapTestAccess::getOverlaps(
        labeller, overlapBand, anchor, g, 150.0);
    TEST(overlaps.lineOverlaps, ==, 2);
  }

  // Station line overlap counting per distinct line when enabled.
  {
    Config cfg;
    cfg.stationLineOverlapPerLine = true;
    RenderGraph g;
    std::vector<std::unique_ptr<Line>> lines;
    MultiLine<double> overlapBand;
    auto* anchor = buildOverlapScenario(&g, &lines, &overlapBand);
    Labeller labeller(&cfg);
    auto overlaps = transitmapper::label::LabellerOverlapTestAccess::getOverlaps(
        labeller, overlapBand, anchor, g, 150.0);
    TEST(overlaps.lineOverlaps, ==, 4);
  }
}

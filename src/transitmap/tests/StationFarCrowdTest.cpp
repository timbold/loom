#include "transitmap/tests/StationFarCrowdTest.h"

#include <cmath>
#include <memory>
#include <vector>

#include "shared/linegraph/Line.h"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/label/Labeller.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"

using shared::linegraph::Line;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNode;
using shared::linegraph::LineNodePL;
using shared::rendergraph::Landmark;
using shared::rendergraph::RenderGraph;
using transitmapper::config::Config;
using transitmapper::label::Labeller;
using transitmapper::label::Overlaps;
using transitmapper::label::StationLabel;
using util::geo::DPoint;
using util::geo::MultiLine;
using util::geo::PolyLine;

namespace transitmapper {
namespace label {

class LabellerFarCrowdTestAccess {
 public:
  static Labeller::StationCrowdContext computeContext(
      Labeller& labeller, const util::geo::MultiLine<double>& band,
      const shared::linegraph::LineNode* node, double searchRadius,
      const shared::rendergraph::RenderGraph& g) {
    return labeller.computeStationFarCrowd(band, node, searchRadius, g);
  }

  static size_t addSyntheticLabel(
      Labeller& labeller, const util::geo::MultiLine<double>& band,
      double fontSize) {
    Overlaps overlaps{0, 0, 0, 0, 0};
    shared::linegraph::Station st("synthetic", "synthetic",
                                  util::geo::DPoint());
    StationLabel label(util::geo::PolyLine<double>(band[0]), band, fontSize,
                       false, 0, 0, overlaps, 0,
                       labeller._cfg->stationLineOverlapPenalty, 0, 0, false,
                       labeller._cfg->clusterPenScale,
                       labeller._cfg->outsidePenalty, nullptr, st);
    labeller._stationLabels.push_back(label);
    size_t idx = labeller._stationLabels.size() - 1;
    labeller._statLblIdx.add(band, idx);
    return idx;
  }
};

}  // namespace label
}  // namespace transitmapper

using transitmapper::label::LabellerFarCrowdTestAccess;

namespace {

MultiLine<double> makeTestBand() {
  MultiLine<double> band;
  util::geo::Line<double> base;
  base.push_back(DPoint(0.0, 0.0));
  base.push_back(DPoint(80.0, 0.0));
  band.push_back(base);
  band.push_back(base);
  return band;
}

LineNode* addCentralStation(RenderGraph* g) {
  auto* node = g->addNd(LineNodePL(DPoint(0.0, 0.0)));
  node->pl().addStop(shared::linegraph::Station("central", "central",
                                               *node->pl().getGeom()));
  return node;
}

}  // namespace

void StationFarCrowdTest::run() {
  Config cfg;
  cfg.outputResolution = 1.0;
  cfg.stationLabelFarCrowdRadius = 30.0;
  cfg.stationLabelFarCrowdPenalty = 10.0;

  double searchRadius = 200.0;

  // Baseline: no nearby features.
  {
    RenderGraph g;
    auto* station = addCentralStation(&g);
    auto band = makeTestBand();
    Labeller labeller(&cfg);
    auto ctx = LabellerFarCrowdTestAccess::computeContext(labeller, band,
                                                          station, searchRadius, g);
    TEST(ctx.farCrowdPen, ==, 0.0);
  }

  // Edge near the far label end.
  {
    RenderGraph g;
    auto* station = addCentralStation(&g);
    auto band = makeTestBand();
    auto farPoint = band[1].back();

    std::vector<std::unique_ptr<Line>> lines;
    lines.emplace_back(new Line("edge", "edge", "#000"));
    auto* edgeA =
        g.addNd(LineNodePL(DPoint(farPoint.getX(), farPoint.getY() + 5.0)));
    auto* edgeB = g.addNd(
        LineNodePL(DPoint(farPoint.getX() + 20.0, farPoint.getY() + 5.0)));
    PolyLine<double> edgeGeom;
    edgeGeom << *edgeA->pl().getGeom() << *edgeB->pl().getGeom();
    auto* edge = g.addEdg(edgeA, edgeB, LineEdgePL(edgeGeom));
    edge->pl().addLine(lines.back().get(), edgeB);

    Labeller labeller(&cfg);
    auto ctx = LabellerFarCrowdTestAccess::computeContext(labeller, band,
                                                          station, searchRadius, g);
    TEST(std::abs(ctx.farCrowdPen - cfg.stationLabelFarCrowdPenalty) < 1e-9);
  }

  // Synthetic label positioned near the far end.
  {
    RenderGraph g;
    auto* station = addCentralStation(&g);
    auto band = makeTestBand();
    auto farPoint = band[1].back();

    MultiLine<double> neighborBand;
    util::geo::Line<double> neighbor;
    neighbor.push_back(DPoint(farPoint.getX() - 5.0, farPoint.getY() + 4.0));
    neighbor.push_back(DPoint(farPoint.getX() + 5.0, farPoint.getY() + 4.0));
    neighborBand.push_back(neighbor);
    neighborBand.push_back(neighbor);

    Labeller labeller(&cfg);
    LabellerFarCrowdTestAccess::addSyntheticLabel(labeller, neighborBand,
                                                  cfg.stationLabelSize);
    auto ctx = LabellerFarCrowdTestAccess::computeContext(labeller, band,
                                                          station, searchRadius, g);
    TEST(std::abs(ctx.farCrowdPen - cfg.stationLabelFarCrowdPenalty) < 1e-9);
  }

  // Station hull from a nearby stop.
  {
    RenderGraph g;
    auto* station = addCentralStation(&g);
    auto band = makeTestBand();
    auto farPoint = band[1].back();

    std::vector<std::unique_ptr<Line>> lines;
    lines.emplace_back(new Line("hull", "hull", "#000"));
    auto* hullNode =
        g.addNd(LineNodePL(DPoint(farPoint.getX() + 5.0, farPoint.getY())));
    auto* hullOther = g.addNd(
        LineNodePL(DPoint(farPoint.getX() + 35.0, farPoint.getY())));
    PolyLine<double> hullGeom;
    hullGeom << *hullNode->pl().getGeom() << *hullOther->pl().getGeom();
    auto* hullEdge = g.addEdg(hullNode, hullOther, LineEdgePL(hullGeom));
    hullEdge->pl().addLine(lines.back().get(), hullOther);

    Labeller labeller(&cfg);
    auto ctxEdgeOnly = LabellerFarCrowdTestAccess::computeContext(
        labeller, band, station, searchRadius, g);
    TEST(std::abs(ctxEdgeOnly.farCrowdPen - cfg.stationLabelFarCrowdPenalty) <
         1e-9);

    hullNode->pl().addStop(shared::linegraph::Station(
        "hull", "hull", *hullNode->pl().getGeom()));
    auto ctxWithHull = LabellerFarCrowdTestAccess::computeContext(
        labeller, band, station, searchRadius, g);
    TEST(std::abs(ctxWithHull.farCrowdPen -
                  2.0 * cfg.stationLabelFarCrowdPenalty) < 1e-9);
  }

  // Landmark close to the far end.
  {
    RenderGraph g;
    auto* station = addCentralStation(&g);
    auto band = makeTestBand();
    auto farPoint = band[1].back();

    Landmark lm;
    lm.coord = DPoint(farPoint.getX(), farPoint.getY() + 6.0);
    lm.label = "LM";
    lm.fontSize = 12.0;
    lm.size = 12.0;
    g.addLandmark(lm);

    Labeller labeller(&cfg);
    auto ctx = LabellerFarCrowdTestAccess::computeContext(labeller, band,
                                                          station, searchRadius, g);
    TEST(std::abs(ctx.farCrowdPen - cfg.stationLabelFarCrowdPenalty) < 1e-9);
  }
}

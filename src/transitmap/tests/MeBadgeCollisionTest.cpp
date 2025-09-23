#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#define private public
#include "transitmap/output/SvgRenderer.h"
#include "transitmap/label/Labeller.h"
#undef private

#include "shared/linegraph/Line.h"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/tests/MeBadgeCollisionTest.h"
#include "transitmap/util/String.h"
#include "util/Misc.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/geo/Line.h"

using shared::linegraph::Line;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNodePL;
using shared::linegraph::Station;
using shared::rendergraph::RenderGraph;
using transitmapper::config::Config;
using transitmapper::label::Labeller;
using transitmapper::label::StationLabel;
using transitmapper::output::RenderParams;
using transitmapper::output::SvgRenderer;
using util::geo::DPoint;
using util::geo::PolyLine;

namespace {

double bandWidth(const util::geo::MultiLine<double> &band) {
  if (band.empty() || band[0].size() < 2) return 0.0;
  return band[0].back().getX() - band[0].front().getX();
}

double bandHeight(const util::geo::MultiLine<double> &band) {
  if (band.size() < 3 || band[0].empty() || band[2].empty()) return 0.0;
  return band[2].front().getY() - band[0].front().getY();
}

}  // namespace

void MeBadgeCollisionTest::run() {
  Config baseCfg;
  baseCfg.outputResolution = 1.0;
  baseCfg.lineWidth = 4.0;
  baseCfg.lineSpacing = 2.0;
  baseCfg.stationLabelSize = 20.0;
  baseCfg.meStarSize = 10.0;

  RenderGraph gBand;
  auto *bandNode = gBand.addNd(LineNodePL(DPoint(0.0, 0.0)));
  bandNode->pl().addStop(Station("id", "Here", *bandNode->pl().getGeom()));

  Config controlCfg = baseCfg;
  controlCfg.meStation.clear();
  controlCfg.meStationWithBg = false;
  controlCfg.highlightMeStationLabel = false;
  Labeller controlLabeller(&controlCfg);
  auto baseBand = controlLabeller.getStationLblBand(
      bandNode, baseCfg.stationLabelSize, 0, gBand);

  Config badgeCfg = baseCfg;
  badgeCfg.meStation = util::sanitizeStationLabel("Here");
  badgeCfg.meStationWithBg = true;
  badgeCfg.highlightMeStationLabel = true;
  Labeller badgeLabeller(&badgeCfg);
  auto badgeBand = badgeLabeller.getStationLblBand(
      bandNode, badgeCfg.stationLabelSize, 0, gBand);

  TEST(!baseBand.empty());
  TEST(!badgeBand.empty());
  double baseStart = baseBand[0].front().getX();
  double badgeStart = badgeBand[0].front().getX();
  double baseWidth = bandWidth(baseBand);
  double badgeWidth = bandWidth(badgeBand);
  double baseH = bandHeight(baseBand);
  double badgeH = bandHeight(badgeBand);

  double textHeightForPadding = std::max(baseH, badgeCfg.meStarSize);
  double expectedLeftShift = textHeightForPadding * 0.6 +
                             badgeCfg.meStarSize * 0.2 + badgeCfg.meStarSize;
  double expectedWidthIncrease = textHeightForPadding * 1.2 +
                                 badgeCfg.meStarSize +
                                 badgeCfg.meStarSize * 0.2;
  double expectedHeight = std::max(baseH, badgeCfg.meStarSize) +
                          textHeightForPadding * 0.28 +
                          textHeightForPadding * 0.12;

  TEST(std::abs((baseStart - badgeStart) - expectedLeftShift) < 1e-6);
  TEST(std::abs((badgeWidth - baseWidth) - expectedWidthIncrease) < 1e-6);
  TEST(std::abs((badgeH - baseH) - (expectedHeight - baseH)) < 1e-6);

  Config solverCfg = badgeCfg;
  solverCfg.renderLabels = true;
  RenderGraph gSolver;
  std::vector<std::unique_ptr<Line>> lines;
  lines.emplace_back(new Line("L1", "Line 1", "#000000"));
  auto *center = gSolver.addNd(LineNodePL(DPoint(0.0, 0.0)));
  center->pl().addStop(Station("center", "Here", *center->pl().getGeom()));
  auto *east = gSolver.addNd(LineNodePL(DPoint(100.0, 0.0)));
  PolyLine eastGeom;
  eastGeom << *center->pl().getGeom() << *east->pl().getGeom();
  auto *edge = gSolver.addEdg(center, east, LineEdgePL(eastGeom));
  edge->pl().addLine(lines.back().get(), east);

  Labeller solverLabeller(&solverCfg);
  solverLabeller.label(gSolver, false);

  const StationLabel *meLabel = nullptr;
  for (const auto &lbl : solverLabeller.getStationLabels()) {
    if (util::sanitizeStationLabel(lbl.s.name) == solverCfg.meStation) {
      meLabel = &lbl;
      break;
    }
  }
  TEST(meLabel != nullptr);
  util::geo::Box<double> badgeBox =
      util::geo::extendBox(meLabel->band, util::geo::Box<double>());
  double sampleX = badgeBox.getLowerLeft().getX() +
                   (badgeBox.getUpperRight().getX() -
                    badgeBox.getLowerLeft().getX()) * 0.1;
  double sampleY = badgeBox.getLowerLeft().getY() +
                   (badgeBox.getUpperRight().getY() -
                    badgeBox.getLowerLeft().getY()) * 0.5;
  util::geo::Box<double> probe(util::geo::DPoint(sampleX, sampleY),
                               util::geo::DPoint(sampleX + 1.0, sampleY + 1.0));
  TEST(solverLabeller.collidesWithLabels(probe));
  util::geo::Box<double> outside(
      util::geo::DPoint(badgeBox.getLowerLeft().getX() - 5.0, sampleY),
      util::geo::DPoint(badgeBox.getLowerLeft().getX() - 4.0, sampleY + 1.0));
  TEST(!solverLabeller.collidesWithLabels(outside));

  Config highlightCfg = solverCfg;
  highlightCfg.renderMe = true;
  highlightCfg.renderMeLabel = true;
  highlightCfg.displacementIterations = 0;
  Labeller highlightLabeller(&highlightCfg);
  highlightLabeller.label(gSolver, false);

  std::ostringstream highlightOut;
  SvgRenderer highlightRenderer(&highlightOut, &highlightCfg);
  RenderParams params;
  params.xOff = 0.0;
  params.yOff = 0.0;
  params.width = 400.0;
  params.height = 400.0;
  highlightRenderer.renderStationLabels(highlightLabeller, params);
  size_t landmarksBefore = highlightLabeller._landmarks.size();
  highlightRenderer.renderMe(gSolver, highlightLabeller, params);
  TEST(highlightLabeller._landmarks.size() == landmarksBefore + 1);

  Config fallbackCfg = badgeCfg;
  fallbackCfg.renderMe = true;
  fallbackCfg.renderMeLabel = true;
  fallbackCfg.highlightMeStationLabel = false;
  fallbackCfg.meLabelSizeExplicit = true;
  fallbackCfg.meLabelSize = 24.0;
  fallbackCfg.meLandmark.coord = DPoint(5.0, 5.0);
  fallbackCfg.meLandmark.label = "Here";
  fallbackCfg.meLandmark.fontSize = fallbackCfg.meLabelSize;
  fallbackCfg.meLandmark.size = fallbackCfg.meStarSize;

  Labeller fallbackLabeller(&fallbackCfg);
  std::ostringstream fallbackOut;
  SvgRenderer fallbackRenderer(&fallbackOut, &fallbackCfg);
  RenderParams fallbackParams;
  fallbackParams.xOff = 0.0;
  fallbackParams.yOff = 0.0;
  fallbackParams.width = 400.0;
  fallbackParams.height = 400.0;

  size_t fallbackBefore = fallbackLabeller._landmarks.size();
  fallbackRenderer.renderMe(RenderGraph(), fallbackLabeller, fallbackParams);
  TEST(fallbackLabeller._landmarks.size() == fallbackBefore + 1);
  const auto &fallbackBox = fallbackLabeller._landmarks.back();
  double labelHeightPx = fallbackCfg.meLandmark.fontSize;
  size_t cpCount = util::toWStr(fallbackCfg.meLandmark.label).size();
  double labelWidthPx = cpCount * (labelHeightPx * 0.6);
  double starPx = fallbackCfg.meStarSize * fallbackCfg.outputResolution;
  double starGapPx = starPx * 0.2;
  double textHeightForPadding = labelHeightPx;
  double padX = textHeightForPadding * 0.6;
  double padTop = textHeightForPadding * 0.28;
  double padBottom = textHeightForPadding * 0.12;
  double contentHeightPx = std::max(starPx, textHeightForPadding);
  double expectedWidthPx =
      padX * 2.0 + starPx + starGapPx + labelWidthPx;
  double expectedHeightPx = padTop + padBottom + contentHeightPx;
  double actualWidth = fallbackBox.getUpperRight().getX() -
                       fallbackBox.getLowerLeft().getX();
  double actualHeight = fallbackBox.getUpperRight().getY() -
                        fallbackBox.getLowerLeft().getY();
  TEST(std::abs(actualWidth - expectedWidthPx) < 1e-6);
  TEST(std::abs(actualHeight - expectedHeightPx) < 1e-6);
}

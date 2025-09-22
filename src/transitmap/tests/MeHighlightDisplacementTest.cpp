#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#define private public
#include "transitmap/output/SvgRenderer.h"
#undef private

#include "shared/linegraph/Line.h"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/label/Labeller.h"
#include "transitmap/tests/MeHighlightDisplacementTest.h"
#include "util/Misc.h"
#include "util/geo/PolyLine.h"

using shared::linegraph::Line;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNodePL;
using shared::linegraph::Station;
using shared::rendergraph::RenderGraph;
using transitmapper::config::Config;
using transitmapper::label::Labeller;
using transitmapper::label::Overlaps;
using transitmapper::label::StationLabel;
using transitmapper::output::RenderParams;
using transitmapper::output::SvgRenderer;
using util::geo::DPoint;
using util::geo::MultiLine;
using util::geo::PolyLine;

void MeHighlightDisplacementTest::run() {
  Config cfg;
  cfg.outputResolution = 1.0;
  cfg.renderMe = true;
  cfg.renderMeLabel = true;
  cfg.highlightMeStationLabel = true;
  cfg.meStationWithBg = true;
  cfg.meStation = "here";
  cfg.meStarSize = 20.0;
  cfg.meStationBgFill = "#112233";
  cfg.meStationBgStroke = "#445566";
  cfg.meStationFill = "#778899";
  cfg.meStationBorder = "#aabbcc";
  cfg.meStationTextColor = "#ddeeff";
  cfg.displacementIterations = 50;
  cfg.displacementCooling = 0.8;
  cfg.landmarkSearchRadius = 4.0;
  cfg.lineWidth = 4.0;
  cfg.lineSpacing = 2.0;

  RenderGraph g;
  std::vector<std::unique_ptr<Line>> lines;
  lines.emplace_back(new Line("L1", "L1", "#000"));
  auto* center = g.addNd(LineNodePL(DPoint(0.0, 0.0)));
  center->pl().addStop(Station("center", "center", *center->pl().getGeom()));
  auto* east = g.addNd(LineNodePL(DPoint(20.0, 0.0)));
  auto* west = g.addNd(LineNodePL(DPoint(-20.0, 0.0)));
  PolyLine<double> eastGeom;
  eastGeom << *center->pl().getGeom() << *east->pl().getGeom();
  auto* eastEdge = g.addEdg(center, east, LineEdgePL(eastGeom));
  eastEdge->pl().addLine(lines.back().get(), east);
  PolyLine<double> westGeom;
  westGeom << *west->pl().getGeom() << *center->pl().getGeom();
  auto* westEdge = g.addEdg(west, center, LineEdgePL(westGeom));
  westEdge->pl().addLine(lines.back().get(), center);

  lines.emplace_back(new Line("L2", "L2", "#000"));
  auto* north = g.addNd(LineNodePL(DPoint(0.0, 20.0)));
  PolyLine<double> northGeom;
  northGeom << *center->pl().getGeom() << *north->pl().getGeom();
  auto* northEdge = g.addEdg(center, north, LineEdgePL(northGeom));
  northEdge->pl().addLine(lines.back().get(), north);

  PolyLine<double> geom;
  geom << DPoint(0.0, 0.0) << DPoint(60.0, 0.0);

  MultiLine<double> band;
  util::geo::Line<double> bandLine;
  bandLine.push_back(DPoint(-2.0, -2.0));
  bandLine.push_back(DPoint(62.0, -2.0));
  band.push_back(bandLine);
  band.push_back(bandLine);

  Overlaps overlaps{0, 0, 0, 0, 0};
  Station station("s1", "Here", DPoint(0.0, 0.0));
  StationLabel label(geom, band, 20.0, false, 0, 0, overlaps, 0.0, 15.0, 0.0,
                     0.0, false, 1.0, 0.0, nullptr, station);

  SvgRenderer::StationLabelVisual highlight;
  highlight.label = &label;
  highlight.pathId = "path-1";
  highlight.shift = "0";
  highlight.textAnchor = "start";
  highlight.startOffset = "0";
  highlight.fontSizePx = 18.0;
  highlight.bold = false;

  std::ostringstream out;
  SvgRenderer renderer(&out, &cfg);
  renderer._meStationLabelVisual = highlight;

  Labeller labeller(&cfg);
  RenderParams params;
  params.xOff = 0;
  params.yOff = 0;
  params.width = 400;
  params.height = 400;

  renderer.renderMe(g, labeller, params);

  std::string svg = out.str();
  size_t rectPos = svg.find("<rect");
  TEST(rectPos != std::string::npos);
  size_t groupStart = svg.rfind("<g", rectPos);
  TEST(groupStart != std::string::npos);
  size_t groupEnd = svg.find('>', groupStart);
  TEST(groupEnd != std::string::npos);
  std::string groupTag = svg.substr(groupStart, groupEnd - groupStart + 1);
  size_t translatePos = groupTag.find("translate(");
  TEST(translatePos != std::string::npos);
  translatePos += std::string("translate(").size();
  size_t translateEnd = groupTag.find(')', translatePos);
  TEST(translateEnd != std::string::npos);
  std::string translate = groupTag.substr(translatePos, translateEnd - translatePos);
  size_t spacePos = translate.find(' ');
  TEST(spacePos != std::string::npos);
  double translateX = std::stod(translate.substr(0, spacePos));
  double translateY = std::stod(translate.substr(spacePos + 1));

  double baseX = (0.0 - params.xOff) * cfg.outputResolution;
  double baseY = params.height - (0.0 - params.yOff) * cfg.outputResolution;

  TEST(std::abs(translateX - baseX) > 1e-3 || std::abs(translateY - baseY) > 1e-3);
}

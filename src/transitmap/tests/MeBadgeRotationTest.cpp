#include <cmath>
#include <sstream>
#include <string>

#define private public
#include "transitmap/output/SvgRenderer.h"
#undef private

#include "shared/linegraph/LineNodePL.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/label/Labeller.h"
#include "transitmap/tests/MeBadgeRotationTest.h"
#include "util/Misc.h"
#include "util/geo/Line.h"
#include "util/geo/PolyLine.h"

using shared::linegraph::Station;
using shared::rendergraph::RenderGraph;
using transitmapper::config::Config;
using transitmapper::label::Labeller;
using transitmapper::label::Overlaps;
using transitmapper::label::StationLabel;
using transitmapper::output::RenderParams;
using transitmapper::output::SvgRenderer;
using util::geo::DPoint;
using util::geo::Line;
using util::geo::MultiLine;
using util::geo::PolyLine;

void MeBadgeRotationTest::run() {
  Config cfg;
  cfg.outputResolution = 1.0;
  cfg.renderMe = true;
  cfg.renderMeLabel = true;
  cfg.highlightMeStationLabel = true;
  cfg.meStationWithBg = true;
  cfg.meStation = "here";
  cfg.meStationId = "here";
  cfg.meStationLabel = "Here";
  cfg.meStarSize = 20.0;
  cfg.meStationBgFill = "#112233";
  cfg.meStationBgStroke = "#445566";
  cfg.meStationFill = "#778899";
  cfg.meStationBorder = "#aabbcc";
  cfg.meStationTextColor = "#ddeeff";

  std::ostringstream out;
  SvgRenderer renderer(&out, &cfg);

  PolyLine<double> geom;
  geom << DPoint(20.0, 20.0) << DPoint(80.0, 50.0);

  Line<double> bandLine;
  bandLine.push_back(DPoint(18.0, 18.0));
  bandLine.push_back(DPoint(82.0, 52.0));
  MultiLine<double> band;
  band.push_back(bandLine);

  Overlaps overlaps{0, 0, 0, 0, 0};
  Station station("s1", "Here", DPoint(20.0, 20.0));
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
  renderer._meStationLabelVisual = highlight;

  RenderGraph g;
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
  size_t polygonPos = svg.find("<polygon", rectPos);
  TEST(polygonPos != std::string::npos);

  size_t groupStart = svg.rfind("<g", rectPos);
  TEST(groupStart != std::string::npos);
  size_t groupEnd = svg.find('>', groupStart);
  TEST(groupEnd != std::string::npos);
  std::string groupTag = svg.substr(groupStart, groupEnd - groupStart + 1);
  size_t transformPos = groupTag.find("transform=\"");
  TEST(transformPos != std::string::npos);
  transformPos += std::string("transform=\"").size();
  size_t transformEnd = groupTag.find('"', transformPos);
  TEST(transformEnd != std::string::npos);
  std::string transform = groupTag.substr(transformPos, transformEnd - transformPos);
  size_t rotatePos = transform.find("rotate(");
  TEST(rotatePos != std::string::npos);
  rotatePos += std::string("rotate(").size();
  size_t rotateEnd = transform.find(')', rotatePos);
  TEST(rotateEnd != std::string::npos);
  double angle = std::stod(transform.substr(rotatePos, rotateEnd - rotatePos));
  TEST(std::abs(angle) > 1.0);
}

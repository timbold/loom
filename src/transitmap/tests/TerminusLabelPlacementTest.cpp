#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#define private public
#include "transitmap/label/Labeller.h"
#include "transitmap/output/SvgRenderer.h"
#undef private

#include "shared/linegraph/Line.h"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/tests/TerminusLabelPlacementTest.h"
#include "util/Misc.h"
#include "util/geo/PolyLine.h"

using shared::linegraph::Line;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNode;
using shared::linegraph::LineNodePL;
using shared::linegraph::NodeFront;
using shared::rendergraph::RenderGraph;
using transitmapper::config::Config;
using transitmapper::config::TerminusLabelAnchor;
using transitmapper::label::Labeller;
using transitmapper::label::Overlaps;
using transitmapper::label::StationLabel;
using transitmapper::output::RenderParams;
using transitmapper::output::SvgRenderer;
using util::geo::DPoint;
using util::geo::PolyLine;
using util::approx;

namespace {

struct UniformBoxMetrics {
  double fontSize;
  double padTop;
  double padBottom;
  double padX;
  double charW;
  double boxH;
  double boxGap;

  double uniformBoxWidth() const { return 5 * charW + padX * 2; }
  double shiftDistance() const { return uniformBoxWidth() + boxGap; }
};

double parseRectAttribute(const std::string &svg, size_t rectIndex,
                          const std::string &attr) {
  size_t pos = svg.find("<rect");
  for (size_t i = 0; i < rectIndex && pos != std::string::npos; ++i) {
    pos = svg.find("<rect", pos + 1);
  }
  if (pos == std::string::npos)
    return std::numeric_limits<double>::quiet_NaN();

  std::string key = attr + "=\"";
  size_t attrPos = svg.find(key, pos);
  if (attrPos == std::string::npos)
    return std::numeric_limits<double>::quiet_NaN();
  attrPos += key.size();
  size_t end = svg.find('"', attrPos);
  if (end == std::string::npos)
    return std::numeric_limits<double>::quiet_NaN();
  return atof(svg.substr(attrPos, end - attrPos).c_str());
}

util::geo::Line<double> makeLine(std::initializer_list<DPoint> pts) {
  util::geo::Line<double> line;
  for (const auto &p : pts) {
    line.push_back(p);
  }
  return line;
}

} // namespace

void TerminusLabelPlacementTest::run() {
  auto computeUniformBoxMetrics = [](const Config &cfg) {
    double fontSize = cfg.lineLabelSize * cfg.outputResolution;
    double padTop = fontSize * 0.28;
    double padBottom = fontSize * 0.12;
    double padX = fontSize * 0.2;
    double charW = fontSize * 0.6;
    double boxH = fontSize + padTop + padBottom;
    double boxGap = cfg.routeLabelBoxGap * cfg.outputResolution;
    return UniformBoxMetrics{fontSize, padTop, padBottom, padX, charW, boxH,
                             boxGap};
  };

  auto runLabelCollisionScenario = [&](double labelLeft, double labelRight,
                                       int expectedShiftSign) {
    Config cfg;
    cfg.outputResolution = 1.0;
    cfg.lineLabelSize = 10.0;
    cfg.routeLabelBoxGap = 4.0;
    cfg.routeLabelTerminusGap = 0.0;
    cfg.stationLabelSize = 60.0;
    cfg.terminusLabelMaxLateralShift = 2;

    RenderGraph g;
    LineNode *term = g.addNd(LineNodePL(DPoint(0.0, 0.0)));
    LineNode *other = g.addNd(LineNodePL(DPoint(0.0, -100.0)));

    PolyLine<double> edgeGeom;
    edgeGeom << *term->pl().getGeom() << *other->pl().getGeom();
    auto *edge = g.addEdg(term, other, LineEdgePL(edgeGeom));

    Line line("L1", "L1", "f00");
    edge->pl().addLine(&line, other);

    shared::linegraph::Station stop("S1", "S1", *term->pl().getGeom());
    term->pl().addStop(stop);

    PolyLine<double> labelGeom;
    labelGeom << DPoint(labelLeft, 10.0) << DPoint(labelRight, 10.0);
    util::geo::MultiLine<double> band;
    band.push_back(makeLine({DPoint(labelLeft, 10.0), DPoint(labelRight, 10.0)}));
    band.push_back(makeLine({DPoint(labelRight, 10.0), DPoint(labelRight, 14.0)}));
    band.push_back(makeLine({DPoint(labelRight, 14.0), DPoint(labelLeft, 14.0)}));
    band.push_back(makeLine({DPoint(labelLeft, 14.0), DPoint(labelLeft, 10.0)}));

    Overlaps overlaps{0, 0, 0, 0, 0};
    Labeller labeller(&cfg);
    labeller._stationLabels.clear();
    std::vector<const shared::linegraph::Line*> lines;
    labeller._stationLabels.emplace_back(
        labelGeom, band, cfg.stationLabelSize, false, lines, 0, 0, overlaps, 0.0,
        cfg.stationLineOverlapPenalty, 0.0, 0.0, false, 1.0, 0.0, nullptr,
        stop);

    RenderParams params{400.0, 400.0, 0, 0};
    std::ostringstream out;
    SvgRenderer renderer(&out, &cfg);
    renderer.renderTerminusLabels(g, labeller, params);

    std::string svg = out.str();
    double rectX = parseRectAttribute(svg, 0, "x");
    double rectY = parseRectAttribute(svg, 0, "y");
    TEST(!std::isnan(rectX));
    TEST(!std::isnan(rectY));

    UniformBoxMetrics metrics = computeUniformBoxMetrics(cfg);
    double uniformBoxW = metrics.uniformBoxWidth();
    double totalW = uniformBoxW;
    double labelCenterX = (labelLeft + labelRight) / 2.0;
    double baseStartX = (labelCenterX - params.xOff) * cfg.outputResolution -
                        totalW / 2.0;
    double shiftDistance = metrics.shiftDistance();
    double expectedRectX =
        baseStartX + static_cast<double>(expectedShiftSign) * shiftDistance;
    TEST(rectX, ==, approx(expectedRectX));

    double labelCenterY = (10.0 + 14.0) / 2.0;
    double labelVExtent = 2.0;
    double anchorY = labelCenterY + labelVExtent;
    double yPx = params.height - (anchorY - params.yOff) * cfg.outputResolution;
    double stationHalfHeight =
        std::abs(labelVExtent * cfg.outputResolution);
    double stackCenterOffset =
        stationHalfHeight + cfg.routeLabelTerminusGap * cfg.outputResolution;
    double stackCenterY = yPx - stackCenterOffset;
    double expectedRectY = stackCenterY - metrics.boxH / 2.0;
    TEST(rectY, ==, approx(expectedRectY));
  };

  auto runFootprintCollisionScenario = [&]() {
    Config cfg;
    cfg.outputResolution = 1.0;
    cfg.lineLabelSize = 10.0;
    cfg.routeLabelBoxGap = 4.0;
    cfg.routeLabelTerminusGap = 0.0;
    cfg.terminusLabelAnchor = TerminusLabelAnchor::Node;
    cfg.terminusLabelMaxLateralShift = 2;

    RenderGraph g;
    LineNode *term = g.addNd(LineNodePL(DPoint(0.0, 0.0)));
    LineNode *other = g.addNd(LineNodePL(DPoint(10.0, 0.0)));

    PolyLine<double> edgeGeom;
    edgeGeom << *term->pl().getGeom() << *other->pl().getGeom();
    auto *edge = g.addEdg(term, other, LineEdgePL(edgeGeom));

    Line line("F1", "F1", "0f0");
    edge->pl().addLine(&line, other);

    NodeFront front(term, edge);
    PolyLine<double> frontGeom;
    frontGeom << DPoint(-10.0, 0.0) << DPoint(10.0, 0.0);
    front.setInitialGeom(frontGeom);
    term->pl().addFront(front);

    Labeller labeller(&cfg);
    RenderParams params{400.0, 400.0, 0, 0};
    std::ostringstream out;
    SvgRenderer renderer(&out, &cfg);
    renderer.renderTerminusLabels(g, labeller, params);

    std::string svg = out.str();
    double rectX = parseRectAttribute(svg, 0, "x");
    double rectY = parseRectAttribute(svg, 0, "y");
    TEST(!std::isnan(rectX));
    TEST(!std::isnan(rectY));

    UniformBoxMetrics metrics = computeUniformBoxMetrics(cfg);
    double uniformBoxW = metrics.uniformBoxWidth();
    double totalW = uniformBoxW;
    double baseStartX = -totalW / 2.0;
    double shiftDistance = metrics.shiftDistance();
    TEST(rectX, ==, approx(baseStartX - shiftDistance));

    double yPx = params.height - (0.0 - params.yOff) * cfg.outputResolution;
    double expectedRectY = yPx - metrics.boxH / 2.0;
    TEST(rectY, ==, approx(expectedRectY));
  };

  auto runStackAvoidanceScenario = [&]() {
    Config cfg;
    cfg.outputResolution = 1.0;
    cfg.lineLabelSize = 10.0;
    cfg.routeLabelBoxGap = 4.0;
    cfg.routeLabelTerminusGap = 0.0;
    cfg.terminusLabelAnchor = TerminusLabelAnchor::Node;
    cfg.terminusLabelMaxLateralShift = 2;

    RenderGraph g;
    LineNode *aTerm = g.addNd(LineNodePL(DPoint(0.0, 0.0)));
    LineNode *aOther = g.addNd(LineNodePL(DPoint(0.0, -10.0)));
    LineNode *bTerm = g.addNd(LineNodePL(DPoint(0.1, 0.0)));
    LineNode *bOther = g.addNd(LineNodePL(DPoint(0.1, -10.0)));

    PolyLine<double> aGeom;
    aGeom << *aTerm->pl().getGeom() << *aOther->pl().getGeom();
    auto *aEdge = g.addEdg(aTerm, aOther, LineEdgePL(aGeom));
    PolyLine<double> bGeom;
    bGeom << *bTerm->pl().getGeom() << *bOther->pl().getGeom();
    auto *bEdge = g.addEdg(bTerm, bOther, LineEdgePL(bGeom));

    Line lineA("A", "A", "00f");
    Line lineB("B", "B", "f0f");
    aEdge->pl().addLine(&lineA, aOther);
    bEdge->pl().addLine(&lineB, bOther);

    Labeller labeller(&cfg);
    RenderParams params{400.0, 400.0, 0, 0};
    std::ostringstream out;
    SvgRenderer renderer(&out, &cfg);
    renderer.renderTerminusLabels(g, labeller, params);

    std::string svg = out.str();
    double rectX0 = parseRectAttribute(svg, 0, "x");
    double rectX2 = parseRectAttribute(svg, 2, "x");
    TEST(!std::isnan(rectX0));
    TEST(!std::isnan(rectX2));

    UniformBoxMetrics metrics = computeUniformBoxMetrics(cfg);
    double uniformBoxW = metrics.uniformBoxWidth();
    double totalW = uniformBoxW;
    double baseStartA =
        (aTerm->pl().getGeom()->getX() - params.xOff) * cfg.outputResolution -
        totalW / 2.0;
    double baseStartB =
        (bTerm->pl().getGeom()->getX() - params.xOff) * cfg.outputResolution -
        totalW / 2.0;
    double shiftDistance = metrics.shiftDistance();

    TEST(rectX0, ==, approx(baseStartA));
    TEST(rectX2, ==, approx(baseStartB - shiftDistance));
  };

  runLabelCollisionScenario(-15.0, -5.0, -1);
  runLabelCollisionScenario(5.0, 15.0, 1);
  runFootprintCollisionScenario();
  runStackAvoidanceScenario();
}


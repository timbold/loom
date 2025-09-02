#include <sstream>

#define private public
#include "transitmap/output/SvgRenderer.h"
#undef private

#include "transitmap/tests/DirMarkerTest.h"
#include "util/Misc.h"

using transitmapper::config::Config;
using transitmapper::output::RenderParams;
using transitmapper::output::SvgRenderer;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNode;
using shared::linegraph::LineNodePL;
using shared::linegraph::NodeFront;
using shared::rendergraph::RenderGraph;
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

void DirMarkerTest::run() {
  Config cfg;
  cfg.renderDirMarkers = true;
  cfg.lineWidth = 20;
  cfg.outlineWidth = 1;
  cfg.lineSpacing = 10;
  cfg.outputResolution = 1.0;

  std::ostringstream out;
  SvgRenderer renderer(&out, &cfg);
  RenderParams params{};

  Line line("L1", "L1", "#000");
  RenderGraph g;

  for (int i = 0; i < 10; ++i) {
    LineEdge* e = makeEdge(g, &line, 10.0);
    renderer.renderEdgeTripGeom(g, e, params);
  }

  TEST(renderer._edgesSinceMarker[&line], ==, 1);

  LineEdge* longE = makeEdge(g, &line, 200.0);
  renderer.renderEdgeTripGeom(g, longE, params);

  TEST(renderer._edgesSinceMarker[&line], ==, 0);
}


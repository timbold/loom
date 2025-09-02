#include <cmath>
#include <sstream>

#define private public
#include "transitmap/output/SvgRenderer.h"
#undef private

#include "transitmap/tests/ArrowHeadDirectionTest.h"
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
LineEdge* makeEdge(RenderGraph& g, const Line* line, const DPoint& a,
                   const DPoint& b) {
  LineNode* na = g.addNd(LineNodePL(a));
  LineNode* nb = g.addNd(LineNodePL(b));

  PolyLine<double> pl;
  pl << a << b;
  LineEdgePL epl(pl);
  LineEdge* e = g.addEdg(na, nb, epl);
  e->pl().addLine(line, nb);

  NodeFront nfA(na, e);
  PolyLine<double> nfAP;
  nfAP << a << DPoint(a.getX(), a.getY() + 1);
  nfA.setInitialGeom(nfAP);
  na->pl().addFront(nfA);

  NodeFront nfB(nb, e);
  PolyLine<double> nfBP;
  nfBP << b << DPoint(b.getX(), b.getY() + 1);
  nfB.setInitialGeom(nfBP);
  nb->pl().addFront(nfB);

  return e;
}
}  // namespace

void ArrowHeadDirectionTest::run() {
  Config cfg;
  cfg.renderDirMarkers = true;
  cfg.renderBiDirMarker = true;
  cfg.lineWidth = 20;
  cfg.outlineWidth = 1;
  cfg.lineSpacing = 10;
  cfg.outputResolution = 1.0;

  std::ostringstream out;
  SvgRenderer renderer(&out, &cfg);
  RenderParams params{};

  Line line("L1", "L1", "#000");
  RenderGraph g;

  LineEdge* e1 = makeEdge(g, &line, DPoint(0, 0), DPoint(10, 0));
  LineEdge* e2 = makeEdge(g, &line, DPoint(10, 0), DPoint(20, 0));

  renderer.renderEdgeTripGeom(g, e1, params);
  renderer.renderEdgeTripGeom(g, e2, params);

  std::vector<const SvgRenderer::ArrowHead*> nearShared;
  for (const auto& ah : renderer._arrowHeads) {
    double baseX = (ah.pts[0].getX() + ah.pts[1].getX()) / 2.0;
    if (std::abs(baseX - 10.0) < 1e-6) {
      nearShared.push_back(&ah);
    }
  }

  TEST(nearShared.size(), ==, 2);

  auto isForward = [](const SvgRenderer::ArrowHead* ah) {
    double baseX = (ah->pts[0].getX() + ah->pts[1].getX()) / 2.0;
    double tipX = ah->pts[3].getX();
    return tipX > baseX;
  };

  bool dir0 = isForward(nearShared[0]);
  bool dir1 = isForward(nearShared[1]);

  TEST(dir0, ==, dir1);
}


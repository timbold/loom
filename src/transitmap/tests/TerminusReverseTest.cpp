#include "transitmap/tests/TerminusReverseTest.h"

#include "shared/rendergraph/RenderGraph.h"
#include "dot/Parser.h"
#include "shared/linegraph/Line.h"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineNodePL.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"

using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNode;
using shared::linegraph::LineNodePL;
using shared::linegraph::NodeFront;
using shared::rendergraph::RenderGraph;
using util::geo::DPoint;
using util::geo::PolyLine;

namespace dot {
namespace parser {
[[gnu::weak]] Parser::Parser(std::istream* stream) { UNUSED(stream); }

[[gnu::weak]] const Entity& Parser::get() {
  static Entity entity;
  entity.type = EMPTY;
  entity.ids.clear();
  entity.attrs.clear();
  entity.graphType = GRAPH;
  entity.graphName.clear();
  entity.level = 0;
  return entity;
}

[[gnu::weak]] bool Parser::has() { return false; }
}  // namespace parser
}  // namespace dot

namespace {
void addFront(LineNode* node, LineEdge* edge) {
  NodeFront front(node, edge);
  PolyLine<double> geom;
  const auto& pos = *node->pl().getGeom();
  geom << DPoint(pos.getX(), pos.getY());
  geom << DPoint(pos.getX(), pos.getY() + 1);
  front.setInitialGeom(geom);
  node->pl().addFront(front);
}

LineEdge* addEdge(RenderGraph& graph, LineNode* from, LineNode* to) {
  PolyLine<double> geom;
  const auto& fromPos = *from->pl().getGeom();
  const auto& toPos = *to->pl().getGeom();
  geom << DPoint(fromPos.getX(), fromPos.getY());
  geom << DPoint(toPos.getX(), toPos.getY());
  LineEdgePL pl(geom);
  LineEdge* edge = graph.addEdg(from, to, pl);
  addFront(from, edge);
  addFront(to, edge);
  return edge;
}
}  // namespace

void TerminusReverseTest::run() {
  RenderGraph graph;
  Line line("L1", "L1", "#000");

  LineNode* a = graph.addNd(LineNodePL(DPoint(0, 0)));
  LineNode* b = graph.addNd(LineNodePL(DPoint(1, 0)));
  LineNode* c = graph.addNd(LineNodePL(DPoint(2, 0)));
  LineNode* d = graph.addNd(LineNodePL(DPoint(1, -1)));

  LineEdge* ab = addEdge(graph, a, b);
  LineEdge* bc = addEdge(graph, b, c);
  LineEdge* bd = addEdge(graph, b, d);

  ab->pl().addLine(&line, b);
  bc->pl().addLine(&line, c);
  bc->pl().addLine(&line, b);
  bd->pl().addLine(&line, d);

  TEST(RenderGraph::isTerminus(a), ==, true);
  TEST(RenderGraph::isTerminus(b), ==, false);
  TEST(RenderGraph::isTerminus(c), ==, false);
  TEST(RenderGraph::isTerminus(d), ==, true);

  TEST(shared::linegraph::LineGraph::terminatesAt(ab, a, &line), ==, true);
  TEST(shared::linegraph::LineGraph::terminatesAt(bc, c, &line), ==, false);
  TEST(shared::linegraph::LineGraph::terminatesAt(bd, d, &line), ==, true);
}

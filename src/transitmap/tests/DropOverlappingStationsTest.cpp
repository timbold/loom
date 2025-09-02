#include "transitmap/tests/DropOverlappingStationsTest.h"

#include "transitmap/graph/GraphBuilder.h"
#include "transitmap/config/TransitMapConfig.h"
#include "util/Misc.h"

using transitmapper::config::Config;
using transitmapper::graph::GraphBuilder;
using shared::rendergraph::RenderGraph;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNode;
using shared::linegraph::LineNodePL;
using util::geo::DPoint;
using util::geo::PolyLine;

void DropOverlappingStationsTest::run() {
  Config cfg;
  GraphBuilder builder(&cfg);
  RenderGraph g;

  Line line("L1", "L1", "#000");

  LineNode* nonTerm = g.addNd(LineNodePL(DPoint(0, 0)));
  LineNode* term = g.addNd(LineNodePL(DPoint(0.5, 0)));
  LineNode* other = g.addNd(LineNodePL(DPoint(10, 0)));

  PolyLine<double> plAB;
  plAB << DPoint(0, 0) << DPoint(0.5, 0);
  LineEdgePL eplAB(plAB);
  LineEdge* eAB = g.addEdg(nonTerm, term, eplAB);
  eAB->pl().addLine(&line, term);

  PolyLine<double> plAC;
  plAC << DPoint(0, 0) << DPoint(10, 0);
  LineEdgePL eplAC(plAC);
  LineEdge* eAC = g.addEdg(nonTerm, other, eplAC);
  eAC->pl().addLine(&line, other);

  builder.writeNodeFronts(&g);

  nonTerm->pl().addStop(shared::linegraph::Station("A", "A", *nonTerm->pl().getGeom()));
  term->pl().addStop(shared::linegraph::Station("B", "B", *term->pl().getGeom()));

  TEST(nonTerm->pl().stops().size(), ==, 1);
  TEST(term->pl().stops().size(), ==, 1);

  builder.dropOverlappingStations(&g);

  TEST(nonTerm->pl().stops().size(), ==, 0);
  TEST(term->pl().stops().size(), ==, 1);
}


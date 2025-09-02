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

  Line wide1("W1", "W1", "#000");
  Line wide2("W2", "W2", "#000");
  Line wide3("W3", "W3", "#000");
  Line narrow("N1", "N1", "#000");

  LineNode* wide = g.addNd(LineNodePL(DPoint(0, 0)));
  LineNode* termWide1 = g.addNd(LineNodePL(DPoint(0, 10)));
  LineNode* termWide2 = g.addNd(LineNodePL(DPoint(-10, 0)));
  LineNode* narrowTerm = g.addNd(LineNodePL(DPoint(20, 0)));
  LineNode* narrowOther = g.addNd(LineNodePL(DPoint(20, 10)));

  PolyLine<double> plWide1;
  plWide1 << DPoint(0, 0) << DPoint(0, 10);
  LineEdge* eWide1 = g.addEdg(wide, termWide1, LineEdgePL(plWide1));
  eWide1->pl().addLine(&wide1, termWide1);
  eWide1->pl().addLine(&wide2, termWide1);
  eWide1->pl().addLine(&wide3, termWide1);

  PolyLine<double> plWide2;
  plWide2 << DPoint(0, 0) << DPoint(-10, 0);
  LineEdge* eWide2 = g.addEdg(wide, termWide2, LineEdgePL(plWide2));
  eWide2->pl().addLine(&wide1, termWide2);
  eWide2->pl().addLine(&wide2, termWide2);
  eWide2->pl().addLine(&wide3, termWide2);

  PolyLine<double> plNarrow;
  plNarrow << DPoint(20, 0) << DPoint(20, 10);
  LineEdge* eNarrow = g.addEdg(narrowTerm, narrowOther, LineEdgePL(plNarrow));
  eNarrow->pl().addLine(&narrow, narrowOther);

  builder.writeNodeFronts(&g);

  wide->pl().addStop(shared::linegraph::Station("W", "W", *wide->pl().getGeom()));
  narrowTerm->pl().addStop(shared::linegraph::Station("N", "N", *narrowTerm->pl().getGeom()));

  double padWide = 0;
  for (auto e : wide->getAdjList())
    padWide = std::max(padWide, g.getTotalWidth(e) + g.getSpacing(e));
  double padNarrow = 0;
  for (auto e : narrowTerm->getAdjList())
    padNarrow = std::max(padNarrow, g.getTotalWidth(e) + g.getSpacing(e));

  TEST(padWide > padNarrow, ==, 1);
}

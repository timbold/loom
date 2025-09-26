#include <memory>
#include <string>
#include <vector>

#include "transitmap/tests/SingleRouteLabelTest.h"

#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/label/Labeller.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"

using shared::linegraph::Line;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNodePL;
using shared::rendergraph::RenderGraph;
using transitmapper::config::Config;
using transitmapper::label::Labeller;
using util::geo::DPoint;
using util::geo::PolyLine;

namespace {

void buildEdgeWithLines(RenderGraph* g,
                        std::vector<std::unique_ptr<Line>>* lines,
                        size_t lineCount) {
  lines->clear();

  auto* nodeA = g->addNd(LineNodePL(DPoint(0.0, 0.0)));
  auto* nodeB = g->addNd(LineNodePL(DPoint(200.0, 0.0)));

  PolyLine<double> geom;
  geom << *nodeA->pl().getGeom() << *nodeB->pl().getGeom();
  auto* edge = g->addEdg(nodeA, nodeB, LineEdgePL(geom));

  for (size_t i = 0; i < lineCount; ++i) {
    std::string id = "L" + std::to_string(i + 1);
    lines->emplace_back(new Line(id, id, "#000"));
    edge->pl().addLine(lines->back().get(), nodeB);
  }
}

}  // namespace

void SingleRouteLabelTest::run() {
  // Default behaviour labels single-route edges.
  {
    Config cfg;
    RenderGraph g;
    std::vector<std::unique_ptr<Line>> lines;
    buildEdgeWithLines(&g, &lines, 1);

    Labeller labeller(&cfg);
    labeller.label(g, false);
    TEST(!labeller.getLineLabels().empty());
  }

  // When disabled, single-route edges are not labelled.
  {
    Config cfg;
    cfg.renderSingleRouteLabels = false;
    RenderGraph g;
    std::vector<std::unique_ptr<Line>> lines;
    buildEdgeWithLines(&g, &lines, 1);

    Labeller labeller(&cfg);
    labeller.label(g, false);
    TEST(labeller.getLineLabels().empty());
  }

  // Multi-route edges remain eligible when the flag is disabled.
  {
    Config cfg;
    cfg.renderSingleRouteLabels = false;
    RenderGraph g;
    std::vector<std::unique_ptr<Line>> lines;
    buildEdgeWithLines(&g, &lines, 2);

    Labeller labeller(&cfg);
    labeller.label(g, false);
    TEST(!labeller.getLineLabels().empty());
    TEST(labeller.getLineLabels()[0].lines.size(), ==, 2);
  }
}

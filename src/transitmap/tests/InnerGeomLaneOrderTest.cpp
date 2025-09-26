#include <algorithm>
#include <cstdint>
#include <sstream>
#include <utility>
#include <vector>

#define private public
#include "transitmap/output/SvgRenderer.h"
#undef private

#include "transitmap/tests/InnerGeomLaneOrderTest.h"
#include "util/Misc.h"

using transitmapper::config::Config;
using transitmapper::output::InnerClique;
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

void addFront(LineNode* node, LineEdge* edge, double dx, double dy) {
  NodeFront front(node, edge);
  PolyLine<double> geom;
  geom << DPoint(node->pl().getGeom()->getX(), node->pl().getGeom()->getY());
  geom << DPoint(node->pl().getGeom()->getX() + dx,
                 node->pl().getGeom()->getY() + dy);
  front.setInitialGeom(geom);
  node->pl().addFront(front);
}

LineEdge* addEdge(RenderGraph& graph, LineNode* from, LineNode* to,
                  const PolyLine<double>& geometry) {
  LineEdge* edge = graph.addEdg(from, to, LineEdgePL(geometry));
  addFront(from, edge, 0, 1);
  addFront(to, edge, 0, -1);
  return edge;
}

}  // namespace

void InnerGeomLaneOrderTest::run() {
  Config cfg;
  cfg.innerGeometryPrecision = 0.0;

  std::ostringstream sink;
  SvgRenderer renderer(&sink, &cfg);

  RenderGraph graph;
  Line lineA("A", "A", "#111111");
  Line lineB("B", "B", "#222222");

  LineNode* left = graph.addNd(LineNodePL(DPoint(0, 0)));
  LineNode* right = graph.addNd(LineNodePL(DPoint(10, 0)));

  PolyLine<double> upperGeom;
  upperGeom << DPoint(0, 0) << DPoint(10, 0);
  LineEdge* upper = addEdge(graph, left, right, upperGeom);

  PolyLine<double> lowerGeom;
  lowerGeom << DPoint(0, 0) << DPoint(5, -2) << DPoint(10, 0);
  LineEdge* lower = addEdge(graph, left, right, lowerGeom);

  upper->pl().addLine(&lineA, right);
  upper->pl().addLine(&lineB, right);
  lower->pl().addLine(&lineA, left);
  lower->pl().addLine(&lineB, left);

  auto geoms = graph.innerGeoms(left, cfg.innerGeometryPrecision);
  auto cliques = renderer.getInnerCliques(left, geoms, 9999);

  std::vector<std::pair<size_t, size_t>> slotPairs;
  const InnerClique* targetClique = nullptr;
  for (const auto& clique : cliques) {
    for (const auto& geom : clique.geoms) {
      if (geom.from.edge == upper && geom.to.edge == lower) {
        slotPairs.emplace_back(geom.slotFrom, geom.slotTo);
        if (!targetClique) {
          targetClique = &clique;
        }
      }
    }
  }

  TEST(slotPairs.size(), ==, 2);
  TEST(targetClique, !=, nullptr);

  renderer.renderClique(*targetClique, graph, left);

  TEST(renderer._innerDelegates.empty(), ==, false);
  const auto& delegateMap = renderer._innerDelegates.back();

  auto itA = delegateMap.find(reinterpret_cast<uintptr_t>(&lineA));
  auto itB = delegateMap.find(reinterpret_cast<uintptr_t>(&lineB));
  TEST(itA, !=, delegateMap.end());
  TEST(itB, !=, delegateMap.end());
  TEST(itA->second.empty(), ==, false);
  TEST(itB->second.empty(), ==, false);

  const auto& polyA = itA->second.front().front.second;
  const auto& polyB = itB->second.front().front.second;

  const auto& aFirst = polyA.front().getX() <= polyA.back().getX()
                           ? polyA.front()
                           : polyA.back();
  const auto& bFirst = polyB.front().getX() <= polyB.back().getX()
                           ? polyB.front()
                           : polyB.back();

  TEST(aFirst.getY(), >=, bFirst.getY());

  std::sort(slotPairs.begin(), slotPairs.end());
  TEST(slotPairs[0].first, ==, 0u);
  TEST(slotPairs[1].first, ==, 1u);
  TEST(slotPairs[0].second, ==, 0u);
  TEST(slotPairs[1].second, ==, 1u);

  renderer._innerDelegates.clear();
}

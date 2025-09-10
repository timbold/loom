#include <fstream>
#include <sstream>
#include <vector>

#define private public
#include "transitmap/output/SvgRenderer.h"
#undef private

#include "transitmap/tests/BgMapTest.h"
#include "transitmap/config/ConfigReader.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"

using transitmapper::config::Config;
using transitmapper::config::ConfigReader;
using transitmapper::output::SvgRenderer;
using shared::rendergraph::RenderGraph;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineNode;
using shared::linegraph::LineNodePL;
using shared::linegraph::NodeFront;
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

void BgMapTest::run() {
  std::string path = "bgmap_test.geojson";
  {
    std::ofstream out(path);
    out << "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[0,0],[10,0]]}}]}";
  }

  Config cfg;
  const char* argv[] = {"prog", "--bg-map", path.c_str()};
  ConfigReader reader;
  reader.read(&cfg, 3, const_cast<char**>(argv));

  RenderGraph g;
  Line line("L1", "L1", "#000");
  makeEdge(g, &line, 10.0);

  std::ostringstream svg;
  SvgRenderer s(&svg, &cfg);
  s.print(g);
  TEST(svg.str().find("bg-map") != std::string::npos, ==, true);

  // Provide lat/long coordinates and ensure they are converted to Web Mercator
  std::string path2 = "bgmap_test_latlng.geojson";
  {
    std::ofstream out(path2);
    out << "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[0,0],[0,1]]}}]}";
  }

  Config cfg2;
  const char* argv2[] = {"prog", "--bg-map", path2.c_str()};
  reader.read(&cfg2, 3, const_cast<char**>(argv2));

  std::ostringstream svg2;
  SvgRenderer s2(&svg2, &cfg2);
  s2.print(g);

  auto expected = util::geo::latLngToWebMerc(DPoint(0, 1)).getY() *
                  cfg2.outputResolution;

  std::string outStr = svg2.str();
  size_t pos = outStr.find("points=\"");
  TEST(pos != std::string::npos, ==, true);
  pos += 8;  // skip 'points="'
  size_t end = outStr.find("\"", pos);
  std::stringstream pts(outStr.substr(pos, end - pos));
  double x1, y1, x2, y2;
  char comma;
  pts >> x1 >> comma >> y1 >> x2 >> comma >> y2;
  TEST(x1, ==, 0);
  TEST(std::abs(y1 - expected) < 1e-6);
  TEST(x2, ==, 0);
  TEST(std::abs(y2) < 1e-6);
  // Verify that background map influences overall SVG dimensions.
  std::string path3 = "bgmap_bbox.geojson";
  {
    std::ofstream out(path3);
    out << "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[0,0],[100,100]]}}]}";
  }

  Config cfgBase;
  const char* argvBase[] = {"prog"};
  reader.read(&cfgBase, 1, const_cast<char**>(argvBase));
  std::ostringstream svgBase;
  SvgRenderer sBase(&svgBase, &cfgBase);
  sBase.print(g);
  std::string baseStr = svgBase.str();
  auto parseDim = [](const std::string& s, const std::string& attr) {
    size_t p = s.find(attr + "=\"");
    if (p == std::string::npos) return 0.0;
    p += attr.size() + 2;
    size_t q = s.find("\"", p);
    return std::stod(s.substr(p, q - p));
  };
  double baseW = parseDim(baseStr, "width");
  double baseH = parseDim(baseStr, "height");

  Config cfgNoExt;
  const char* argvNoExt[] = {"prog", "--bg-map", path3.c_str(),
                             "--bg-map-webmerc"};
  reader.read(&cfgNoExt, 4, const_cast<char**>(argvNoExt));
  std::ostringstream svgNoExt;
  SvgRenderer sNoExt(&svgNoExt, &cfgNoExt);
  sNoExt.print(g);
  std::string outNoExt = svgNoExt.str();
  double wNoExt = parseDim(outNoExt, "width");
  double hNoExt = parseDim(outNoExt, "height");
  TEST(std::abs(wNoExt - baseW) < 1e-6);
  TEST(std::abs(hNoExt - baseH) < 1e-6);

  Config cfg3;
  const char* argv3[] = {"prog",        "--bg-map", path3.c_str(),
                         "--bg-map-webmerc", "--extend-with-bgmap"};
  reader.read(&cfg3, 5, const_cast<char**>(argv3));
  std::ostringstream svg3;
  SvgRenderer s3(&svg3, &cfg3);
  s3.print(g);
  std::string out3 = svg3.str();
  double w = parseDim(out3, "width");
  double h = parseDim(out3, "height");
  auto bgBox = s3.computeBgMapBBox();
  auto netBox = util::geo::pad(g.getBBox(),
                               g.getMaxLineNum() *
                                   (cfg3.lineWidth + cfg3.lineSpacing));
  auto merged = util::geo::extendBox(bgBox, netBox);
  double expW =
      (merged.getUpperRight().getX() - merged.getLowerLeft().getX() +
       cfg3.paddingLeft + cfg3.paddingRight) * cfg3.outputResolution;
  double expH =
      (merged.getUpperRight().getY() - merged.getLowerLeft().getY() +
       cfg3.paddingTop + cfg3.paddingBottom) * cfg3.outputResolution;
  TEST(std::abs(w - expW) < 1e-6);
  TEST(std::abs(h - expH) < 1e-6);
  TEST(w > baseW);
  TEST(h > baseH);

}

#include <fstream>
#include <sstream>
#include <stdexcept>
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

  Config cfgOpacity;
  const char* argvOp[] = {"prog", "--bg-map", path.c_str(),
                          "--bg-map-opacity", "0.5"};
  reader.read(&cfgOpacity, 5, const_cast<char**>(argvOp));
  std::ostringstream svgOp;
  SvgRenderer sOp(&svgOp, &cfgOpacity);
  sOp.print(g);
  std::string svgOpStr = svgOp.str();
  TEST(svgOpStr.find("stroke-opacity:0.5") != std::string::npos, ==, true);

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

  // Render polygons with explicit style properties and ensure they are preserved.
  std::string stylePath = "bgmap_style.geojson";
  {
    std::ofstream out(stylePath);
    out << R"({"type":"FeatureCollection","features":[
      {"type":"Feature","properties":{"stroke":"#f00","stroke-width":2,"fill":"#0f0","opacity":0.5,"class":"my-poly"},"geometry":{"type":"Polygon","coordinates":[[[0,0],[0,1],[1,1],[1,0],[0,0]]]}},
      {"type":"Feature","properties":{"stroke":"#00f","stroke-width":1,"fill":"#fff","opacity":0.3},"geometry":{"type":"MultiPolygon","coordinates":[[[[2,0],[3,0],[3,1],[2,1],[2,0]]],[[[4,0],[5,0],[5,1],[4,1],[4,0]]]]}}
    ]})";
  }

  Config cfgStyle;
  const char *argvStyle[] = {"prog", "--bg-map", stylePath.c_str()};
  reader.read(&cfgStyle, 3, const_cast<char **>(argvStyle));
  cfgStyle.outputResolution = 1.0;
  RenderGraph gStyle;
  makeEdge(gStyle, &line, 1.0);
  std::ostringstream svgStyle;
  SvgRenderer sStyle(&svgStyle, &cfgStyle);
  sStyle.print(gStyle);
  std::string outStyle = svgStyle.str();
  TEST(outStyle.find("class=\"bg-map my-poly\"") != std::string::npos, ==, true);
  TEST(outStyle.find("stroke:#f00") != std::string::npos, ==, true);
  TEST(outStyle.find("stroke-width:2") != std::string::npos, ==, true);
  TEST(outStyle.find("fill:#0f0") != std::string::npos, ==, true);
  TEST(outStyle.find("stroke-opacity:0.5") != std::string::npos, ==, true);
  TEST(outStyle.find("fill-opacity:0.5") != std::string::npos, ==, true);
  int mpCount = 0;
  for (size_t pos = outStyle.find("stroke:#00f"); pos != std::string::npos;
       pos = outStyle.find("stroke:#00f", pos + 1)) {
    mpCount++;
  }
  TEST(mpCount, ==, 2);

  // Render polygons without properties and ensure default styling is applied.
  std::string noPropPath = "bgmap_noprop.geojson";
  {
    std::ofstream out(noPropPath);
    out << R"({"type":"FeatureCollection","features":[
      {"type":"Feature","geometry":{"type":"Polygon","coordinates":[[[0,0],[0,1],[1,1],[1,0],[0,0]]]}},
      {"type":"Feature","geometry":{"type":"MultiPolygon","coordinates":[[[[2,0],[3,0],[3,1],[2,1],[2,0]]]]}}
    ]})";
  }

  Config cfgDefault;
  const char *argvDefault[] = {"prog", "--bg-map", noPropPath.c_str()};
  reader.read(&cfgDefault, 3, const_cast<char **>(argvDefault));
  cfgDefault.outputResolution = 1.0;
  RenderGraph gDefault;
  makeEdge(gDefault, &line, 1.0);
  std::ostringstream svgDefault;
  SvgRenderer sDefault(&svgDefault, &cfgDefault);
  sDefault.print(gDefault);
  std::string outDefault = svgDefault.str();
  TEST(outDefault.find("stroke:#ccc") != std::string::npos, ==, true);
  TEST(outDefault.find("stroke-width:20") != std::string::npos, ==, true);
  TEST(outDefault.find("fill:none") != std::string::npos, ==, true);
  TEST(outDefault.find("stroke-opacity:1") != std::string::npos, ==, true);
  TEST(outDefault.find("fill-opacity:1") != std::string::npos, ==, true);
  int defaultCount = 0;
  for (size_t pos = outDefault.find("stroke:#ccc"); pos != std::string::npos;
       pos = outDefault.find("stroke:#ccc", pos + 1)) {
    defaultCount++;
  }
  TEST(defaultCount, ==, 2);

  // Render a feature with a very large coordinate list and ensure streaming
  // keeps memory bounded.
  std::string hugePath = "bgmap_huge_line.geojson";
  {
    std::ofstream out(hugePath);
    out << "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"geometry\":{\"type\":\"LineString\",\"coordinates\":[";
    const size_t numPoints = 200000;
    for (size_t i = 0; i < numPoints; ++i) {
      if (i)
        out << ",";
      out << "[" << (i % 1024) << "," << (i / 1024) << "]";
    }
    out << "]}}]}";
  }

  Config cfgHuge;
  const char *argvHuge[] = {"prog", "--bg-map", hugePath.c_str()};
  reader.read(&cfgHuge, 3, const_cast<char **>(argvHuge));
  std::ostringstream svgHuge;
  SvgRenderer sHuge(&svgHuge, &cfgHuge);
  bool lengthErrorThrown = false;
  try {
    sHuge.print(g);
  } catch (const std::length_error &) {
    lengthErrorThrown = true;
  }
  TEST(lengthErrorThrown, ==, false);

}

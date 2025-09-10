// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_SVGRENDERER_H_
#define TRANSITMAP_OUTPUT_SVGRENDERER_H_

#include <ostream>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "Renderer.h"
#include "shared/linegraph/Line.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/label/Labeller.h"
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"
#include "util/xml/XmlWriter.h"

using util::Nullable;

namespace transitmapper {
namespace output {

class SvgRenderer : public Renderer {
 public:
  SvgRenderer(std::ostream* o, const config::Config* cfg);
  virtual ~SvgRenderer(){};

  virtual void print(const shared::rendergraph::RenderGraph& outputGraph);

  void printLine(const util::geo::PolyLine<double>& l,
                 const std::map<std::string, std::string>& ps,
                 const RenderParams& params);
  void printLine(const util::geo::PolyLine<double>& l, const std::string& style,
                 const RenderParams& params);
  void printPoint(const util::geo::DPoint& p, const std::string& style,
                  const RenderParams& params);
  void printPolygon(const util::geo::Polygon<double>& g,
                    const std::map<std::string, std::string>& ps,
                    const RenderParams& params);
  void printCircle(const util::geo::DPoint& center, double rad,
                   const std::map<std::string, std::string>& ps,
                   const RenderParams& rparams);
  void printCircle(const util::geo::DPoint& center, double rad,
                   const std::string& style, const RenderParams& rparams);

 private:
  std::ostream* _o;
  util::xml::XmlWriter _w;

  const config::Config* _cfg;

  std::map<uintptr_t, std::vector<OutlinePrintPair>> _delegates;
  std::vector<std::map<uintptr_t, std::vector<OutlinePrintPair>>>
      _innerDelegates;
  struct ArrowHead {
    std::vector<util::geo::DPoint> pts;
  };
  std::vector<ArrowHead> _arrowHeads;
  mutable std::map<std::string, int> lineClassIds;
  mutable int lineClassId = 0;
  std::unordered_map<const shared::linegraph::Line*, int> _edgesSinceMarker;
  std::unordered_map<const shared::linegraph::Line*,
                     std::unordered_set<const shared::linegraph::LineEdge*>>
      _forceDirMarker;

  util::geo::Box<double> computeBgMapBBox() const;

  void outputNodes(const shared::rendergraph::RenderGraph& outputGraph,
                   const RenderParams& params);
  void outputEdges(const shared::rendergraph::RenderGraph& outputGraph,
                   const RenderParams& params);

  void renderEdgeTripGeom(const shared::rendergraph::RenderGraph& outG,
                          const shared::linegraph::LineEdge* e,
                          const RenderParams& params);

  bool needsDirMarker(const shared::linegraph::LineEdge* e,
                      const util::geo::PolyLine<double>& center,
                      const shared::linegraph::Line* line);

  bool hasSharpAngle(const shared::linegraph::LineEdge* e,
                     const util::geo::PolyLine<double>& center,
                     const shared::linegraph::Line* line);

  bool edgeHasSharpAngle(const util::geo::PolyLine<double>& center,
                         const shared::linegraph::LineEdge* e,
                         const shared::linegraph::Line* line,
                         bool markAdjacent);

  void renderNodeConnections(const shared::rendergraph::RenderGraph& outG,
                             const shared::linegraph::LineNode* n,
                             const RenderParams& params);

  void renderLinePart(const util::geo::PolyLine<double> p, double width,
                      const shared::linegraph::Line& line,
                      const std::string& css,
                      const std::string& oCss);

  void renderArrowHead(const util::geo::PolyLine<double>& p, double width,
                       bool flipDir = false, bool atStart = false);

  void renderDelegates(const shared::rendergraph::RenderGraph& outG,
                       const RenderParams& params);

  void renderNodeFronts(const shared::rendergraph::RenderGraph& outG,
                        const RenderParams& params);
  void renderBackground(const RenderParams& params);

  void renderLandmarks(
      const shared::rendergraph::RenderGraph& g,
      const std::vector<shared::rendergraph::Landmark>& landmarks,
      const RenderParams& params);
  void renderMe(const shared::rendergraph::RenderGraph& g,
                label::Labeller& labeller,
                const RenderParams& params);

  void renderLineLabels(const label::Labeller& lbler,
                        const RenderParams& params);

  void renderStationLabels(const label::Labeller& lbler,
                           const RenderParams& params);

  void renderTerminusLabels(const shared::rendergraph::RenderGraph& g,
                            const label::Labeller& lbler,
                            const RenderParams& params);

  std::multiset<InnerClique> getInnerCliques(
      const shared::linegraph::LineNode* n,
      std::vector<shared::rendergraph::InnerGeom> geoms, size_t level) const;

  void renderClique(const InnerClique& c,
                    const shared::linegraph::LineNode* node);

  bool isNextTo(const shared::rendergraph::InnerGeom& a,
                const shared::rendergraph::InnerGeom& b) const;
  bool hasSameOrigin(const shared::rendergraph::InnerGeom& a,
                     const shared::rendergraph::InnerGeom& b) const;

  size_t getNextPartner(const InnerClique& forGeom,
                        const std::vector<shared::rendergraph::InnerGeom>& pool,
                        size_t level) const;

  std::string getLineClass(const std::string& id) const;
};
}  // namespace output
}  // namespace transitmapper

#endif  // TRANSITMAP_OUTPUT_SVGRENDERER_H_

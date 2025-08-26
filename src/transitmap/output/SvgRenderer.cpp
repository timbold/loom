// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdint.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <ostream>
#include <unordered_map>

#include "shared/linegraph/Line.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/label/Labeller.h"
#include "transitmap/output/SvgRenderer.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using shared::linegraph::Line;
using shared::linegraph::LineNode;
using shared::rendergraph::InnerGeom;
using shared::rendergraph::RenderGraph;
using transitmapper::label::Labeller;
using transitmapper::label::StationLabel;
using transitmapper::output::InnerClique;
using transitmapper::output::SvgRenderer;
using util::geo::DPoint;
using util::geo::DPolygon;
using util::geo::LinePoint;
using util::geo::LinePointCmp;
using util::geo::Polygon;
using util::geo::PolyLine;

// _____________________________________________________________________________
SvgRenderer::SvgRenderer(std::ostream *o, const config::Config *cfg)
    : _o(o), _w(o, true), _cfg(cfg) {}

// _____________________________________________________________________________
void SvgRenderer::print(const RenderGraph &outG) {
  std::map<std::string, std::string> params;
  RenderParams rparams;

  auto box = outG.getBBox();

  box = util::geo::pad(box, outG.getMaxLineNum() *
                                (_cfg->lineWidth + _cfg->lineSpacing));

  Labeller labeller(_cfg);
  for (const auto &lm : outG.getLandmarks()) {
    double half = (lm.size / _cfg->outputResolution) / 2.0;
    util::geo::Box<double> lmBox(
        DPoint(lm.pos.getX() - half, lm.pos.getY() - half),
        DPoint(lm.pos.getX() + half, lm.pos.getY() + half));
    box = util::geo::extendBox(lmBox, box);
    labeller.addLandmark(lmBox);
  }
  if (_cfg->renderLabels) {
    LOGTO(DEBUG, std::cerr) << "Rendering labels...";
    labeller.label(outG, _cfg->dontLabelDeg2);
    box = util::geo::extendBox(labeller.getBBox(), box);
  }

  DPoint ll(box.getLowerLeft().getX() - _cfg->paddingLeft,
            box.getLowerLeft().getY() - _cfg->paddingBottom);
  DPoint ur(box.getUpperRight().getX() + _cfg->paddingRight,
            box.getUpperRight().getY() + _cfg->paddingTop);
  box = util::geo::Box<double>(ll, ur);

  if (_cfg->ratio > 0) {
    double curWidth = box.getUpperRight().getX() - box.getLowerLeft().getX();
    double curHeight = box.getUpperRight().getY() - box.getLowerLeft().getY();
    double desiredWidth = curHeight * _cfg->ratio;
    if (desiredWidth > curWidth) {
      double pad = (desiredWidth - curWidth) / 2.0;
      DPoint nll(box.getLowerLeft().getX() - pad, box.getLowerLeft().getY());
      DPoint nur(box.getUpperRight().getX() + pad, box.getUpperRight().getY());
      box = util::geo::Box<double>(nll, nur);
    } else if (desiredWidth < curWidth) {
      double desiredHeight = curWidth / _cfg->ratio;
      double pad = (desiredHeight - curHeight) / 2.0;
      DPoint nll(box.getLowerLeft().getX(), box.getLowerLeft().getY() - pad);
      DPoint nur(box.getUpperRight().getX(), box.getUpperRight().getY() + pad);
      box = util::geo::Box<double>(nll, nur);
    }
  }

  if (_cfg->tlRatio > 0) {
    double curWidth = box.getUpperRight().getX() - box.getLowerLeft().getX();
    double curHeight = box.getUpperRight().getY() - box.getLowerLeft().getY();
    double desiredWidth = curHeight * _cfg->tlRatio;
    if (desiredWidth > curWidth) {
      double pad = desiredWidth - curWidth;
      DPoint nll(box.getLowerLeft().getX() - pad, box.getLowerLeft().getY());
      DPoint nur(box.getUpperRight().getX(), box.getUpperRight().getY());
      box = util::geo::Box<double>(nll, nur);
    } else if (desiredWidth < curWidth) {
      double desiredHeight = curWidth / _cfg->tlRatio;
      double pad = desiredHeight - curHeight;
      DPoint nll(box.getLowerLeft().getX(), box.getLowerLeft().getY());
      DPoint nur(box.getUpperRight().getX(), box.getUpperRight().getY() + pad);
      box = util::geo::Box<double>(nll, nur);
    }
  }

  if (!_cfg->worldFilePath.empty()) {
    std::ofstream file;
    file.open(_cfg->worldFilePath);
    if (file) {
      file << 1 / _cfg->outputResolution << std::endl
           << 0 << std::endl
           << 0 << std::endl
           << -1 / _cfg->outputResolution << std::endl
           << std::fixed << box.getLowerLeft().getX() << std::endl
           << box.getUpperRight().getY() << std::endl;
      file.close();
    }
  }

  rparams.xOff = box.getLowerLeft().getX();
  rparams.yOff = box.getLowerLeft().getY();

  rparams.width = box.getUpperRight().getX() - rparams.xOff;
  rparams.height = box.getUpperRight().getY() - rparams.yOff;

  rparams.width *= _cfg->outputResolution;
  rparams.height *= _cfg->outputResolution;

  auto latLngLL = util::geo::webMercToLatLng<double>(box.getLowerLeft().getX(),
                                                     box.getLowerLeft().getY());
  auto latLngUR = util::geo::webMercToLatLng<double>(
      box.getUpperRight().getX(), box.getUpperRight().getY());

  params["latlng-box"] = std::to_string(latLngLL.getX()) + "," +
                         std::to_string(latLngLL.getY()) + "," +
                         std::to_string(latLngUR.getX()) + "," +
                         std::to_string(latLngUR.getY());

  params["width"] = std::to_string(rparams.width);
  params["height"] = std::to_string(rparams.height);
  params["viewBox"] = "0 0 " + std::to_string(rparams.width) + " " +
                      std::to_string(rparams.height);
  params["xmlns"] = "http://www.w3.org/2000/svg";
  params["xmlns:xlink"] = "http://www.w3.org/1999/xlink";

  *_o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  *_o << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
         "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">";

  LOGTO(DEBUG, std::cerr) << "Rendering edges...";
  if (_cfg->renderEdges) {
    outputEdges(outG, rparams);
  }
  _w.openTag("svg", params);

  _w.openTag("defs");

  LOGTO(DEBUG, std::cerr) << "Rendering markers...";
  for (auto const &m : _markers) {
    params.clear();
    params["id"] = m.name;
    params["orient"] = "auto";
    params["markerWidth"] = "20";
    params["markerHeight"] = "4";
    params["refY"] = "0.5";
    params["refX"] = "0";

    _w.openTag("marker", params);

    params.clear();
    params["d"] = m.path;
    params["fill"] = m.color;
    ;

    _w.openTag("path", params);

    _w.closeTag();
    _w.closeTag();
  }

  _w.closeTag();

  LOGTO(DEBUG, std::cerr) << "Rendering nodes...";
  for (auto n : outG.getNds()) {
    if (_cfg->renderNodeConnections) {
      renderNodeConnections(outG, n, rparams);
    }
  }

  LOGTO(DEBUG, std::cerr) << "Writing edges...";
  renderDelegates(outG, rparams);

  LOGTO(DEBUG, std::cerr) << "Writing nodes...";
  outputNodes(outG, rparams);
  if (_cfg->renderNodeFronts) {
    renderNodeFronts(outG, rparams);
  }

  renderLandmarks(outG, rparams);

  LOGTO(DEBUG, std::cerr) << "Writing labels...";
  if (_cfg->renderLabels) {
    renderLineLabels(labeller, rparams);

    if (_cfg->renderRouteLabels) {
      renderTerminusLabels(outG, labeller, rparams);
    }

    renderStationLabels(labeller, rparams);
  }

  _w.closeTags();
}

// _____________________________________________________________________________
void SvgRenderer::outputNodes(const RenderGraph &outG,
                              const RenderParams &rparams) {
  _w.openTag("g");
  for (auto n : outG.getNds()) {
    std::map<std::string, std::string> params;

    if (_cfg->renderStations && n->pl().stops().size() > 0 &&
        n->pl().fronts().size() > 0) {
      params["stroke"] =
          (_cfg->highlightTerminals && RenderGraph::isTerminus(n)) ? "#9E9898"
                                                                   : "black";
      params["stroke-width"] =
          util::toString((_cfg->lineWidth / 2) * _cfg->outputResolution);
      params["fill"] = (_cfg->highlightTerminals && RenderGraph::isTerminus(n))
                           ? "black"
                           : "white";

      for (const auto &geom : outG.getStopGeoms(n, _cfg->tightStations, 32)) {
        printPolygon(geom, params, rparams);
      }
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::renderNodeFronts(const RenderGraph &outG,
                                   const RenderParams &rparams) {
  _w.openTag("g");
  for (auto n : outG.getNds()) {
    std::string color = n->pl().stops().size() > 0 ? "red" : "black";
    for (auto &f : n->pl().fronts()) {
      const PolyLine<double> p = f.geom;
      std::stringstream style;
      style << "fill:none;stroke:" << color
            << ";stroke-linejoin: "
               "miter;stroke-linecap:round;stroke-opacity:0.9;stroke-width:1";
      std::map<std::string, std::string> params;
      params["style"] = style.str();
      printLine(p, params, rparams);

      DPoint a = p.getPointAt(.5).p;

      std::stringstream styleA;
      styleA << "fill:none;stroke:" << color
             << ";stroke-linejoin: "
                "miter;stroke-linecap:round;stroke-opacity:1;stroke-width:.5";
      params["style"] = styleA.str();

      printLine(PolyLine<double>(*n->pl().getGeom(), a), params, rparams);
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::renderLandmarks(const RenderGraph &g,
                                  const RenderParams &rparams) {
  std::map<std::string, std::string> iconIds;
  size_t id = 0;

  _w.openTag("defs");
  for (const auto &lm : g.getLandmarks()) {
    auto it = iconIds.find(lm.icon);
    if (it == iconIds.end()) {
      std::string idStr = "lmk" + util::toString(id++);
      iconIds[lm.icon] = idStr;
      std::string svg = lm.icon;
      size_t pos = svg.find("<svg");
      if (pos != std::string::npos) {
        size_t end = svg.find('>', pos);
        if (end != std::string::npos) {
          svg.insert(end, " id=\"" + idStr + "\"");
        }
        *_o << svg;
      } else {
        *_o << "<svg id=\"" << idStr
            << "\" xmlns=\"http://www.w3.org/2000/svg\">" << svg << "</svg>";
      }
    }
  }
  _w.closeTag();

  _w.openTag("g");
  for (const auto &lm : g.getLandmarks()) {
    auto it = iconIds.find(lm.icon);
    if (it == iconIds.end())
      continue;

    double half = lm.size / 2.0;
    double x = (lm.pos.getX() - rparams.xOff) * _cfg->outputResolution - half;
    double y = rparams.height -
               (lm.pos.getY() - rparams.yOff) * _cfg->outputResolution - half;

    _w.openTag("use", {{"xlink:href", "#" + it->second},
                       {"x", util::toString(x)},
                       {"y", util::toString(y)},
                       {"width", util::toString(lm.size)},
                       {"height", util::toString(lm.size)}});
    _w.closeTag();
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::outputEdges(const RenderGraph &outG,
                              const RenderParams &rparams) {
  struct cmp {
    bool operator()(const LineNode *lhs, const LineNode *rhs) const {
      return lhs->getAdjList().size() > rhs->getAdjList().size() ||
             (lhs->getAdjList().size() == rhs->getAdjList().size() &&
              RenderGraph::getConnCardinality(lhs) >
                  RenderGraph::getConnCardinality(rhs)) ||
             (lhs->getAdjList().size() == rhs->getAdjList().size() &&
              lhs > rhs);
    }
  };

  struct cmpEdge {
    bool operator()(const shared::linegraph::LineEdge *lhs,
                    const shared::linegraph::LineEdge *rhs) const {
      return lhs->pl().getLines().size() < rhs->pl().getLines().size() ||
             (lhs->pl().getLines().size() == rhs->pl().getLines().size() &&
              lhs < rhs);
    }
  };

  std::set<const LineNode *, cmp> nodesOrdered;
  std::set<const shared::linegraph::LineEdge *, cmpEdge> edgesOrdered;
  for (auto nd : outG.getNds())
    nodesOrdered.insert(nd);

  std::set<const shared::linegraph::LineEdge *> rendered;

  for (const auto n : nodesOrdered) {
    edgesOrdered.insert(n->getAdjList().begin(), n->getAdjList().end());

    for (const auto *e : edgesOrdered) {
      if (rendered.insert(e).second)
        renderEdgeTripGeom(outG, e, rparams);
    }
  }
}

// _____________________________________________________________________________
void SvgRenderer::renderNodeConnections(const RenderGraph &outG,
                                        const LineNode *n,
                                        const RenderParams &rparams) {
  UNUSED(rparams);
  auto geoms = outG.innerGeoms(n, _cfg->innerGeometryPrecision);

  for (auto &clique : getInnerCliques(n, geoms, 9999))
    renderClique(clique, n);
}

// _____________________________________________________________________________
std::multiset<InnerClique>
SvgRenderer::getInnerCliques(const shared::linegraph::LineNode *n,
                             std::vector<InnerGeom> pool, size_t level) const {
  std::multiset<InnerClique> ret;

  // start with the first geom in pool
  while (!pool.empty()) {
    InnerClique cur(n, pool.front());
    pool.erase(pool.begin());

    size_t p;
    while ((p = getNextPartner(cur, pool, level)) < pool.size()) {
      cur.geoms.push_back(pool[p]);
      pool.erase(pool.begin() + p);
    }

    ret.insert(cur);
  }

  return ret;
}

// _____________________________________________________________________________
size_t SvgRenderer::getNextPartner(const InnerClique &forClique,
                                   const std::vector<InnerGeom> &pool,
                                   size_t level) const {
  for (size_t i = 0; i < pool.size(); i++) {
    const auto &ic = pool[i];
    for (auto &ciq : forClique.geoms) {
      if (isNextTo(ic, ciq) || (level > 1 && hasSameOrigin(ic, ciq))) {
        return i;
      }
    }
  }

  return pool.size();
}

// _____________________________________________________________________________
bool SvgRenderer::isNextTo(const InnerGeom &a, const InnerGeom &b) const {
  double THRESHOLD = 0.5 * M_PI + 0.1;

  if (!a.from.edge)
    return false;
  if (!b.from.edge)
    return false;
  if (!a.to.edge)
    return false;
  if (!b.to.edge)
    return false;

  auto nd = RenderGraph::sharedNode(a.from.edge, a.to.edge);

  assert(a.from.edge);
  assert(b.from.edge);
  assert(a.to.edge);
  assert(b.to.edge);

  bool aFromInv = a.from.edge->getTo() == nd;
  bool bFromInv = b.from.edge->getTo() == nd;
  bool aToInv = a.to.edge->getTo() == nd;
  bool bToInv = b.to.edge->getTo() == nd;

  int aSlotFrom = !aFromInv
                      ? a.slotFrom
                      : (a.from.edge->pl().getLines().size() - 1 - a.slotFrom);
  int aSlotTo =
      !aToInv ? a.slotTo : (a.to.edge->pl().getLines().size() - 1 - a.slotTo);
  int bSlotFrom = !bFromInv
                      ? b.slotFrom
                      : (b.from.edge->pl().getLines().size() - 1 - b.slotFrom);
  int bSlotTo =
      !bToInv ? b.slotTo : (b.to.edge->pl().getLines().size() - 1 - b.slotTo);

  if (a.from.edge == b.from.edge && a.to.edge == b.to.edge) {
    if ((aSlotFrom - bSlotFrom == 1 && bSlotTo - aSlotTo == 1) ||
        (bSlotFrom - aSlotFrom == 1 && aSlotTo - bSlotTo == 1)) {
      return true;
      double ang1 = fabs(util::geo::angBetween(a.geom.front(), a.geom.back()));
      double ang2 = fabs(util::geo::angBetween(b.geom.front(), b.geom.back()));

      return ang1 > THRESHOLD && ang2 > THRESHOLD;
    }
  }

  if (a.to.edge == b.from.edge && a.from.edge == b.to.edge) {
    if ((aSlotFrom - bSlotTo == 1 && bSlotFrom - aSlotTo == 1) ||
        (bSlotTo - aSlotFrom == 1 && aSlotTo - bSlotFrom == 1)) {
      return true;
      double ang1 = fabs(util::geo::angBetween(a.geom.front(), a.geom.back()));
      double ang2 = fabs(util::geo::angBetween(b.geom.front(), b.geom.back()));

      return ang1 > THRESHOLD && ang2 > THRESHOLD;
    }
  }

  return false;
}

// _____________________________________________________________________________
bool SvgRenderer::hasSameOrigin(const InnerGeom &a, const InnerGeom &b) const {
  if (a.from.edge == b.from.edge) {
    return a.slotFrom == b.slotFrom;
  }
  if (a.to.edge == b.from.edge) {
    return a.slotTo == b.slotFrom;
  }
  if (a.to.edge == b.to.edge) {
    return a.slotTo == b.slotTo;
  }
  if (a.from.edge == b.to.edge) {
    return a.slotFrom == b.slotTo;
  }

  return false;
}

// _____________________________________________________________________________
void SvgRenderer::renderClique(const InnerClique &cc, const LineNode *n) {
  _innerDelegates.push_back(
      std::map<uintptr_t, std::vector<OutlinePrintPair>>());
  std::multiset<InnerClique> renderCliques = getInnerCliques(n, cc.geoms, 0);
  for (const auto &c : renderCliques) {
    // the longest geom will be the ref geom
    InnerGeom ref = c.geoms[0];
    for (size_t i = 1; i < c.geoms.size(); i++) {
      if (c.geoms[i].geom.getLength() > ref.geom.getLength())
        ref = c.geoms[i];
    }

    for (size_t i = 0; i < c.geoms.size(); i++) {
      PolyLine<double> pl = c.geoms[i].geom;

      if (ref.geom.getLength() >
          (_cfg->lineWidth + 2 * _cfg->outlineWidth + _cfg->lineSpacing) * 4) {
        double off =
            -(_cfg->lineWidth + _cfg->lineSpacing + 2 * _cfg->outlineWidth) *
            (static_cast<int>(c.geoms[i].slotFrom) -
             static_cast<int>(ref.slotFrom));

        if (ref.from.edge->getTo() == n)
          off = -off;

        pl = ref.geom.offsetted(off);

        if (pl.getLength() / c.geoms[i].geom.getLength() > 1.5)
          pl = c.geoms[i].geom;

        std::set<LinePoint<double>, LinePointCmp<double>> a;
        std::set<LinePoint<double>, LinePointCmp<double>> b;

        if (ref.from.edge)
          a = n->pl().frontFor(ref.from.edge)->geom.getIntersections(pl);
        if (ref.to.edge)
          b = n->pl().frontFor(ref.to.edge)->geom.getIntersections(pl);

        if (a.size() > 0 && b.size() > 0) {
          pl = pl.getSegment(a.begin()->totalPos, b.begin()->totalPos);
        } else if (a.size() > 0) {
          pl = pl.getSegment(a.begin()->totalPos, 1);
        } else if (b.size() > 0) {
          pl = pl.getSegment(0, b.begin()->totalPos);
        }
      }

      std::stringstream styleOutlineCropped;
      styleOutlineCropped << "fill:none;stroke:#000000";

      styleOutlineCropped << ";stroke-linecap:butt;stroke-width:"
                          << (_cfg->lineWidth + _cfg->outlineWidth) *
                                 _cfg->outputResolution;
      Params paramsOutlineCropped;
      paramsOutlineCropped["style"] = styleOutlineCropped.str();
      paramsOutlineCropped["class"] += " inner-geom-outline";
      paramsOutlineCropped["class"] +=
          " " + getLineClass(c.geoms[i].from.line->id());

      std::stringstream styleStr;
      styleStr << "fill:none;stroke:#" << c.geoms[i].from.line->color();

      styleStr << ";stroke-linecap:round;stroke-opacity:1;stroke-width:"
               << _cfg->lineWidth * _cfg->outputResolution;
      Params params;
      params["style"] = styleStr.str();
      params["class"] += " inner-geom ";
      params["class"] += " " + getLineClass(c.geoms[i].from.line->id());

      _innerDelegates.back()[(uintptr_t)c.geoms[i].from.line].push_back(
          OutlinePrintPair(PrintDelegate(params, pl),
                           PrintDelegate(paramsOutlineCropped, pl)));
    }
  }
}

// _____________________________________________________________________________
void SvgRenderer::renderLinePart(const PolyLine<double> p, double width,
                                 const Line &line, const std::string &css,
                                 const std::string &oCss) {
  renderLinePart(p, width, line, css, oCss, "", "");
}

// _____________________________________________________________________________
void SvgRenderer::renderLinePart(const PolyLine<double> p, double width,
                                 const Line &line, const std::string &css,
                                 const std::string &oCss,
                                 const std::string &endMarker,
                                 const std::string &startMarker) {
  std::stringstream styleOutline;
  styleOutline << "fill:none;stroke:#000000;stroke-linecap:round;stroke-width:"
               << (width + _cfg->outlineWidth) * _cfg->outputResolution << ";"
               << oCss;
  Params paramsOutline;
  paramsOutline["style"] = styleOutline.str();
  paramsOutline["class"] = "transit-edge-outline " + getLineClass(line.id());

  std::stringstream styleStr;
  styleStr << "fill:none;stroke:#" << line.color() << ";" << css;

  if (!endMarker.empty()) {
    styleStr << ";marker-end:url(#" << endMarker << ")";
  }
  if (!startMarker.empty()) {
    styleStr << ";marker-start:url(#" << startMarker << ")";
  }

  styleStr << ";stroke-linecap:round;stroke-opacity:1;stroke-width:"
           << width * _cfg->outputResolution;
  Params params;
  params["style"] = styleStr.str();
  params["class"] = "transit-edge " + getLineClass(line.id());

  _delegates[0].insert(_delegates[0].begin(),
                       OutlinePrintPair(PrintDelegate(params, p),
                                        PrintDelegate(paramsOutline, p)));
}

// _____________________________________________________________________________
bool SvgRenderer::needsDirMarker(const shared::linegraph::LineEdge *e,
                                 const PolyLine<double> &center,
                                 const shared::linegraph::Line *line) {
  auto it = _forceDirMarker.find(line);
  if (it != _forceDirMarker.end()) {
    auto &s = it->second;
    if (s.find(e) != s.end()) {
      s.erase(e);
      return true;
    }
  }

  if (e->pl().getLines().size() >= _cfg->crowdedLineThresh) {
    return true;
  }

  const auto &pts = center.getLine();

  for (size_t i = 1; i + 1 < pts.size(); ++i) {
    const DPoint &a = pts[i - 1];
    const DPoint &b = pts[i];
    const DPoint &c = pts[i + 1];
    double ux = b.getX() - a.getX();
    double uy = b.getY() - a.getY();
    double vx = c.getX() - b.getX();
    double vy = c.getY() - b.getY();
    double dot = ux * vx + uy * vy;
    double lu = std::sqrt(ux * ux + uy * uy);
    double lv = std::sqrt(vx * vx + vy * vy);
    if (lu == 0 || lv == 0)
      continue;
    double cosang = dot / (lu * lv);
    cosang = std::max(-1.0, std::min(1.0, cosang));
    double ang = std::acos(cosang);
    if (ang > _cfg->sharpTurnAngle) {
      return true;
    }
  }

  bool sharp = false;
  double checkDist = 10.0;
  auto checkNode = [&](const shared::linegraph::LineNode *n, bool fromStart) {
    PolyLine<double> plE(*e->pl().getGeom());
    DPoint base = fromStart ? plE.front() : plE.back();
    double lenE = plE.getLength();
    DPoint otherE;
    if (fromStart) {
      otherE = plE.getPointAtDist(std::min(checkDist, lenE)).p;
    } else {
      otherE = plE.getPointAtDist(std::max(0.0, lenE - checkDist)).p;
    }
    double ux = otherE.getX() - base.getX();
    double uy = otherE.getY() - base.getY();
    for (auto ne : n->getAdjList()) {
      if (ne == e)
        continue;
      if (!ne->pl().hasLine(line))
        continue;
      PolyLine<double> plN(*ne->pl().getGeom());
      double lenN = plN.getLength();
      DPoint baseN = (ne->getFrom() == n) ? plN.front() : plN.back();
      DPoint otherN;
      if (ne->getFrom() == n) {
        otherN = plN.getPointAtDist(std::min(checkDist, lenN)).p;
      } else {
        otherN = plN.getPointAtDist(std::max(0.0, lenN - checkDist)).p;
      }
      double vx = otherN.getX() - baseN.getX();
      double vy = otherN.getY() - baseN.getY();
      double dot = ux * vx + uy * vy;
      double lu = std::sqrt(ux * ux + uy * uy);
      double lv = std::sqrt(vx * vx + vy * vy);
      if (lu == 0 || lv == 0)
        continue;
      double cosang = dot / (lu * lv);
      cosang = std::max(-1.0, std::min(1.0, cosang));
      double ang = std::acos(cosang);
      if (ang > _cfg->sharpTurnAngle) {
        sharp = true;
        _forceDirMarker[line].insert(ne);
      }
    }
  };

  checkNode(e->getFrom(), true);
  checkNode(e->getTo(), false);
  if (sharp) {
    return true;
  }

  if (_edgesSinceMarker[line] >= 3) {
    return true;
  }

  return false;
}

// _____________________________________________________________________________
void SvgRenderer::renderEdgeTripGeom(const RenderGraph &outG,
                                     const shared::linegraph::LineEdge *e,
                                     const RenderParams &rparams) {
  UNUSED(rparams);
  const shared::linegraph::NodeFront *nfTo = e->getTo()->pl().frontFor(e);
  const shared::linegraph::NodeFront *nfFrom = e->getFrom()->pl().frontFor(e);

  assert(nfTo);
  assert(nfFrom);

  PolyLine<double> center(*e->pl().getGeom());

  double lineW = _cfg->lineWidth;
  double outlineW = _cfg->outlineWidth;
  double lineSpc = _cfg->lineSpacing;
  double offsetStep = lineW + 2.0 * outlineW + lineSpc;
  double oo = outG.getTotalWidth(e);

  double o = oo;

  for (size_t i = 0; i < e->pl().getLines().size(); i++) {
    const auto &lo = e->pl().lineOccAtPos(i);

    const Line *line = lo.line;
    PolyLine<double> p = center;

    if (p.getLength() < 0.01)
      continue;

    double offset = -(o - oo / 2.0 - (2.0 * outlineW + _cfg->lineWidth) / 2.0);

    p.offsetPerp(offset);

    auto iSects = nfTo->geom.getIntersections(p);
    if (iSects.size() > 0) {
      p = p.getSegment(0, iSects.begin()->totalPos);
    } else {
      p << nfTo->geom.projectOn(p.back()).p;
    }

    auto iSects2 = nfFrom->geom.getIntersections(p);
    if (iSects2.size() > 0) {
      p = p.getSegment(iSects2.begin()->totalPos, 1);
    } else {
      p >> nfFrom->geom.projectOn(p.front()).p;
    }

    double arrowLength = (_cfg->lineWidth * 2.5);

    std::string css, oCss;

    if (!lo.style.isNull()) {
      css = lo.style.get().getCss();
      oCss = lo.style.get().getOutlineCss();
    }

    bool needMarker = _cfg->renderDirMarkers && needsDirMarker(e, center, line);
    bool drawMarker = needMarker && center.getLength() > arrowLength * 3;

    if (drawMarker) {
      _edgesSinceMarker[line] = 0;
    } else {
      _edgesSinceMarker[line]++;
    }

    if (drawMarker) {
      if (lo.direction == 0 && !_cfg->renderBiDirMarker) {
        renderLinePart(p, lineW, *line, css, oCss);
      } else {
        std::stringstream markerName;
        markerName << e << ":" << line << ":" << i;

        std::string markerPathMale = getMarkerPathMale(lineW);
        EndMarker emm(markerName.str() + "_m", "white", markerPathMale, lineW,
                      lineW);

        _markers.push_back(emm);

        PolyLine<double> firstPart = p.getSegmentAtDist(0, p.getLength() / 2);
        PolyLine<double> secondPart =
            p.getSegmentAtDist(p.getLength() / 2, p.getLength());
        PolyLine<double> revSecond = secondPart.reversed();

        if (lo.direction == 0) {
          double mid = p.getLength() / 2;
          double tailWorld = 15.0 / _cfg->outputResolution;
          double tailStart = mid - tailWorld / 2;
          double tailEnd = mid + tailWorld / 2;

          PolyLine<double> firstHalf = p.getSegmentAtDist(0, mid);
          PolyLine<double> secondHalf = p.getSegmentAtDist(mid, p.getLength());
          PolyLine<double> revFirstHalf = firstHalf.reversed();

          if (_cfg->renderMarkersTail) {
            EndMarker emmTail(markerName.str() + "_mt", "white", markerPathMale,
                              lineW, lineW);
            _markers.push_back(emmTail);

            PolyLine<double> tailToStart =
                p.getSegmentAtDist(tailStart, mid).reversed();
            PolyLine<double> tailToEnd = p.getSegmentAtDist(mid, tailEnd);
            renderLinePart(tailToStart, lineW, *line, "stroke:black",
                           "stroke:none", markerName.str() + "_mt");
            renderLinePart(tailToEnd, lineW, *line, "stroke:black",
                           "stroke:none", markerName.str() + "_mt");
          }

          renderLinePart(revFirstHalf, lineW, *line, css, oCss,
                         markerName.str() + "_m");
          renderLinePart(secondHalf, lineW, *line, css, oCss,
                         markerName.str() + "_m");
        } else if (lo.direction == e->getTo()) {
          double tailWorld = 15.0 / _cfg->outputResolution;
          if (_cfg->renderMarkersTail) {
            double tailStart = std::max(0.0, firstPart.getLength() - tailWorld);
            PolyLine<double> tail =
                firstPart.getSegmentAtDist(tailStart, firstPart.getLength());

            renderLinePart(tail, lineW, *line, "stroke:black", "stroke:none");
          }
          renderLinePart(firstPart, lineW, *line, css, oCss,
                         markerName.str() + "_m");
          renderLinePart(revSecond, lineW, *line, css, oCss);
        } else {
          double tailWorld = 15.0 / _cfg->outputResolution;
          if (_cfg->renderMarkersTail) {
            double tailStart = std::max(0.0, revSecond.getLength() - tailWorld);
            PolyLine<double> tail =
                revSecond.getSegmentAtDist(tailStart, revSecond.getLength());
            renderLinePart(tail, lineW, *line, "stroke:black", "stroke:none");
          }
          renderLinePart(revSecond, lineW, *line, css, oCss,
                         markerName.str() + "_m");
          renderLinePart(firstPart, lineW, *line, css, oCss);
        }
      }
    } else {
      renderLinePart(p, lineW, *line, css, oCss);
    }

    o -= offsetStep;
  }
}

// _____________________________________________________________________________
std::string SvgRenderer::getMarkerPathMale(double w) const {
  UNUSED(w);
  return "M0,0 V1 H.5 L1.3,.5 L.5,0 Z";
}

// _____________________________________________________________________________
void SvgRenderer::renderDelegates(const RenderGraph &outG,
                                  const RenderParams &rparams) {
  UNUSED(outG);
  for (auto &a : _delegates) {
    _w.openTag("g");
    for (auto &pd : a.second) {
      if (_cfg->outlineWidth > 0) {
        printLine(pd.back.second, pd.back.first, rparams);
      }
      printLine(pd.front.second, pd.front.first, rparams);
    }
    _w.closeTag();
  }

  for (auto &a : _innerDelegates) {
    _w.openTag("g");
    for (auto &b : a) {
      for (auto &pd : b.second) {
        if (_cfg->outlineWidth > 0) {
          printLine(pd.back.second, pd.back.first, rparams);
        }
      }
      for (auto &pd : b.second) {
        printLine(pd.front.second, pd.front.first, rparams);
      }
    }
    _w.closeTag();
  }
}

// _____________________________________________________________________________
void SvgRenderer::printPoint(const DPoint &p, const std::string &style,
                             const RenderParams &rparams) {
  std::map<std::string, std::string> params;
  params["cx"] =
      std::to_string((p.getX() - rparams.xOff) * _cfg->outputResolution);
  params["cy"] = std::to_string(rparams.height - (p.getY() - rparams.yOff) *
                                                     _cfg->outputResolution);
  params["r"] = "2";
  params["style"] = style;
  _w.openTag("circle", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::printLine(const PolyLine<double> &l, const std::string &style,
                            const RenderParams &rparams) {
  std::map<std::string, std::string> params;
  params["style"] = style;
  printLine(l, params, rparams);
}

// _____________________________________________________________________________
void SvgRenderer::printLine(const PolyLine<double> &l,
                            const std::map<std::string, std::string> &ps,
                            const RenderParams &rparams) {
  std::map<std::string, std::string> params = ps;
  std::stringstream points;

  for (auto &p : l.getLine()) {
    points << " " << (p.getX() - rparams.xOff) * _cfg->outputResolution << ","
           << rparams.height -
                  (p.getY() - rparams.yOff) * _cfg->outputResolution;
  }

  params["points"] = points.str();

  _w.openTag("polyline", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::printPolygon(const Polygon<double> &g,
                               const std::map<std::string, std::string> &ps,
                               const RenderParams &rparams) {
  std::map<std::string, std::string> params = ps;
  std::stringstream points;

  for (auto &p : g.getOuter()) {
    points << " " << (p.getX() - rparams.xOff) * _cfg->outputResolution << ","
           << rparams.height -
                  (p.getY() - rparams.yOff) * _cfg->outputResolution;
  }

  params["points"] = points.str();
  params["class"] = "station-poly";

  _w.openTag("polygon", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::printCircle(const DPoint &center, double rad,
                              const std::string &style,
                              const RenderParams &rparams) {
  std::map<std::string, std::string> params;
  params["style"] = style;
  printCircle(center, rad, params, rparams);
}

// _____________________________________________________________________________
void SvgRenderer::printCircle(const DPoint &center, double rad,
                              const std::map<std::string, std::string> &ps,
                              const RenderParams &rparams) {
  std::map<std::string, std::string> params = ps;
  std::stringstream points;

  params["cx"] =
      std::to_string((center.getX() - rparams.xOff) * _cfg->outputResolution);
  params["cy"] = std::to_string(
      rparams.height - (center.getY() - rparams.yOff) * _cfg->outputResolution);
  params["r"] = std::to_string(rad * _cfg->outputResolution);

  _w.openTag("circle", params);
  _w.closeTag();
}

// _____________________________________________________________________________
size_t
InnerClique::getNumBranchesIn(const shared::linegraph::LineEdge *edg) const {
  std::set<size_t> slots;
  size_t ret = 0;
  for (const auto &ig : geoms) {
    if (ig.from.edge == edg && !slots.insert(ig.slotFrom).second)
      ret++;
    if (ig.to.edge == edg && !slots.insert(ig.slotTo).second)
      ret++;
  }

  return ret;
}

// _____________________________________________________________________________
void SvgRenderer::renderStationLabels(const Labeller &labeller,
                                      const RenderParams &rparams) {
  _w.openTag("g");
  size_t id = 0;
  for (auto label : labeller.getStationLabels()) {
    std::string shift = "0em";
    std::string textAnchor = "start";
    std::string startOffset = "0";
    auto textPath = label.geom;
    double ang = util::geo::angBetween(textPath.front(), textPath.back());

    if ((fabs(ang) < (3 * M_PI / 2)) && (fabs(ang) > (M_PI / 2))) {
      shift = ".75em";
      textAnchor = "end";
      startOffset = "100%";
      textPath.reverse();
    }

    std::stringstream points;
    std::map<std::string, std::string> pathPars;

    points << "M"
           << (textPath.front().getX() - rparams.xOff) * _cfg->outputResolution
           << " "
           << rparams.height - (textPath.front().getY() - rparams.yOff) *
                                   _cfg->outputResolution;

    for (auto &p : textPath.getLine()) {
      points << " L" << (p.getX() - rparams.xOff) * _cfg->outputResolution
             << " "
             << rparams.height -
                    (p.getY() - rparams.yOff) * _cfg->outputResolution;
    }

    std::string idStr = "stlblp" + util::toString(id);

    pathPars["d"] = points.str();
    pathPars["id"] = idStr;
    id++;

    _w.openTag("defs");
    _w.openTag("path", pathPars);
    _w.closeTag();
    _w.closeTag();

    std::map<std::string, std::string> params;
    params["class"] = "station-label";
    params["font-weight"] = label.bold ? "bold" : "normal";
    params["font-family"] = "TT Norms Pro";
    params["dy"] = shift;
    params["font-size"] =
        util::toString(label.fontSize * _cfg->outputResolution);

    _w.openTag("text", params);
    _w.openTag("textPath", {{"dy", shift},
                            {"xlink:href", "#" + idStr},
                            {"startOffset", startOffset},
                            {"text-anchor", textAnchor}});

    _w.writeText(label.s.name);
    _w.closeTag();
    _w.closeTag();
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::renderLineLabels(const Labeller &labeller,
                                   const RenderParams &rparams) {
  _w.openTag("g");
  size_t id = 0;
  for (auto label : labeller.getLineLabels()) {
    std::string shift = "0em";
    auto textPath = label.geom;
    double ang = util::geo::angBetween(textPath.front(), textPath.back());

    if ((fabs(ang) < (3 * M_PI / 2)) && (fabs(ang) > (M_PI / 2))) {
      shift = ".75em";
      textPath.reverse();
    }

    std::stringstream points;
    std::map<std::string, std::string> pathPars;

    points << "M"
           << (textPath.front().getX() - rparams.xOff) * _cfg->outputResolution
           << " "
           << rparams.height - (textPath.front().getY() - rparams.yOff) *
                                   _cfg->outputResolution;

    for (auto &p : textPath.getLine()) {
      points << " L" << (p.getX() - rparams.xOff) * _cfg->outputResolution
             << " "
             << rparams.height -
                    (p.getY() - rparams.yOff) * _cfg->outputResolution;
    }

    std::string idStr = "textp" + util::toString(id);

    pathPars["d"] = points.str();
    pathPars["id"] = idStr;
    id++;

    _w.openTag("defs");
    _w.openTag("path", pathPars);
    _w.closeTag();
    _w.closeTag();

    std::map<std::string, std::string> params;
    params["class"] = "line-label";
    params["font-weight"] = "bold";
    params["font-family"] = "TT Norms Pro";
    params["dy"] = shift;
    params["font-size"] =
        util::toString(label.fontSize * _cfg->outputResolution);

    _w.openTag("text", params);
    _w.openTag("textPath", {{"dy", shift},
                            {"xlink:href", "#" + idStr},
                            {"text-anchor", "middle"},
                            {"startOffset", "50%"}});

    double dy = 0;
    for (auto line : label.lines) {
      _w.openTag("tspan",
                 {{"fill", "#" + line->color()}, {"dx", util::toString(dy)}});
      dy = (label.fontSize * _cfg->outputResolution) / 3;
      _w.writeText(line->label());
      _w.closeTag();
    }
    _w.closeTag();
    _w.closeTag();
  }
  _w.closeTag();
}

static bool isLightColor(const std::string &hex) {
  if (hex.size() != 6)
    return false;
  auto hexToInt = [](char c) {
    if (c >= '0' && c <= '9')
      return c - '0';
    if (c >= 'a' && c <= 'f')
      return 10 + (c - 'a');
    if (c >= 'A' && c <= 'F')
      return 10 + (c - 'A');
    return 0;
  };
  int r = hexToInt(hex[0]) * 16 + hexToInt(hex[1]);
  int g = hexToInt(hex[2]) * 16 + hexToInt(hex[3]);
  int b = hexToInt(hex[4]) * 16 + hexToInt(hex[5]);

  // Relative luminance formula
  double luminance = 0.299 * r + 0.587 * g + 0.114 * b;
  return luminance > 186; // threshold, tweak if needed
}

void SvgRenderer::renderTerminusLabels(const RenderGraph &g,
                                       const label::Labeller &labeller,
                                       const RenderParams &rparams) {
  _w.openTag("g");
  std::unordered_map<std::string, const StationLabel *> stationLabelMap;
  for (const auto &lbl : labeller.getStationLabels()) {
    stationLabelMap[lbl.s.id] = &lbl;
  }

  for (auto n : g.getNds()) {
    std::set<const Line *> lines;
    std::set<const Line *> seen;
    for (auto e : n->getAdjList()) {
      for (const auto &lo : e->pl().getLines()) {
        if (seen.insert(lo.line).second && g.lineTerminatesAt(n, lo.line)) {
          lines.insert(lo.line);
        }
      }
    }
    if (lines.empty())
      continue;

    double nodeX = n->pl().getGeom()->getX();
    double nodeY = n->pl().getGeom()->getY();

    const StationLabel *sLbl = nullptr;
    if (!n->pl().stops().empty()) {
      const std::string &sid = n->pl().stops().front().id;
      auto it = stationLabelMap.find(sid);
      if (it != stationLabelMap.end()) {
        sLbl = it->second;
      }
    }

    double anchorX = nodeX;
    double anchorY = nodeY;
    bool above = true;
    if (sLbl) {
      double minX = std::numeric_limits<double>::max();
      double maxX = std::numeric_limits<double>::lowest();
      double minY = std::numeric_limits<double>::max();
      double maxY = std::numeric_limits<double>::lowest();
      for (const auto &ln : sLbl->band) {
        for (const auto &p : ln) {
          minX = std::min(minX, p.getX());
          maxX = std::max(maxX, p.getX());
          minY = std::min(minY, p.getY());
          maxY = std::max(maxY, p.getY());
        }
      }

      double centerX = (minX + maxX) / 2;
      double centerY = (minY + maxY) / 2;
      above = centerY > nodeY;

      // Determine label dimensions and rotation to derive a rotation-aware
      // distance from the label center to its outer edge along the vertical
      // axis. This avoids using the axis-aligned bounding box which leads to
      // inconsistent gaps for rotated labels.
      const auto &base = sLbl->band[0];
      const auto &top = sLbl->band[2];
      double dx = base[1].getX() - base[0].getX();
      double dy = base[1].getY() - base[0].getY();
      double width = std::sqrt(dx * dx + dy * dy);
      double angle = std::atan2(dy, dx);
      double hdx = top[0].getX() - base[0].getX();
      double hdy = top[0].getY() - base[0].getY();
      double height = std::sqrt(hdx * hdx + hdy * hdy);

      double vExtent;
      double absTan = std::abs(std::tan(angle));
      if (absTan <= width / height) {
        vExtent = std::abs(height / (2.0 * std::cos(angle)));
      } else {
        vExtent = std::abs(width / (2.0 * std::sin(angle)));
      }

      anchorX = centerX;
      anchorY = above ? centerY + vExtent : centerY - vExtent;
    }

    double x = (anchorX - rparams.xOff) * _cfg->outputResolution;
    double y =
        rparams.height - (anchorY - rparams.yOff) * _cfg->outputResolution;

    double fontSize = _cfg->lineLabelSize * _cfg->outputResolution;
    double pad = fontSize * 0.2;
    double boxH = fontSize + pad * 2;
    double charW = fontSize * 0.6;
    double boxW = 5 * charW + pad * 2; // uniform width for up to 4 chars

    size_t idx = 0;
    // Use a uniform gap to achieve consistent spacing regardless of the
    // orientation of the station label. The gap is configurable to allow
    // tuning without recompilation.
    double gap = _cfg->routeLabelGap * _cfg->outputResolution;
    double startY = above ? y - boxH - gap : y + gap;
    double step = boxH + gap;

    for (auto line : lines) {
      std::string label = line->label();
      double rectX = x - boxW / 2;
      double rectY = above ? startY - idx * step : startY + idx * step;

      std::string fillColor = line->color(); // e.g. "ffcc00"
      std::string textColor = isLightColor(fillColor) ? "black" : "white";

      _w.openTag("rect", {{"x", util::toString(rectX)},
                          {"y", util::toString(rectY)},
                          {"width", util::toString(boxW)},
                          {"height", util::toString(boxH)},
                          {"fill", "#" + fillColor}});
      _w.closeTag();

      _w.openTag("text", {{"class", "line-label"},
                          {"font-weight", "bold"},
                          {"font-family", "TT Norms Pro"},
                          {"text-anchor", "middle"},
                          {"dominant-baseline", "middle"},
                          {"alignment-baseline", "middle"},
                          {"font-size", util::toString(fontSize)},
                          {"fill", textColor},
                          {"x", util::toString(x)},
                          {"y", util::toString(rectY + boxH / 2)}});

      _w.writeText(label);
      _w.closeTag();
      idx++;
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
double InnerClique::getZWeight() const {
  // more weight = more to the bottom

  double BRANCH_WEIGHT = 4;

  double ret = 0;

  ret = geoms.size(); // baseline: threads with more lines to the bottom,
                      // because they are easier to follow

  for (const auto &nf : n->pl().fronts()) {
    ret -= getNumBranchesIn(nf.edge) * BRANCH_WEIGHT;
  }

  return ret;
}

// _____________________________________________________________________________
std::string SvgRenderer::getLineClass(const std::string &id) const {
  auto i = lineClassIds.find(id);
  if (i != lineClassIds.end())
    return "line-" + std::to_string(i->second);

  lineClassIds[id] = ++lineClassId;
  return "line-" + std::to_string(lineClassId);
}

// _____________________________________________________________________________
bool InnerClique::operator<(const InnerClique &rhs) const {
  // more weight = more to the bottom
  return getZWeight() > rhs.getZWeight();
}
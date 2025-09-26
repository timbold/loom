// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdint.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <ostream>
#include <regex>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <vector>
#include <array>

#include "3rdparty/json.hpp"
#include "shared/linegraph/Line.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/label/Labeller.h"
#include "transitmap/output/SvgRenderer.h"
#include "transitmap/util/String.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"
#include "util/geo/Polygon.h"
#include "util/log/Log.h"

using shared::linegraph::Line;
using shared::linegraph::LineNode;
using shared::rendergraph::InnerGeom;
using shared::rendergraph::Landmark;
using shared::rendergraph::RenderGraph;
using transitmapper::config::Config;
using transitmapper::label::Labeller;
using transitmapper::label::StationLabel;
using transitmapper::config::TerminusLabelAnchor;
using transitmapper::output::InnerClique;
using transitmapper::output::SvgRenderer;
using util::geo::DPoint;
using util::geo::DPolygon;
using util::geo::LinePoint;
using util::geo::LinePointCmp;
using util::geo::Polygon;
using util::geo::PolyLine;

static const std::regex scriptRe(
    R"(<\s*(script|foreignObject|iframe)[^>]*>[\s\S]*?<\s*/\s*(script|foreignObject|iframe)\s*>)",
    std::regex::icase);
static const std::regex onAttrRe(R"(\son[\w:-]+\s*=\s*(\"[^\"]*\"|'[^']*'))",
                                 std::regex::icase);
static const std::regex jsHrefRe(
    R"((xlink:href|href)\s*=\s*(\"javascript:[^\"]*\"|'javascript:[^']*'))",
    std::regex::icase);
static const std::regex
    styleTagRe(R"(<\s*style[^>]*>[\s\S]*?<\s*/\s*style\s*>)",
               std::regex::icase);
static const std::regex styleAttrRe(R"(\sstyle\s*=\s*(\"[^\"]*\"|'[^']*'))",
                                    std::regex::icase);
// Allow data:image/... URIs but strip all other data:* URIs
static const std::regex dataUriAttrRe(
    R"(\s[\w:-]+\s*=\s*(\"data:(?!image/)[^\"]*\"|'data:(?!image/)[^']*'))",
    std::regex::icase);

namespace {

constexpr double kBadgePadXFactor = 0.2;
constexpr double kBadgePadTopFactor = 0.04;
constexpr double kBadgePadBottomFactor = 0.16;
constexpr double kBadgeStarGapFactor = 0.14;
constexpr const char kBadgeStarPath[] =
    "M12.231 7.79537C12.777 6.73487 14.307 6.73487 14.852 7.79537L16.691 "
    "11.3687C16.722 11.4279 16.779 11.4692 16.845 11.4797L20.842 12.1101C22.0"
    "28 12.2971 22.5 13.7357 21.652 14.5775L18.791 17.4188C18.744 17.4659 18.7"
    "22 17.5325 18.732 17.5981L19.363 21.5635C19.55 22.7388 18.3131 23.6283 17"
    ".2421 23.0883L13.637 21.2695C13.577 21.2393 13.506 21.2393 13.446 21.2695"
    "L9.84105 23.0883C8.77105 23.6283 7.53402 22.7388 7.72102 21.5635L8.35103 "
    "17.5981C8.36203 17.5325 8.34004 17.4659 8.29204 17.4188L5.43105 14.5775C4"
    ".58305 13.7357 5.05604 12.2971 6.24104 12.1101L10.238 11.4797C10.305 11.4"
    "692 10.362 11.4279 10.392 11.3687L12.231 7.79537Z";
constexpr double kBadgeStarPathCenterX = 13.541525;
constexpr double kBadgeStarPathCenterY = 15.181585;
constexpr double kBadgeStarPathWidth = 17.91695;
constexpr double kBadgeStarPathHeight = 16.89343;

}  // namespace

// Remove XML or DOCTYPE declarations and strip potentially dangerous
// constructs. Returns true when unsafe content was found.
bool sanitizeSvg(std::string &s) {
  bool unsafe = false;
  size_t p;
  while ((p = s.find("<?xml")) != std::string::npos) {
    size_t q = s.find("?>", p);
    if (q == std::string::npos)
      break;
    s.erase(p, q - p + 2);
    unsafe = true;
  }
  while ((p = s.find("<!DOCTYPE")) != std::string::npos) {
    size_t q = s.find('>', p);
    if (q == std::string::npos)
      break;
    s.erase(p, q - p + 1);
    unsafe = true;
  }

  auto replaceAndCheck = [&s, &unsafe](const std::regex &re,
                                       const std::string &rep) {
    std::string replaced = std::regex_replace(s, re, rep);
    if (replaced != s) {
      s = std::move(replaced);
      unsafe = true;
    }
  };

  replaceAndCheck(scriptRe, "");
  replaceAndCheck(onAttrRe, " ");
  replaceAndCheck(jsHrefRe, "");
  replaceAndCheck(styleTagRe, "");
  replaceAndCheck(styleAttrRe, " ");
  replaceAndCheck(dataUriAttrRe, " ");

  return unsafe;
}

// Compute the size of a landmark in pixels while respecting a maximum
// allowed width. The width cap is determined by measuring the rendered
// width of a placeholder string (ten underscores) at the configured
// station label size. If an icon or text would exceed this width, it is
// scaled down proportionally.
std::pair<double, double> getLandmarkSizePx(const Landmark &lm,
                                            const Config *cfg) {
  // Compute the maximum allowed width in pixels.
  double maxWidth = cfg->stationLabelSize * cfg->outputResolution * 0.6 *
                    10.0; // "__________"

  if (!lm.iconPath.empty()) {
    // lm.size is stored in map units, convert to pixels first
    double targetH = lm.size * cfg->outputResolution;
    std::ifstream iconFile(lm.iconPath);
    if (iconFile.good()) {
      std::stringstream buf;
      buf << iconFile.rdbuf();
      std::string svg = buf.str();

      auto extractAttr = [](const std::string &s,
                            const std::string &attr) -> double {
        size_t p = s.find(attr);
        if (p == std::string::npos)
          return std::numeric_limits<double>::quiet_NaN();
        p = s.find('"', p);
        if (p == std::string::npos)
          return std::numeric_limits<double>::quiet_NaN();
        size_t q = s.find('"', ++p);
        if (q == std::string::npos)
          return std::numeric_limits<double>::quiet_NaN();
        std::string val = s.substr(p, q - p);
        size_t e = 0;
        while (e < val.size() &&
               (std::isdigit(val[e]) || val[e] == '.' || val[e] == '-'))
          ++e;
        val = val.substr(0, e);
        try {
          return std::stod(val);
        } catch (...) {
          return std::numeric_limits<double>::quiet_NaN();
        }
      };

      double svgW = extractAttr(svg, "width");
      double svgH = extractAttr(svg, "height");
      if (!std::isnan(svgW) && !std::isnan(svgH) && svgH > 0) {
        double scale = targetH / svgH;
        double w = svgW * scale;
        double h = targetH;
        if (w > maxWidth) {
          double f = maxWidth / w;
          w = maxWidth;
          h *= f;
        }
        return {w, h};
      }
      size_t vbPos = svg.find("viewBox");
      if (vbPos != std::string::npos) {
        vbPos = svg.find('"', vbPos);
        if (vbPos != std::string::npos) {
          size_t vbEnd = svg.find('"', vbPos + 1);
          if (vbEnd != std::string::npos) {
            std::string vb = svg.substr(vbPos + 1, vbEnd - vbPos - 1);
            std::stringstream ss(vb);
            double minx, miny, vbW, vbH;
            if (ss >> minx >> miny >> vbW >> vbH && vbH > 0) {
              double scale = targetH / vbH;
              double w = vbW * scale;
              double h = targetH;
              if (w > maxWidth) {
                double f = maxWidth / w;
                w = maxWidth;
                h *= f;
              }
              return {w, h};
            }
          }
        }
      }
    }
    double w = targetH;
    double h = targetH;
    if (w > maxWidth) {
      double f = maxWidth / w;
      w = maxWidth;
      h *= f;
    }
    return {w, h};
  } else if (!lm.label.empty()) {
    // Desired label height is given directly in pixels.
    double h = lm.fontSize;
    // Use UTF-8 aware character counting for width estimation
    size_t cpCount = util::toWStr(lm.label).size();
    double w = cpCount * (h * 0.6);
    return {w, h};
  }
  // Fallback square size, again converting from map units to pixels
  double w = lm.size * cfg->outputResolution;
  double h = w;
  if (w > maxWidth) {
    double f = maxWidth / w;
    w = maxWidth;
    h *= f;
  }
  return {w, h};
}

namespace {

Landmark getAdjustedMeLandmark(const Config *cfg, double starPx,
                               bool badgeMode) {
  Landmark lm = cfg->meLandmark;
  if (badgeMode && !cfg->meLabelSizeExplicit) {
    Config defaults;
    double baseStarSize =
        cfg->meStarSizeExplicit ? defaults.meStarSize : cfg->meStarSize;
    double baseStarPx = baseStarSize * cfg->outputResolution;
    if (baseStarPx > 0.0) {
      double scaledStarPx = std::max(starPx, 0.0);
      double derivedFontSize =
          defaults.meLabelSize * (scaledStarPx / baseStarPx);
      lm.fontSize = derivedFontSize;
    } else {
      lm.fontSize = defaults.meLabelSize;
    }
  }
  return lm;
}

}  // namespace

// _____________________________________________________________________________
SvgRenderer::SvgRenderer(std::ostream *o, const Config *cfg)
    : _o(o), _w(o, true), _cfg(cfg) {}

// _____________________________________________________________________________
util::geo::Box<double> SvgRenderer::computeBgMapBBox() const {
  util::geo::Box<double> box;
  if (_cfg->bgMapPath.empty())
    return box;
  std::ifstream in(_cfg->bgMapPath);
  if (!in.good())
    return box;
  nlohmann::json j;
  try {
    in >> j;
  } catch (...) {
    return box;
  }
  if (!j.contains("features"))
    return box;
  std::function<void(const nlohmann::json &)> collect =
      [&](const nlohmann::json &coords) {
        if (!coords.is_array())
          return;
        if (!coords.empty() && coords[0].is_array()) {
          for (const auto &sub : coords)
            collect(sub);
        } else if (coords.size() >= 2 && coords[0].is_number() &&
                   coords[1].is_number()) {
          DPoint p(coords[0].get<double>(), coords[1].get<double>());
          if (!_cfg->bgMapWebmerc) {
            p = util::geo::latLngToWebMerc(p);
          }
          box = util::geo::extendBox(p, box);
        }
      };
  for (const auto &f : j["features"]) {
    if (!f.contains("geometry"))
      continue;
    const auto &geom = f["geometry"];
    if (!geom.contains("coordinates"))
      continue;
    collect(geom["coordinates"]);
  }
  return box;
}

// _____________________________________________________________________________
util::geo::DPoint SvgRenderer::findFreeLandmarkPosition(
    util::geo::DPoint base, double halfW, double halfH,
    const util::geo::Box<double> &renderBox,
    const std::vector<util::geo::Box<double>> &usedBoxes,
    double radius) const {
  util::geo::DPoint pos = base;
  double step = std::max(halfW, halfH);
  auto clampToRender = [&](const util::geo::DPoint &p) {
    double minX = renderBox.getLowerLeft().getX() + halfW;
    double maxX = renderBox.getUpperRight().getX() - halfW;
    double minY = renderBox.getLowerLeft().getY() + halfH;
    double maxY = renderBox.getUpperRight().getY() - halfH;
    double x = std::min(std::max(p.getX(), minX), maxX);
    double y = std::min(std::max(p.getY(), minY), maxY);
    return util::geo::DPoint(x, y);
  };

  for (int i = 0; i < _cfg->displacementIterations; ++i) {
    util::geo::Box<double> box(util::geo::DPoint(pos.getX() - halfW,
                                                pos.getY() - halfH),
                               util::geo::DPoint(pos.getX() + halfW,
                                                pos.getY() + halfH));
    bool overlap = false;
    double fx = 0.0, fy = 0.0;
    auto c1 = util::geo::centroid(box);
    for (const auto &b : usedBoxes) {
      if (util::geo::intersects(box, b)) {
        overlap = true;
        auto c2 = util::geo::centroid(b);
        double dx = c1.getX() - c2.getX();
        double dy = c1.getY() - c2.getY();
        double dist = std::sqrt(dx * dx + dy * dy);
        if (dist < 1e-6) {
          dx = 1.0;
          dy = 0.0;
          dist = 1.0;
        }
        fx += dx / dist;
        fy += dy / dist;
      }
    }
    if (!overlap) {
      break;
    }
    double norm = std::sqrt(fx * fx + fy * fy);
    if (norm > 0) {
      pos = util::geo::DPoint(pos.getX() + (fx / norm) * step,
                              pos.getY() + (fy / norm) * step);
      pos = clampToRender(pos);
      double dx = pos.getX() - base.getX();
      double dy = pos.getY() - base.getY();
      double dist = std::sqrt(dx * dx + dy * dy);
      if (dist > radius) {
        double f = radius / dist;
        pos = util::geo::DPoint(base.getX() + dx * f,
                                base.getY() + dy * f);
      }
    }
    step *= _cfg->displacementCooling;
  }
  return pos;
}

// _____________________________________________________________________________
void SvgRenderer::print(const RenderGraph &outG) {
  std::map<std::string, std::string> params;
  RenderParams rparams;
  _arrowHeads.clear();
  _meStationLabelVisual = Nullable<StationLabelVisual>();

  auto box = outG.getBBox();
  box = util::geo::pad(box, outG.getMaxLineNum() *
                                (_cfg->lineWidth + _cfg->lineSpacing));
  if (_cfg->geoLock) {
    box = util::geo::extendBox(_cfg->geoLockBox, box);
  }
  auto initialBox = box;

  if (_cfg->extendWithBgMap && !_cfg->bgMapPath.empty()) {
    auto bgBox = computeBgMapBBox();
    box = util::geo::extendBox(bgBox, box);
  }

  Labeller labeller(_cfg);
  std::vector<Landmark> acceptedLandmarks;
  for (const auto &lm : outG.getLandmarks()) {
    auto dims = ::getLandmarkSizePx(lm, _cfg);
    double halfW = (dims.first / _cfg->outputResolution) / 2.0;
    double halfH = (dims.second / _cfg->outputResolution) / 2.0;
    util::geo::Box<double> lmBox(
        DPoint(lm.coord.getX() - halfW, lm.coord.getY() - halfH),
        DPoint(lm.coord.getX() + halfW, lm.coord.getY() + halfH));

    if (!util::geo::intersects(lmBox, initialBox))
      continue;

    if (lm.label.empty()) {
      labeller.addLandmark(lmBox);
    }
    box = util::geo::extendBox(lmBox, box);
    acceptedLandmarks.push_back(lm);
  }
  if (_cfg->renderMe || _cfg->forceMeStar) {
    double starPx = _cfg->meStarSize * _cfg->outputResolution;
    bool highlightIntent =
        _cfg->highlightMeStationLabel && !_cfg->meStationId.empty();
    bool showLabel = _cfg->renderMeLabel || _cfg->meStationWithBg || highlightIntent;
    bool badgeMode = (_cfg->meStationWithBg || highlightIntent) && showLabel;
    Landmark meLm = getAdjustedMeLandmark(_cfg, starPx, badgeMode);
    if (meLm.label.empty() && !_cfg->meStationLabel.empty()) {
      meLm.label = _cfg->meStationLabel;
    }
    double labelWpx = 0.0;
    double labelHpx = 0.0;
    if (showLabel) {
      auto dims = ::getLandmarkSizePx(meLm, _cfg);
      labelWpx = dims.first;
      labelHpx = dims.second;
    }
    double starGap = showLabel ? starPx * 0.2 : 0.0;
    double boxWpx = 0.0;
    double boxHpx = 0.0;
    if (badgeMode) {
      double textHeightForPadding =
          labelHpx > 0.0 ? labelHpx : starPx;
      double padX = textHeightForPadding * 0.6;
      double padTop = textHeightForPadding * 0.28;
      double padBottom = textHeightForPadding * 0.12;
      double contentHeightPx = std::max(starPx, textHeightForPadding);
      boxWpx = padX * 2.0 + starPx + starGap + labelWpx;
      boxHpx = padTop + padBottom + contentHeightPx;
    } else {
      boxWpx = std::max(labelWpx, starPx);
      boxHpx = starPx + starGap + labelHpx;
    }
    double halfW = (boxWpx / _cfg->outputResolution) / 2.0;
    double halfH = (boxHpx / _cfg->outputResolution) / 2.0;
    util::geo::Box<double> lmBox(DPoint(_cfg->meLandmark.coord.getX() - halfW,
                                        _cfg->meLandmark.coord.getY() - halfH),
                                 DPoint(_cfg->meLandmark.coord.getX() + halfW,
                                        _cfg->meLandmark.coord.getY() + halfH));
    box = util::geo::extendBox(lmBox, box);
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

  renderBackground(rparams);

  // Landmarks collected above already lie within the padded network box.
  LOGTO(DEBUG, std::cerr) << "[DEBUG] acceptedLandmarks.size() = "
                          << acceptedLandmarks.size() << '\n';

  LOGTO(DEBUG, std::cerr) << "Rendering nodes...";
  for (auto n : outG.getNds()) {
    if (_cfg->renderNodeConnections) {
      renderNodeConnections(outG, n, rparams);
    }
  }

  LOGTO(DEBUG, std::cerr) << "Writing edges...";
  renderDelegates(outG, rparams);

  for (const auto &ah : _arrowHeads) {
    if (ah.pts.empty())
      continue;
    std::stringstream d;
    auto pt = ah.pts.begin();
    double x = (pt->getX() - rparams.xOff) * _cfg->outputResolution;
    double y =
        rparams.height - (pt->getY() - rparams.yOff) * _cfg->outputResolution;
    d << "M" << x << " " << y;
    for (++pt; pt != ah.pts.end(); ++pt) {
      x = (pt->getX() - rparams.xOff) * _cfg->outputResolution;
      y = rparams.height - (pt->getY() - rparams.yOff) * _cfg->outputResolution;
      d << " L" << x << " " << y;
    }
    d << " Z";
    params.clear();
    params["d"] = d.str();
    params["style"] = "fill:white;stroke:none";
    _w.openTag("path", params);
    _w.closeTag();
  }

  LOGTO(DEBUG, std::cerr) << "Writing nodes...";
  outputNodes(outG, rparams);
  if (_cfg->renderNodeFronts) {
    renderNodeFronts(outG, rparams);
  }

  // Render landmarks on top of edges and nodes but below labels.
  LOGTO(DEBUG, std::cerr) << "Writing landmarks...";
  renderLandmarks(outG, acceptedLandmarks, rparams);

  LOGTO(DEBUG, std::cerr) << "Writing labels...";
  if (_cfg->renderLabels) {
    renderLineLabels(labeller, rparams);

    if (_cfg->renderRouteLabels) {
      renderTerminusLabels(outG, labeller, rparams);
    }

    renderStationLabels(labeller, rparams);
  }

  if (_cfg->renderMe || _cfg->forceMeStar) {
    renderMe(outG, labeller, rparams);
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
      bool isTerminus = RenderGraph::isTerminus(n);
      std::string stroke = "black";
      std::string fill = "white";
      if (_cfg->highlightTerminals && isTerminus) {
        stroke = _cfg->terminusHighlightStroke;
        fill = _cfg->terminusHighlightFill;
      }
      params["stroke"] = stroke;
      params["stroke-width"] =
          util::toString((_cfg->lineWidth / 2) * _cfg->outputResolution);
      params["fill"] = fill;
      const auto &st = n->pl().stops().front();
      if (st.labelDeg != std::numeric_limits<size_t>::max()) {
        params["labelDeg"] = util::toString(st.labelDeg);
      }

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
void SvgRenderer::renderBackground(const RenderParams &rparams) {
  if (_cfg->bgMapPath.empty())
    return;
  std::ifstream in(_cfg->bgMapPath);
  if (!in.good())
    return;
  nlohmann::json j;
  try {
    in >> j;
  } catch (...) {
    return;
  }
  if (!j.contains("features"))
    return;
  Params baseParams;
  baseParams["class"] = "bg-map";
  for (const auto &f : j["features"]) {
    if (!f.contains("geometry"))
      continue;
    const auto &geom = f["geometry"];
    if (!geom.contains("type") || !geom.contains("coordinates"))
      continue;

    std::map<std::string, std::string> params = baseParams;
    std::string stroke = "#ccc";
    double strokeWidth = _cfg->lineWidth;
    std::string fill = "none";
    double opacity = _cfg->bgMapOpacity;

    if (f.contains("properties") && f["properties"].is_object()) {
      const auto &props = f["properties"];
      auto getStr = [](const nlohmann::json &v) {
        return v.is_string() ? v.get<std::string>()
                             : std::to_string(v.get<double>());
      };
      auto getDouble = [](const nlohmann::json &v) {
        return v.is_number() ? v.get<double>()
                             : atof(v.get<std::string>().c_str());
      };
      if (props.contains("stroke") && !props["stroke"].is_null())
        stroke = getStr(props["stroke"]);
      if (props.contains("stroke-width") && !props["stroke-width"].is_null())
        strokeWidth = getDouble(props["stroke-width"]);
      if (props.contains("fill") && !props["fill"].is_null())
        fill = getStr(props["fill"]);
      if (props.contains("opacity") && !props["opacity"].is_null())
        opacity = getDouble(props["opacity"]);
      if (props.contains("class") && !props["class"].is_null())
        params["class"] += " " + getStr(props["class"]);
    }

    std::stringstream style;
    style << "fill:" << fill << ";stroke:" << stroke
          << ";stroke-width:" << strokeWidth * _cfg->outputResolution
          << ";stroke-opacity:" << opacity << ";fill-opacity:" << opacity;
    params["style"] = style.str();

    std::string type = geom["type"].get<std::string>();
    if (type == "LineString") {
      PolyLine<double> pl;
      for (const auto &c : geom["coordinates"]) {
        if (c.size() < 2)
          continue;
        DPoint p(c[0].get<double>(), c[1].get<double>());
        if (!_cfg->bgMapWebmerc) {
          p = util::geo::latLngToWebMerc(p);
        }
        pl << p;
      }
      if (pl.getLine().size() > 1)
        printLine(pl, params, rparams);
    } else if (type == "MultiLineString") {
      for (const auto &line : geom["coordinates"]) {
        PolyLine<double> pl;
        for (const auto &c : line) {
          if (c.size() < 2)
            continue;
          DPoint p(c[0].get<double>(), c[1].get<double>());
          if (!_cfg->bgMapWebmerc) {
            p = util::geo::latLngToWebMerc(p);
          }
          pl << p;
        }
        if (pl.getLine().size() > 1)
          printLine(pl, params, rparams);
      }
    } else if (type == "Polygon") {
      const auto &coords = geom["coordinates"];
      if (!coords.empty()) {
        util::geo::Line<double> outer;
        for (const auto &c : coords[0]) {
          if (c.size() < 2)
            continue;
          DPoint p(c[0].get<double>(), c[1].get<double>());
          if (!_cfg->bgMapWebmerc) {
            p = util::geo::latLngToWebMerc(p);
          }
          outer.push_back(p);
        }
        if (outer.size() > 2) {
          util::geo::Polygon<double> poly(outer);
          printPolygon(poly, params, rparams);
        }
      }
    } else if (type == "MultiPolygon") {
      for (const auto &polyCoords : geom["coordinates"]) {
        if (polyCoords.empty())
          continue;
        util::geo::Line<double> outer;
        for (const auto &c : polyCoords[0]) {
          if (c.size() < 2)
            continue;
          DPoint p(c[0].get<double>(), c[1].get<double>());
          if (!_cfg->bgMapWebmerc) {
            p = util::geo::latLngToWebMerc(p);
          }
          outer.push_back(p);
        }
        if (outer.size() > 2) {
          util::geo::Polygon<double> poly(outer);
          printPolygon(poly, params, rparams);
        }
      }
    }
  }
}

// _____________________________________________________________________________
void SvgRenderer::renderLandmarks(const RenderGraph &g,
                                  const std::vector<Landmark> &landmarks,
                                  const RenderParams &rparams) {
  std::map<std::string, std::string> iconIds;
  size_t id = 0;

  auto logBox = [](const char *name, const util::geo::Box<double> &box) {
    LOGTO(DEBUG, std::cerr)
        << name << " ll=(" << box.getLowerLeft().getX() << ", "
        << box.getLowerLeft().getY() << ") ur=(" << box.getUpperRight().getX()
        << ", " << box.getUpperRight().getY() << ")";
  };

  // collect existing geometry bounding boxes (nodes and edges) to avoid
  // drawing landmarks on top of them
  std::vector<util::geo::Box<double>> usedBoxes;
  std::set<const shared::linegraph::LineEdge *> processedEdges;
  for (auto n : g.getNds()) {
    for (const auto &poly : g.getStopGeoms(n, _cfg->tightStations, 32)) {
      usedBoxes.push_back(util::geo::extendBox(poly, util::geo::Box<double>()));
    }
    for (auto e : n->getAdjList()) {
      if (processedEdges.insert(e).second) {
        util::geo::Box<double> b = util::geo::extendBox(
            e->pl().getPolyline().getLine(), util::geo::Box<double>());
        b = util::geo::pad(b, g.getTotalWidth(e) / 2.0);
        usedBoxes.push_back(b);
      }
    }
  }

  util::geo::Box<double> renderBox(
      DPoint(rparams.xOff, rparams.yOff),
      DPoint(rparams.xOff + rparams.width / _cfg->outputResolution,
             rparams.yOff + rparams.height / _cfg->outputResolution));

  logBox("renderBox", renderBox);

  _w.openTag("defs");
  _w.writeText("");

  for (const auto &lm : landmarks) {
    LOGTO(DEBUG, std::cerr)
        << "Adding landmark "
        << (!lm.iconPath.empty() ? "icon=" + lm.iconPath : "label=" + lm.label)
        << " color=" << lm.color << " size=" << lm.size
        << " fontSize=" << lm.fontSize << " coord=(" << lm.coord.getX() << ","
        << lm.coord.getY() << ")";

    if (lm.iconPath.empty())
      continue;
    auto it = iconIds.find(lm.iconPath);
    if (it == iconIds.end()) {
      std::ifstream iconFile(lm.iconPath);
      if (!iconFile.good()) {
        LOGTO(DEBUG, std::cerr)
            << "Cannot read icon file \"" << lm.iconPath << "\"";
        continue;
      }
      std::stringstream buf;
      buf << iconFile.rdbuf();
      std::string svg = buf.str();
      std::string idStr = "lmk" + util::toString(id++);
      iconIds[lm.iconPath] = idStr;
      // Remove XML or DOCTYPE declarations and strip potentially dangerous
      // constructs. Returns true when unsafe content was found.
      bool unsafe = sanitizeSvg(svg);

      size_t pos = svg.find("<svg");
      if (pos != std::string::npos) {
        svg = svg.substr(pos);
        unsafe |= sanitizeSvg(svg);
        size_t end = svg.find('>');
        if (end != std::string::npos) {
          svg.insert(end, " id=\"" + idStr + "\"");
        }
        size_t close = svg.rfind("</svg>");
        if (close != std::string::npos) {
          svg = svg.substr(0, close + 6);
        }
        *_o << svg;
      } else {
        unsafe |= sanitizeSvg(svg);
        *_o << "<svg id=\"" << idStr
            << "\" xmlns=\"http://www.w3.org/2000/svg\">" << svg << "</svg>";
      }

      if (unsafe) {
        LOGTO(DEBUG, std::cerr)
            << "Unsafe SVG content removed from icon '" << lm.iconPath << "'";
      }
    }
  }
  _w.closeTag();

  _w.openTag("g");
  for (const auto &lm : landmarks) {
    auto dimsPx = ::getLandmarkSizePx(lm, _cfg);
    double wPx = dimsPx.first;
    double hPx = !lm.label.empty() ? lm.fontSize : dimsPx.second;
    double fontSizePx = lm.fontSize;
    double halfW = (wPx / _cfg->outputResolution) / 2.0;
    double halfH = (hPx / _cfg->outputResolution) / 2.0;
    util::geo::DPoint coord = lm.coord;
    util::geo::Box<double> lmBox(
        DPoint(coord.getX() - halfW, coord.getY() - halfH),
        DPoint(coord.getX() + halfW, coord.getY() + halfH));

    LOGTO(DEBUG, std::cerr)
        << "Landmark "
        << (!lm.iconPath.empty() ? "icon"
                                 : (!lm.label.empty() ? "label" : "unknown"))
        << " at (" << coord.getX() << ", " << coord.getY() << ") dimsPx=("
        << dimsPx.first << ", " << dimsPx.second << ")";

    if (!util::geo::contains(lmBox, renderBox)) {
      LOGTO(DEBUG, std::cerr) << "Skipping landmark at (" << coord.getX()
                              << ", " << coord.getY() << ") outside render box";
      continue;
    }

    bool overlaps = false;
    if (!_cfg->renderOverlappingLandmarks && lm.label.empty()) {
      for (const auto &b : usedBoxes) {
        if (util::geo::intersects(lmBox, b)) {
          overlaps = true;
          break;
        }
      }
      if (overlaps && !lm.iconPath.empty()) {
        double searchRadius =
            _cfg->landmarkSearchRadius * std::max(halfW, halfH);
        util::geo::DPoint cand = findFreeLandmarkPosition(
            coord, halfW, halfH, renderBox, usedBoxes, searchRadius);
        util::geo::Box<double> box(
            DPoint(cand.getX() - halfW, cand.getY() - halfH),
            DPoint(cand.getX() + halfW, cand.getY() + halfH));
        coord = cand;
        lmBox = box;
        overlaps = false;
        for (const auto &b : usedBoxes) {
          if (util::geo::intersects(lmBox, b)) {
            overlaps = true;
            break;
          }
        }
      }
      if (overlaps && lm.label.empty() && lm.iconPath.empty()) {
        LOGTO(DEBUG, std::cerr) << "Skipping landmark at (" << coord.getX()
                                << ", " << coord.getY() << ") due to overlap";
        continue;
      }
    }

    if (!lm.iconPath.empty()) {
      auto it = iconIds.find(lm.iconPath);

      double x = (coord.getX() - rparams.xOff) * _cfg->outputResolution -
                 dimsPx.first / 2.0;
      double y = rparams.height -
                 (coord.getY() - rparams.yOff) * _cfg->outputResolution -
                 dimsPx.second / 2.0;

      if (it == iconIds.end()) {
        LOGTO(DEBUG, std::cerr) << "Missing icon '" << lm.iconPath
                                << "', drawing placeholder rectangle";
      } else {
        std::map<std::string, std::string> attrs;
        attrs["xlink:href"] = "#" + it->second;
        attrs["x"] = util::toString(x);
        attrs["y"] = util::toString(y);
        attrs["width"] = util::toString(dimsPx.first);
        attrs["height"] = util::toString(dimsPx.second);
        attrs["class"] = util::toString(lm.cssClass);
        _w.openTag("use", attrs);
        _w.closeTag();
      }
      usedBoxes.push_back(lmBox);
    }
    if (!lm.label.empty()) {
      double x =
          (coord.getX() - rparams.xOff) * _cfg->outputResolution - wPx / 2.0;
      double y = rparams.height -
                 (coord.getY() - rparams.yOff) * _cfg->outputResolution -
                 fontSizePx / 2.0;

      std::map<std::string, std::string> params;
      params["x"] = util::toString(x + wPx / 2.0);
      params["y"] = util::toString(y + fontSizePx / 2.0);
      params["font-size"] = util::toString(fontSizePx);
      params["font-weight"] = "bold";
      params["text-anchor"] = "middle";
      params["fill"] = lm.color;
      params["font-family"] = "TT Norms Pro";
      params["class"] = util::toString(lm.cssClass);
      params["opacity"] = util::toString(lm.opacity);
      _w.openTag("text", params);
      _w.writeText(lm.label);
      _w.closeTag();
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::renderMe(const RenderGraph &g, Labeller &labeller,
                           const RenderParams &rparams) {
  double starPx = _cfg->meStarSize * _cfg->outputResolution;
  bool highlightAvailable = false;
  StationLabelVisual highlightInfo;
  if (_cfg->highlightMeStationLabel && !_cfg->meStationId.empty() &&
      !_meStationLabelVisual.isNull()) {
    const StationLabelVisual &stored = _meStationLabelVisual.get();
    if (stored.label && !stored.pathId.empty()) {
      highlightAvailable = true;
      highlightInfo = stored;
    }
  }

  if (highlightAvailable) {
    const StationLabel *label = highlightInfo.label;
    auto isValidBox = [](const util::geo::Box<double> &b) {
      return b.getLowerLeft().getX() <= b.getUpperRight().getX() &&
             b.getLowerLeft().getY() <= b.getUpperRight().getY();
    };
    double res = _cfg->outputResolution;

    auto textPath = label->geom;
    double pathAng = util::geo::angBetween(textPath.front(), textPath.back());
    if ((fabs(pathAng) < (3 * M_PI / 2)) && (fabs(pathAng) > (M_PI / 2))) {
      textPath.reverse();
    }

    bool sampleEnd = highlightInfo.startOffset == "100%" ||
                     highlightInfo.textAnchor == "end";
    double sampleSpan = 0.05;
    double tStart = sampleEnd ? std::max(0.0, 1.0 - sampleSpan) : 0.0;
    double tEnd = sampleEnd ? 1.0 : std::min(1.0, sampleSpan);
    if (std::fabs(tEnd - tStart) < 1e-6) {
      if (sampleEnd) {
        tStart = std::max(0.0, 1.0 - 1e-3);
      } else {
        tEnd = std::min(1.0, 1e-3);
      }
    }

    util::geo::LinePoint<double> anchorPt =
        textPath.getPointAt(sampleEnd ? 1.0 : 0.0);
    auto slope = textPath.getSlopeBetween(tStart, tEnd);
    double dirX = slope.first;
    double dirY = slope.second;
    double dirNorm = std::hypot(dirX, dirY);
    if (!(dirNorm > 0)) {
      dirX = 1.0;
      dirY = 0.0;
      dirNorm = 1.0;
    }
    dirX /= dirNorm;
    dirY /= dirNorm;

    double dirScreenX = dirX;
    double dirScreenY = -dirY;
    double dirScreenNorm = std::hypot(dirScreenX, dirScreenY);
    if (!(dirScreenNorm > 0)) {
      dirScreenX = 1.0;
      dirScreenY = 0.0;
      dirScreenNorm = 1.0;
    }
    dirScreenX /= dirScreenNorm;
    dirScreenY /= dirScreenNorm;

    double angleRad = std::atan2(dirScreenY, dirScreenX);
    double angleDeg = angleRad * 180.0 / M_PI;
    if (!std::isfinite(angleDeg)) {
      angleDeg = 0.0;
    }

    double anchorXPx = (anchorPt.p.getX() - rparams.xOff) * res;
    double anchorYPx =
        rparams.height - (anchorPt.p.getY() - rparams.yOff) * res;

    struct ScreenPoint {
      double x;
      double y;
    };

    std::vector<ScreenPoint> extentPoints;
    extentPoints.reserve(16);
    for (const auto &line : label->band) {
      for (const auto &pt : line) {
        double px = (pt.getX() - rparams.xOff) * res;
        double py = rparams.height - (pt.getY() - rparams.yOff) * res;
        if (std::isfinite(px) && std::isfinite(py)) {
          extentPoints.push_back({px, py});
        }
      }
    }

    size_t cpCount = util::toWStr(label->s.name).size();
    double estimatedWidth = cpCount * (highlightInfo.fontSizePx * 0.6);
    if (!(estimatedWidth > 0)) {
      estimatedWidth = highlightInfo.fontSizePx;
    }
    double estimatedHalfWidth = estimatedWidth / 2.0;
    double estimatedHalfHeight = highlightInfo.fontSizePx / 2.0;

    bool usedFallbackGeometry = false;
    auto addFallbackCorners = [&]() {
      util::geo::Box<double> labelBox =
          util::geo::extendBox(label->band, util::geo::Box<double>());
      double labelLeftPx = 0.0;
      double labelRightPx = 0.0;
      double labelTopPx = 0.0;
      double labelBottomPx = 0.0;
      if (isValidBox(labelBox)) {
        labelLeftPx = (labelBox.getLowerLeft().getX() - rparams.xOff) * res;
        labelRightPx = (labelBox.getUpperRight().getX() - rparams.xOff) * res;
        labelTopPx = rparams.height -
                     (labelBox.getUpperRight().getY() - rparams.yOff) * res;
        labelBottomPx = rparams.height -
                        (labelBox.getLowerLeft().getY() - rparams.yOff) * res;
      } else {
        util::geo::LinePoint<double> mid = label->geom.getPointAt(0.5);
        double centerXPx = (mid.p.getX() - rparams.xOff) * res;
        double centerYPx =
            rparams.height - (mid.p.getY() - rparams.yOff) * res;
        labelLeftPx = centerXPx - estimatedHalfWidth;
        labelRightPx = centerXPx + estimatedHalfWidth;
        labelTopPx = centerYPx - estimatedHalfHeight;
        labelBottomPx = centerYPx + estimatedHalfHeight;
      }
      std::array<ScreenPoint, 4> corners = {{{labelLeftPx, labelTopPx},
                                             {labelRightPx, labelTopPx},
                                             {labelRightPx, labelBottomPx},
                                             {labelLeftPx, labelBottomPx}}};
      extentPoints.insert(extentPoints.end(), corners.begin(), corners.end());
    };

    if (extentPoints.empty()) {
      addFallbackCorners();
      usedFallbackGeometry = true;
    }

    double textAlongMin = 0.0;
    double textAlongMax = 0.0;
    double textPerpMin = 0.0;
    double textPerpMax = 0.0;
    auto computeExtents = [&]() {
      double alongMin = std::numeric_limits<double>::infinity();
      double alongMax = -std::numeric_limits<double>::infinity();
      double perpMin = std::numeric_limits<double>::infinity();
      double perpMax = -std::numeric_limits<double>::infinity();
      for (const auto &pt : extentPoints) {
        double relX = pt.x - anchorXPx;
        double relY = pt.y - anchorYPx;
        double along = relX * dirScreenX + relY * dirScreenY;
        double perp = relX * (-dirScreenY) + relY * dirScreenX;
        alongMin = std::min(alongMin, along);
        alongMax = std::max(alongMax, along);
        perpMin = std::min(perpMin, perp);
        perpMax = std::max(perpMax, perp);
      }
      if (!std::isfinite(alongMin) || !std::isfinite(alongMax) ||
          !std::isfinite(perpMin) || !std::isfinite(perpMax)) {
        return false;
      }
      textAlongMin = alongMin;
      textAlongMax = alongMax;
      textPerpMin = perpMin;
      textPerpMax = perpMax;
      return true;
    };

    bool extentsValid = computeExtents();
    if ((!extentsValid || !(textAlongMax > textAlongMin) ||
         !(textPerpMax > textPerpMin)) &&
        !usedFallbackGeometry) {
      extentPoints.clear();
      addFallbackCorners();
      usedFallbackGeometry = true;
      extentsValid = computeExtents();
    }

    if (!extentsValid || !(textAlongMax > textAlongMin) ||
        !(textPerpMax > textPerpMin)) {
      textAlongMin = -estimatedHalfWidth;
      textAlongMax = estimatedHalfWidth;
      textPerpMin = -estimatedHalfHeight;
      textPerpMax = estimatedHalfHeight;
    }

    double textWidthAlong = textAlongMax - textAlongMin;
    if (!(textWidthAlong > 0)) {
      textWidthAlong = std::max(estimatedWidth, highlightInfo.fontSizePx);
      textAlongMin = -textWidthAlong / 2.0;
      textAlongMax = textAlongMin + textWidthAlong;
    }

    double textPerpSpan = textPerpMax - textPerpMin;
    if (!(textPerpSpan > 0)) {
      textPerpSpan = highlightInfo.fontSizePx;
      if (!(textPerpSpan > 0)) {
        textPerpSpan = 1.0;
      }
      textPerpMin = -textPerpSpan / 2.0;
      textPerpMax = textPerpSpan / 2.0;
    }
    double textPerpCenter = (textPerpMin + textPerpMax) / 2.0;

    double starSizePx = starPx;
    if (textPerpSpan > 0.0) {
      starSizePx = std::min(starSizePx, textPerpSpan);
    }
    double starGapPx = starSizePx * kBadgeStarGapFactor;
    double textHeightForPadding = std::max(textPerpSpan, 1.0);
    double padX = textHeightForPadding * kBadgePadXFactor;
    double padTop = textHeightForPadding * kBadgePadTopFactor;
    double padBottom =
        textHeightForPadding * kBadgePadBottomFactor;
    double contentHeightPx = textHeightForPadding;
    double rectHeight = padTop + padBottom + contentHeightPx;

    double rectWidth = padX * 2.0 + starSizePx + starGapPx + textWidthAlong;
    double rectStartAlong = textAlongMin - padX - starGapPx - starSizePx;
    double rectPerpTop = textPerpCenter - contentHeightPx / 2.0 - padTop;
    double badgeCenterPerp = textPerpCenter + (padBottom - padTop) / 2.0;
    double starCenterAlong = textAlongMin - starGapPx - starSizePx / 2.0;
    double starCenterPerp = badgeCenterPerp;

    std::stringstream transform;
    transform << "translate(" << anchorXPx << " " << anchorYPx
              << ") rotate(" << angleDeg << ")";
    std::map<std::string, std::string> groupAttrs;
    groupAttrs["transform"] = transform.str();
    _w.openTag("g", groupAttrs);

    std::map<std::string, std::string> rectAttrs;
    rectAttrs["x"] = util::toString(rectStartAlong);
    rectAttrs["y"] = util::toString(rectPerpTop);
    rectAttrs["width"] = util::toString(rectWidth);
    rectAttrs["height"] = util::toString(rectHeight);
    double radius = std::min(rectHeight / 2.0, textHeightForPadding);
    rectAttrs["rx"] = util::toString(radius);
    rectAttrs["ry"] = util::toString(radius);
    rectAttrs["fill"] = _cfg->meStationBgFill;
    rectAttrs["stroke"] = _cfg->meStationBgStroke;
    _w.openTag("rect", rectAttrs);
    _w.closeTag();

    double scaleX = starSizePx / kBadgeStarPathWidth;
    double scaleY = starSizePx / kBadgeStarPathHeight;
    std::stringstream starTransform;
    starTransform << "translate(" << starCenterAlong << ' ' << starCenterPerp
                  << ") scale(" << scaleX << ' ' << scaleY << ") translate("
                  << -kBadgeStarPathCenterX << ' ' << -kBadgeStarPathCenterY
                  << ")";
    std::map<std::string, std::string> starAttrs;
    starAttrs["d"] = kBadgeStarPath;
    starAttrs["transform"] = starTransform.str();
    starAttrs["fill-rule"] = "evenodd";
    starAttrs["clip-rule"] = "evenodd";
    starAttrs["fill"] = _cfg->meStationFill;
    starAttrs["stroke"] = _cfg->meStationBorder;
    _w.openTag("path", starAttrs);
    _w.closeTag();

    std::map<std::string, std::string> textAttrs;
    double textCenterAlong = (textAlongMin + textAlongMax) / 2.0;
    textAttrs["class"] = "station-label";
    textAttrs["x"] = util::toString(textCenterAlong);
    textAttrs["y"] = util::toString(badgeCenterPerp);
    textAttrs["text-anchor"] = "middle";
    textAttrs["dominant-baseline"] = "middle";
    textAttrs["alignment-baseline"] = "middle";
    textAttrs["fill"] = _cfg->meStationTextColor;
    textAttrs["font-family"] = "TT Norms Pro";
    textAttrs["font-size"] =
        util::toString(highlightInfo.fontSizePx > 0.0
                           ? highlightInfo.fontSizePx
                           : textHeightForPadding);
    textAttrs["font-weight"] = highlightInfo.bold ? "bold" : "normal";
    _w.openTag("text", textAttrs);
    _w.writeText(label->s.name);
    _w.closeTag();
    _w.closeTag();

    if (_cfg->outputResolution > 0.0) {
      double res = _cfg->outputResolution;
      double alongStartMu = rectStartAlong / res;
      double alongEndMu = (rectStartAlong + rectWidth) / res;
      double perpTopMu = rectPerpTop / res;
      double perpBottomMu = (rectPerpTop + rectHeight) / res;
      double perpX = -dirY;
      double perpY = dirX;
      double anchorX = anchorPt.p.getX();
      double anchorY = anchorPt.p.getY();
      double minX = std::numeric_limits<double>::infinity();
      double minY = std::numeric_limits<double>::infinity();
      double maxX = -std::numeric_limits<double>::infinity();
      double maxY = -std::numeric_limits<double>::infinity();
      std::array<std::pair<double, double>, 4> localCorners = {
          {{alongStartMu, perpTopMu},
           {alongEndMu, perpTopMu},
           {alongEndMu, perpBottomMu},
           {alongStartMu, perpBottomMu}}};
      for (const auto &corner : localCorners) {
        double mapX = anchorX + corner.first * dirX + corner.second * perpX;
        double mapY = anchorY + corner.first * dirY + corner.second * perpY;
        minX = std::min(minX, mapX);
        minY = std::min(minY, mapY);
        maxX = std::max(maxX, mapX);
        maxY = std::max(maxY, mapY);
      }
      if (std::isfinite(minX) && std::isfinite(minY) && std::isfinite(maxX) &&
          std::isfinite(maxY) && minX <= maxX && minY <= maxY) {
        util::geo::Box<double> badgeBox(util::geo::DPoint(minX, minY),
                                        util::geo::DPoint(maxX, maxY));
        labeller.addLandmark(badgeBox, label);
      }
    }
    return;
  }

  bool showLabel =
      _cfg->renderMeLabel ||
      (_cfg->meStationWithBg &&
       (!_cfg->highlightMeStationLabel || _meStationLabelVisual.isNull()));
  bool badgeMode = _cfg->meStationWithBg && showLabel;
  Landmark lm = getAdjustedMeLandmark(_cfg, starPx, badgeMode);
  if (lm.label.empty() && !_cfg->meStationLabel.empty()) {
    lm.label = _cfg->meStationLabel;
  }
  std::pair<double, double> dims = {0.0, 0.0};
  if (showLabel) {
    dims = ::getLandmarkSizePx(lm, _cfg);
  }
  double labelWidthPx = dims.first;
  double labelHeightPx = dims.second;
  double starRenderSize = starPx;
  if (badgeMode && showLabel && labelHeightPx > 0.0) {
    starRenderSize = std::min(starRenderSize, labelHeightPx);
  }
  double starGapPx = showLabel ? starRenderSize * kBadgeStarGapFactor : 0.0;
  double textHeightForPadding =
      (showLabel && labelHeightPx > 0.0) ? labelHeightPx : starRenderSize;
  double padX = badgeMode ? textHeightForPadding * kBadgePadXFactor : 0.0;
  double padTop = badgeMode ? textHeightForPadding * kBadgePadTopFactor : 0.0;
  double padBottom =
      badgeMode ? textHeightForPadding * kBadgePadBottomFactor : 0.0;
  double boxWpx = 0.0;
  double boxHpx = 0.0;
  double contentHeightPx = 0.0;
  if (badgeMode) {
    contentHeightPx = textHeightForPadding;
    boxWpx = padX * 2.0 + starRenderSize + starGapPx + labelWidthPx;
    boxHpx = padTop + padBottom + contentHeightPx;
  } else {
    boxWpx = std::max(labelWidthPx, starRenderSize);
    boxHpx = labelHeightPx + starRenderSize + starGapPx;
  }
  double halfW = (boxWpx / _cfg->outputResolution) / 2.0;
  double halfH = (boxHpx / _cfg->outputResolution) / 2.0;
  util::geo::DPoint base = lm.coord;
  auto makeLandmarkBox = [&](const util::geo::DPoint &center) {
    return util::geo::Box<double>(
        util::geo::DPoint(center.getX() - halfW, center.getY() - halfH),
        util::geo::DPoint(center.getX() + halfW, center.getY() + halfH));
  };
  util::geo::Box<double> box = makeLandmarkBox(base);
  labeller.addLandmark(box);
  util::geo::DPoint placed = base;

  double x = (placed.getX() - rparams.xOff) * _cfg->outputResolution;
  double y =
      rparams.height - (placed.getY() - rparams.yOff) * _cfg->outputResolution;
  double boxLeftPx = x - boxWpx / 2.0;
  double boxTopPx = y - boxHpx / 2.0;

  if (badgeMode) {
    std::map<std::string, std::string> rectAttrs;
    rectAttrs["x"] = util::toString(boxLeftPx);
    rectAttrs["y"] = util::toString(boxTopPx);
    rectAttrs["width"] = util::toString(boxWpx);
    rectAttrs["height"] = util::toString(boxHpx);
    double radius = std::min(boxHpx / 2.0, textHeightForPadding);
    rectAttrs["rx"] = util::toString(radius);
    rectAttrs["ry"] = util::toString(radius);
    rectAttrs["fill"] = _cfg->meStationBgFill;
    rectAttrs["stroke"] = _cfg->meStationBgStroke;
    _w.openTag("rect", rectAttrs);
    _w.closeTag();
  }

  double starCx = badgeMode ? boxLeftPx + padX + starRenderSize / 2.0 : x;
  double badgeCenterY = boxTopPx + boxHpx / 2.0;
  double starCy = badgeMode ? badgeCenterY
                            : (showLabel ? y - starGapPx - starRenderSize / 2.0 : y);
  double scaleX = starRenderSize / kBadgeStarPathWidth;
  double scaleY = starRenderSize / kBadgeStarPathHeight;
  std::stringstream starTransform;
  starTransform << "translate(" << starCx << ' ' << starCy << ") scale("
                << scaleX << ' ' << scaleY << ") translate("
                << -kBadgeStarPathCenterX << ' ' << -kBadgeStarPathCenterY
                << ")";
  std::map<std::string, std::string> attrs;
  attrs["d"] = kBadgeStarPath;
  attrs["transform"] = starTransform.str();
  attrs["fill-rule"] = "evenodd";
  attrs["clip-rule"] = "evenodd";
  attrs["fill"] = _cfg->meStationFill;
  attrs["stroke"] = _cfg->meStationBorder;
  _w.openTag("path", attrs);
  _w.closeTag();

  if (showLabel) {
    std::map<std::string, std::string> params;
    if (badgeMode) {
      double textAreaLeft = boxLeftPx + padX + starRenderSize + starGapPx;
      double textX = textAreaLeft + labelWidthPx / 2.0;
      params["x"] = util::toString(textX);
      params["y"] = util::toString(badgeCenterY);
      params["text-anchor"] = "middle";
      params["dominant-baseline"] = "middle";
      params["alignment-baseline"] = "middle";
      params["fill"] = _cfg->meStationTextColor;
    } else {
      params["x"] = util::toString(x);
      params["y"] = util::toString(y);
      params["text-anchor"] = "middle";
      params["fill"] = _cfg->meStationFill;
    }
    params["font-size"] =
        util::toString(std::max(labelHeightPx, 1.0));
    params["font-family"] = "TT Norms Pro";

    _w.openTag("text", params);
    _w.writeText(lm.label);
    _w.closeTag();
  }
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
      double ang1 = fabs(util::geo::angBetween(a.geom.front(), a.geom.back()));
      double ang2 = fabs(util::geo::angBetween(b.geom.front(), b.geom.back()));

      return ang1 > THRESHOLD && ang2 > THRESHOLD;
    }
  }

  if (a.to.edge == b.from.edge && a.from.edge == b.to.edge) {
    if ((aSlotFrom - bSlotTo == 1 && bSlotFrom - aSlotTo == 1) ||
        (bSlotTo - aSlotFrom == 1 && aSlotTo - bSlotFrom == 1)) {
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

    const LineNode *nd = n;
    if (ref.from.edge && ref.to.edge) {
      if (const auto *shared =
              RenderGraph::sharedNode(ref.from.edge, ref.to.edge)) {
        nd = shared;
      }
    }

    struct NormalizedSlots {
      int from;
      int to;
    };

    auto normalizeSlots = [nd](const InnerGeom &geom) {
      NormalizedSlots slots{geom.slotFrom, geom.slotTo};

      if (geom.from.edge) {
        bool fromInv = geom.from.edge->getTo() == nd;
        slots.from = !fromInv
                         ? geom.slotFrom
                         : (geom.from.edge->pl().getLines().size() - 1 -
                            geom.slotFrom);
      }

      if (geom.to.edge) {
        bool toInv = geom.to.edge->getTo() == nd;
        slots.to = !toInv
                       ? geom.slotTo
                       : (geom.to.edge->pl().getLines().size() - 1 -
                          geom.slotTo);
      }

      return slots;
    };

    NormalizedSlots refSlots = normalizeSlots(ref);

    for (size_t i = 0; i < c.geoms.size(); i++) {
      PolyLine<double> pl = c.geoms[i].geom;
      NormalizedSlots curSlots = normalizeSlots(c.geoms[i]);

      if (ref.geom.getLength() >
          (_cfg->lineWidth + 2 * _cfg->outlineWidth + _cfg->lineSpacing) * 4) {
        double off =
            -(_cfg->lineWidth + _cfg->lineSpacing + 2 * _cfg->outlineWidth) *
            (static_cast<int>(curSlots.from) -
             static_cast<int>(refSlots.from));

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
  std::stringstream styleOutline;
  styleOutline << "fill:none;stroke:#000000;stroke-linecap:round;stroke-width:"
               << (width + 2 * _cfg->outlineWidth) * _cfg->outputResolution
               << ";"
               << oCss;
  Params paramsOutline;
  paramsOutline["style"] = styleOutline.str();
  paramsOutline["class"] = "transit-edge-outline " + getLineClass(line.id());

  std::stringstream styleStr;
  styleStr << "fill:none;stroke:#" << line.color() << ";" << css;

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
void SvgRenderer::renderArrowHead(const PolyLine<double> &p, double width,
                                  bool flipDir, bool atStart) {
  if (p.getLine().size() < 2)
    return;

  const DPoint *a;
  const DPoint *b;
  if (atStart) {
    a = &(*p.getLine().begin());
    b = &(*(p.getLine().begin() + 1));
  } else {
    a = &(*(p.getLine().end() - 2));
    b = &(*(p.getLine().end() - 1));
  }

  double dx = b->getX() - a->getX();
  double dy = b->getY() - a->getY();
  if (flipDir) {
    dx = -dx;
    dy = -dy;
  }
  double len = std::sqrt(dx * dx + dy * dy);
  if (len == 0)
    return;
  dx /= len;
  dy /= len;

  const DPoint &anchor = atStart ? *a : *b;

  ArrowHead ah;
  auto addPt = [&](double x, double y) {
    double rx = dx * x - dy * y;
    double ry = dy * x + dx * y;
    ah.pts.emplace_back(anchor.getX() + rx, anchor.getY() + ry);
  };

  addPt(0.0, -0.5 * width);
  addPt(0.0, 0.5 * width);
  addPt(0.5 * width, 0.5 * width);
  addPt(1.3 * width, 0.0);
  addPt(0.5 * width, -0.5 * width);

  _arrowHeads.push_back(ah);
}

// _____________________________________________________________________________
bool SvgRenderer::edgeHasSharpAngle(const PolyLine<double> &center,
                                    const shared::linegraph::LineEdge *e,
                                    const shared::linegraph::Line *line,
                                    bool markAdjacent) {
  const auto &pts = center.getLine();
  const double sharpTurnCos = std::cos(_cfg->sharpTurnAngle);

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
    if (cosang < sharpTurnCos) {
      return true;
    }
  }

  double checkDist = 10.0;
  auto checkNode = [&](const shared::linegraph::LineNode *n,
                       bool fromStart) -> bool {
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
    bool sharp = false;
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
      if (cosang < sharpTurnCos) {
        sharp = true;
        if (markAdjacent) {
          _forceDirMarker[line].insert(ne);
        }
      }
    }
    return sharp;
  };

  if (checkNode(e->getFrom(), true))
    return true;
  if (checkNode(e->getTo(), false))
    return true;

  return false;
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

  if (edgeHasSharpAngle(center, e, line, true)) {
    return true;
  }

  if (_edgesSinceMarker[line] >= static_cast<int>(_cfg->dirMarkerSpacing)) {
    return true;
  }

  return false;
}

// _____________________________________________________________________________
bool SvgRenderer::hasSharpAngle(const shared::linegraph::LineEdge *e,
                                const PolyLine<double> &center,
                                const shared::linegraph::Line *line) {
  return edgeHasSharpAngle(center, e, line, false);
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
    double tailWorld = 15.0 / _cfg->outputResolution;
    double minLengthForTail = arrowLength * 3 + tailWorld;
    bool sharpAngle = hasSharpAngle(e, center, line);
    double pLen = p.getLength();
    bool useTail = pLen > minLengthForTail &&
                   (_cfg->tailIgnoreSharpAngle || !sharpAngle);

    std::string css, oCss;

    if (!lo.style.isNull()) {
      css = lo.style.get().getCss();
      oCss = lo.style.get().getOutlineCss();
    }

    _edgesSinceMarker[line]++;

    bool needMarker =
        (_cfg->renderDirMarkers && needsDirMarker(e, center, line)) ||
        sharpAngle;
    bool drawMarker = needMarker && pLen > arrowLength * 3;

    if (drawMarker) {
      _edgesSinceMarker[line] = 0;
    } else if (needMarker) {
      // Edge is too short to draw a marker but one is needed; clamp the
      // counter so we do not exceed the forcing threshold.
      _edgesSinceMarker[line] = std::min(
          _edgesSinceMarker[line], static_cast<int>(_cfg->dirMarkerSpacing));
    }

    if (drawMarker) {
      if (lo.direction == 0 && !_cfg->renderBiDirMarker) {
        renderLinePart(p, lineW, *line, css, oCss);
      } else {
        PolyLine<double> firstPart = p.getSegmentAtDist(0, p.getLength() / 2);
        PolyLine<double> secondPart =
            p.getSegmentAtDist(p.getLength() / 2, p.getLength());
        PolyLine<double> revSecond = secondPart.reversed();

        if (lo.direction == 0) {
          double mid = p.getLength() / 2;
          double tailStart = mid - tailWorld / 2;
          double tailEnd = mid + tailWorld / 2;

          PolyLine<double> firstHalf = p.getSegmentAtDist(0, mid);
          PolyLine<double> secondHalf = p.getSegmentAtDist(mid, p.getLength());

          if (useTail) {
            PolyLine<double> tailToStart = p.getSegmentAtDist(tailStart, mid);
            PolyLine<double> tailToEnd = p.getSegmentAtDist(mid, tailEnd);
            renderLinePart(tailToStart, lineW, *line, "stroke:black",
                           "stroke:none");
            renderArrowHead(tailToStart, lineW, false, true);
            renderLinePart(tailToEnd, lineW, *line, "stroke:black",
                           "stroke:none");
            renderArrowHead(tailToEnd, lineW);
          }

          renderLinePart(firstHalf, lineW, *line, css, oCss);
          renderArrowHead(firstHalf, lineW, false, true);
          renderLinePart(secondHalf, lineW, *line, css, oCss);
          renderArrowHead(secondHalf, lineW);
        } else if (lo.direction == e->getTo()) {
          if (useTail) {
            double tailStart = std::max(0.0, firstPart.getLength() - tailWorld);
            PolyLine<double> tail =
                firstPart.getSegmentAtDist(tailStart, firstPart.getLength());

            renderLinePart(tail, lineW, *line, "stroke:black", "stroke:none");
            renderArrowHead(tail, lineW);
          }
          renderLinePart(firstPart, lineW, *line, css, oCss);
          renderArrowHead(firstPart, lineW);
          renderLinePart(revSecond, lineW, *line, css, oCss);
        } else {
          if (useTail) {
            double tailStart = std::max(0.0, revSecond.getLength() - tailWorld);
            PolyLine<double> tail =
                revSecond.getSegmentAtDist(tailStart, revSecond.getLength());
            renderLinePart(tail, lineW, *line, "stroke:black", "stroke:none");
            renderArrowHead(tail, lineW);
          }
          renderLinePart(revSecond, lineW, *line, css, oCss);
          renderArrowHead(revSecond, lineW);
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
  if (!params.count("class")) {
    params["class"] = "station-poly";
  }
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

  std::vector<std::map<std::string, std::string>> paths;
  std::vector<std::string> pathIds;

  size_t id = 0;
  const auto &labels = labeller.getStationLabels();
  bool wantHighlight =
      _cfg->highlightMeStationLabel && !_cfg->meStationId.empty();

  for (const auto &label : labels) {
    auto textPath = label.geom;
    double ang = util::geo::angBetween(textPath.front(), textPath.back());

    if ((fabs(ang) < (3 * M_PI / 2)) && (fabs(ang) > (M_PI / 2))) {
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

    std::string idStr = "stlblp" + util::toString(id++);

    pathPars["d"] = points.str();
    pathPars["id"] = idStr;

    paths.push_back(pathPars);
    pathIds.push_back(idStr);
  }

  _w.openTag("defs");
  for (auto &pathPars : paths) {
    _w.openTag("path", pathPars);
    _w.closeTag();
  }
  _w.closeTag();

  id = 0;
  for (const auto &label : labels) {
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

    std::map<std::string, std::string> params;
    params["class"] = "station-label";
    params["font-weight"] = label.bold ? "bold" : "normal";
    params["font-family"] = "TT Norms Pro";
    params["dy"] = shift;
    double fontSize = label.fontSize * _cfg->outputResolution;
    if (_cfg->fontSvgMax >= 0 && fontSize > _cfg->fontSvgMax)
      fontSize = _cfg->fontSvgMax;
    params["font-size"] = util::toString(fontSize);

    bool isMeLabel = false;
    if (wantHighlight && _meStationLabelVisual.isNull()) {
      std::string sanitized = util::sanitizeStationLabel(label.s.name);
      if (sanitized == _cfg->meStationId) {
        StationLabelVisual info;
        info.label = &label;
        info.pathId = pathIds[id];
        info.shift = shift;
        info.textAnchor = textAnchor;
        info.startOffset = startOffset;
        info.fontSizePx = fontSize;
        info.bold = label.bold;
        _meStationLabelVisual = info;
        isMeLabel = true;
      }
    }

    if (!isMeLabel) {
      _w.openTag("text", params);
      std::map<std::string, std::string> attrs;
      attrs["dy"] = shift;
      attrs["xlink:href"] = "#" + pathIds[id];
      attrs["startOffset"] = startOffset;
      attrs["text-anchor"] = textAnchor;
      _w.openTag("textPath", attrs);

      _w.writeText(label.s.name);
      _w.closeTag();
      _w.closeTag();
    }
    id++;
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::renderLineLabels(const Labeller &labeller,
                                   const RenderParams &rparams) {
  _w.openTag("g");

  std::vector<std::map<std::string, std::string>> paths;
  std::vector<std::string> pathIds;

  size_t id = 0;
  const auto &labels = labeller.getLineLabels();

  for (const auto &label : labels) {
    auto textPath = label.geom;
    double ang = util::geo::angBetween(textPath.front(), textPath.back());

    if ((fabs(ang) < (3 * M_PI / 2)) && (fabs(ang) > (M_PI / 2))) {
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

    std::string idStr = "textp" + util::toString(id++);

    pathPars["d"] = points.str();
    pathPars["id"] = idStr;

    paths.push_back(pathPars);
    pathIds.push_back(idStr);
  }

  _w.openTag("defs");
  for (auto &pathPars : paths) {
    _w.openTag("path", pathPars);
    _w.closeTag();
  }
  _w.closeTag();

  id = 0;
  for (const auto &label : labels) {
    auto textPath = label.geom;
    double ang = util::geo::angBetween(textPath.front(), textPath.back());
    double shift = 0.0;

    if ((fabs(ang) < (3 * M_PI / 2)) && (fabs(ang) > (M_PI / 2))) {
      shift = 0.75;
      textPath.reverse();
    }

    auto emitRow = [&](size_t start, size_t end, double rowOff) {
      std::map<std::string, std::string> params;
      params["class"] = "line-label";
      params["font-weight"] = "bold";
      params["font-family"] = "TT Norms Pro";
      params["dy"] = util::toString(shift + rowOff) + "em";
      params["font-size"] =
          util::toString(label.fontSize * _cfg->outputResolution);

      _w.openTag("text", params);
      std::map<std::string, std::string> attrs;
      attrs["dy"] = util::toString(shift + rowOff) + "em";
      attrs["xlink:href"] = "#" + pathIds[id];
      attrs["text-anchor"] = "middle";
      attrs["startOffset"] = "50%";
      _w.openTag("textPath", attrs);

      double dx = 0;
      for (size_t i = start; i < end; ++i) {
        auto line = label.lines[i];
        std::map<std::string, std::string> attrs;
        attrs["fill"] = "#" + line->color();
        attrs["dx"] = util::toString(dx);
        _w.openTag("tspan", attrs);
        dx = (label.fontSize * _cfg->outputResolution) / 3;
        _w.writeText(line->label());
        _w.closeTag();
      }
      _w.closeTag();
      _w.closeTag();
    };

    if (_cfg->compactRouteLabel && label.lines.size() > 4) {
      size_t rows = (label.lines.size() + 3) / 4;
      size_t perRow = (label.lines.size() + rows - 1) / rows;
      double offsetStep = 1.2;
      for (size_t r = 0; r < rows; ++r) {
        size_t start = r * perRow;
        size_t end = std::min(start + perRow, label.lines.size());
        double rowOff =
            (static_cast<double>(r) - (rows - 1) / 2.0) * offsetStep;
        emitRow(start, end, rowOff);
      }
    } else {
      emitRow(0, label.lines.size(), 0.0);
    }

    id++;
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

  struct AABB {
    double minX = 0.0;
    double minY = 0.0;
    double maxX = 0.0;
    double maxY = 0.0;

    bool intersects(const AABB &other) const {
      return !(maxX < other.minX || other.maxX < minX || maxY < other.minY ||
               other.maxY < minY);
    }
  };

  std::vector<AABB> placedStackBoxes;

  for (auto n : g.getNds()) {
    std::vector<const Line *> lines;
    std::unordered_set<const Line *> seen;
    for (auto e : n->getAdjList()) {
      for (const auto &lo : e->pl().getLines()) {
        if (seen.insert(lo.line).second && g.lineTerminatesAt(n, lo.line)) {
          lines.push_back(lo.line);
        }
      }
    }
    if (lines.empty())
      continue;

    std::sort(lines.begin(), lines.end(), [](const Line *lhs, const Line *rhs) {
      if (lhs->label() != rhs->label())
        return lhs->label() < rhs->label();
      if (lhs->color() != rhs->color())
        return lhs->color() < rhs->color();
      return lhs->id() < rhs->id();
    });

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

    double footprintMinX = std::numeric_limits<double>::max();
    double footprintMaxX = std::numeric_limits<double>::lowest();
    double footprintMinY = std::numeric_limits<double>::max();
    double footprintMaxY = std::numeric_limits<double>::lowest();
    bool hasFootprint = false;
    auto stopGeoms = g.getStopGeoms(n, _cfg->tightStations, 32);
    for (const auto &poly : stopGeoms) {
      const auto &outer = poly.getOuter();
      if (outer.empty())
        continue;
      for (const auto &p : outer) {
        footprintMinX = std::min(footprintMinX, p.getX());
        footprintMaxX = std::max(footprintMaxX, p.getX());
        footprintMinY = std::min(footprintMinY, p.getY());
        footprintMaxY = std::max(footprintMaxY, p.getY());
        hasFootprint = true;
      }
    }

    double footprintCenterX = nodeX;
    double footprintCenterY = nodeY;
    double footprintHalfHeight = 0.0;
    if (hasFootprint) {
      footprintCenterX = (footprintMinX + footprintMaxX) / 2.0;
      footprintCenterY = (footprintMinY + footprintMaxY) / 2.0;
      footprintHalfHeight = (footprintMaxY - footprintMinY) / 2.0;
    }

    double labelCenterX = 0.0;
    double labelCenterY = 0.0;
    double labelVExtent = 0.0;
    double labelMinX = std::numeric_limits<double>::max();
    double labelMaxX = std::numeric_limits<double>::lowest();
    double labelMinY = std::numeric_limits<double>::max();
    double labelMaxY = std::numeric_limits<double>::lowest();
    bool hasLabelGeom = false;
    if (sLbl) {
      const auto &base = sLbl->band[0];
      const auto &top = sLbl->band[2];
      double baseX = base[0].getX();
      double baseY = base[0].getY();
      double scale = sLbl->fontSize / _cfg->stationLabelSize;

      for (const auto &ln : sLbl->band) {
        for (const auto &p : ln) {
          double x = baseX + (p.getX() - baseX) * scale;
          double y = baseY + (p.getY() - baseY) * scale;
          labelMinX = std::min(labelMinX, x);
          labelMaxX = std::max(labelMaxX, x);
          labelMinY = std::min(labelMinY, y);
          labelMaxY = std::max(labelMaxY, y);
        }
      }

      if (labelMinX <= labelMaxX && labelMinY <= labelMaxY) {
        labelCenterX = (labelMinX + labelMaxX) / 2.0;
        labelCenterY = (labelMinY + labelMaxY) / 2.0;

        double dx = (base[1].getX() - base[0].getX()) * scale;
        double dy = (base[1].getY() - base[0].getY()) * scale;
        double width = std::sqrt(dx * dx + dy * dy);
        double angle = std::atan2(dy, dx);
        double hdx = (top[0].getX() - base[0].getX()) * scale;
        double hdy = (top[0].getY() - base[0].getY()) * scale;
        double height = std::sqrt(hdx * hdx + hdy * hdy);

        double vExtent = 0.0;
        double absTan = std::abs(std::tan(angle));
        double cosAngle = std::cos(angle);
        double sinAngle = std::sin(angle);
        if (height > 0 && std::abs(cosAngle) > 1e-9 &&
            (width == 0.0 || absTan <= width / height)) {
          vExtent = std::abs(height / (2.0 * cosAngle));
        } else if (std::abs(sinAngle) > 1e-9) {
          vExtent = std::abs(width / (2.0 * sinAngle));
        } else {
          vExtent = height / 2.0;
        }

        labelVExtent = vExtent;
        hasLabelGeom = true;
      }
    }

    bool above = true;
    if (hasLabelGeom) {
      above = labelCenterY > nodeY;
    } else if (hasFootprint) {
      above = footprintCenterY > nodeY;
    }

    double anchorX = nodeX;
    double anchorY = nodeY;
    double clearance = 0.0;
    bool anchorIsNode = false;

    switch (_cfg->terminusLabelAnchor) {
    case TerminusLabelAnchor::StationLabel:
      if (hasLabelGeom) {
        anchorX = labelCenterX;
        anchorY = above ? labelCenterY + labelVExtent
                        : labelCenterY - labelVExtent;
        clearance = labelVExtent;
      } else if (hasFootprint) {
        anchorX = footprintCenterX;
        anchorY = above ? footprintCenterY + footprintHalfHeight
                        : footprintCenterY - footprintHalfHeight;
        clearance = footprintHalfHeight;
      }
      break;
    case TerminusLabelAnchor::StopFootprint:
      if (hasFootprint) {
        anchorX = footprintCenterX;
        anchorY = above ? footprintCenterY + footprintHalfHeight
                        : footprintCenterY - footprintHalfHeight;
        clearance = footprintHalfHeight;
      }
      break;
    case TerminusLabelAnchor::Node:
      anchorX = nodeX;
      anchorY = nodeY;
      clearance = 0.0;
      anchorIsNode = true;
      break;
    }

    double x = (anchorX - rparams.xOff) * _cfg->outputResolution;
    double y =
        rparams.height - (anchorY - rparams.yOff) * _cfg->outputResolution;

    double fontSize = _cfg->lineLabelSize * _cfg->outputResolution;
    double padTop = fontSize * 0.28;
    double padBottom = fontSize * 0.12;
    double padX = fontSize * 0.2;
    double boxH = fontSize + padTop + padBottom;
    double charW = fontSize * 0.6;
    double boxR = padX * 2;

    size_t idx = 0;
    // Use a uniform gap to achieve consistent spacing regardless of the
    // orientation of the station label. The gap is configurable to allow
    // tuning without recompilation.
    double boxGap = _cfg->routeLabelBoxGap * _cfg->outputResolution;
    double terminusGap = _cfg->routeLabelTerminusGap * _cfg->outputResolution;
    double step = boxH + boxGap;

    // Use a constant label width based on five characters plus padding
    // so that all route label boxes share uniform dimensions.
    double uniformBoxW = 5 * charW + padX * 2;
    size_t linesPerCol = _cfg->compactTerminusLabel ? 4 : lines.size();
    if (linesPerCol == 0)
      linesPerCol = 1;
    size_t numCols = (lines.size() + linesPerCol - 1) / linesPerCol;
    double totalW = numCols * uniformBoxW + (numCols - 1) * boxGap;
    double baseStartX = x - totalW / 2;

    std::vector<double> colHeights(numCols, 0.0);
    double maxColumnHeight = 0.0;
    for (size_t col = 0; col < numCols; ++col) {
      size_t startIdx = col * linesPerCol;
      if (startIdx >= lines.size())
        break;
      size_t count = std::min(linesPerCol, lines.size() - startIdx);
      double colHeight =
          count > 0 ? count * boxH + (count - 1) * boxGap : 0.0;
      colHeights[col] = colHeight;
      maxColumnHeight = std::max(maxColumnHeight, colHeight);
    }

    double stationHalfHeight =
        anchorIsNode ? 0.0 : std::abs(clearance * _cfg->outputResolution);

    double stackCenterOffset = stationHalfHeight + terminusGap;
    auto computeCenterY = [&](bool placeAbove) {
      return placeAbove ? y - stackCenterOffset : y + stackCenterOffset;
    };

    auto buildStackBox = [&](double startPx, double centerPx) {
      AABB box;
      box.minX = std::min(startPx, startPx + totalW);
      box.maxX = std::max(startPx, startPx + totalW);
      if (maxColumnHeight > 0) {
        box.minY = centerPx - maxColumnHeight / 2.0;
        box.maxY = centerPx + maxColumnHeight / 2.0;
      } else {
        box.minY = centerPx;
        box.maxY = centerPx;
      }
      return box;
    };

    auto toPixelX = [&](double mapX) {
      return (mapX - rparams.xOff) * _cfg->outputResolution;
    };
    auto toPixelY = [&](double mapY) {
      return rparams.height - (mapY - rparams.yOff) * _cfg->outputResolution;
    };

    AABB footprintBox;
    bool hasFootprintBox = false;
    if (hasFootprint) {
      double fx1 = toPixelX(footprintMinX);
      double fx2 = toPixelX(footprintMaxX);
      double fy1 = toPixelY(footprintMinY);
      double fy2 = toPixelY(footprintMaxY);
      footprintBox.minX = std::min(fx1, fx2);
      footprintBox.maxX = std::max(fx1, fx2);
      footprintBox.minY = std::min(fy1, fy2);
      footprintBox.maxY = std::max(fy1, fy2);
      hasFootprintBox = true;
    }

    AABB labelBox;
    bool hasLabelBox = false;
    if (hasLabelGeom && labelMinX <= labelMaxX && labelMinY <= labelMaxY) {
      double lx1 = toPixelX(labelMinX);
      double lx2 = toPixelX(labelMaxX);
      double ly1 = toPixelY(labelMinY);
      double ly2 = toPixelY(labelMaxY);
      labelBox.minX = std::min(lx1, lx2);
      labelBox.maxX = std::max(lx1, lx2);
      labelBox.minY = std::min(ly1, ly2);
      labelBox.maxY = std::max(ly1, ly2);
      hasLabelBox = true;
    }

    int horizontalPreference = 0;
    if (hasLabelGeom) {
      constexpr double tol = 1e-6;
      if (nodeX < labelMinX - tol) {
        horizontalPreference = 1;
      } else if (nodeX > labelMaxX + tol) {
        horizontalPreference = -1;
      } else if (labelCenterX > nodeX + tol) {
        horizontalPreference = 1;
      } else if (labelCenterX < nodeX - tol) {
        horizontalPreference = -1;
      }
    }

    auto collidesWith = [&](const AABB &candidate) {
      if (hasFootprintBox && candidate.intersects(footprintBox))
        return true;
      if (hasLabelBox && candidate.intersects(labelBox))
        return true;
      for (const auto &placed : placedStackBoxes) {
        if (candidate.intersects(placed))
          return true;
      }
      return false;
    };

    double shiftDistance = uniformBoxW + boxGap;
    int maxShiftSteps = std::max(0, _cfg->terminusLabelMaxLateralShift);
    bool selectedAbove = above;
    double selectedStartX = baseStartX;
    double selectedCenterY = computeCenterY(selectedAbove);
    AABB selectedBox = buildStackBox(selectedStartX, selectedCenterY);
    bool foundPlacement = false;

    std::array<bool, 2> orientationOrder = {above, !above};
    for (int step = 0; step <= maxShiftSteps && !foundPlacement; ++step) {
      std::vector<double> shifts;
      if (step == 0) {
        shifts.push_back(0.0);
      } else {
        double delta = shiftDistance * step;
        if (horizontalPreference > 0) {
          shifts.push_back(delta);
          shifts.push_back(-delta);
        } else {
          shifts.push_back(-delta);
          shifts.push_back(delta);
        }
      }
      for (bool orientationOption : orientationOrder) {
        for (double shift : shifts) {
          double candidateStartX = baseStartX + shift;
          double candidateCenterY = computeCenterY(orientationOption);
          AABB candidateBox = buildStackBox(candidateStartX, candidateCenterY);
          if (collidesWith(candidateBox)) {
            continue;
          }
          selectedAbove = orientationOption;
          selectedStartX = candidateStartX;
          selectedCenterY = candidateCenterY;
          selectedBox = candidateBox;
          foundPlacement = true;
          break;
        }
        if (foundPlacement)
          break;
      }
    }

    if (!foundPlacement) {
      selectedAbove = above;
      selectedStartX = baseStartX;
      selectedCenterY = computeCenterY(selectedAbove);
      selectedBox = buildStackBox(selectedStartX, selectedCenterY);
    }

    double startX = selectedStartX;
    double stackCenterY = selectedCenterY;
    above = selectedAbove;

    std::vector<double> colTops(numCols, stackCenterY);
    for (size_t col = 0; col < numCols; ++col) {
      double colHeight = colHeights[col];
      if (colHeight > 0) {
        colTops[col] = stackCenterY - colHeight / 2.0;
      } else if (maxColumnHeight > 0) {
        colTops[col] = stackCenterY - maxColumnHeight / 2.0;
      }
    }

    placedStackBoxes.push_back(selectedBox);

    if (_cfg->compactTerminusLabel) {
      for (auto line : lines) {
        std::string label = line->label();
        double boxW = uniformBoxW;
        size_t col = idx / linesPerCol;
        size_t row = idx % linesPerCol;
        double rectX = startX + col * (uniformBoxW + boxGap);
        double rectY = colTops[col] + row * step;

        std::string fillColor = line->color();
        std::string textColor = isLightColor(fillColor) ? "black" : "white";

        {
          std::map<std::string, std::string> attrs;
          attrs["x"] = util::toString(rectX);
          attrs["y"] = util::toString(rectY);
          attrs["width"] = util::toString(boxW);
          attrs["height"] = util::toString(boxH);
          attrs["rx"] = util::toString(boxR);
          attrs["ry"] = util::toString(boxR);
          attrs["fill"] = "#" + fillColor;
          _w.openTag("rect", attrs);
        }
        _w.closeTag();

        {
          std::map<std::string, std::string> attrs;
          attrs["class"] = "line-label";
          attrs["font-weight"] = "bold";
          attrs["font-family"] = "TT Norms Pro";
          attrs["text-anchor"] = "middle";
          attrs["dominant-baseline"] = "middle";
          attrs["alignment-baseline"] = "middle";
          attrs["font-size"] = util::toString(fontSize);
          attrs["fill"] = textColor;
          attrs["x"] = util::toString(rectX + boxW / 2);
          attrs["y"] = util::toString(rectY + padTop + fontSize / 2);
          _w.openTag("text", attrs);
        }

        _w.writeText(label);
        _w.closeTag();
        idx++;
      }
    } else {
      for (auto line : lines) {
        std::string label = line->label();
        double boxW = uniformBoxW;
        size_t col = idx / linesPerCol;
        size_t row = idx % linesPerCol;
        double rectX = startX + col * (uniformBoxW + boxGap);
        double rectY = colTops[col] + row * step;

        std::string fillColor = line->color();
        std::string textColor = isLightColor(fillColor) ? "black" : "white";

        {
          std::map<std::string, std::string> attrs;
          attrs["x"] = util::toString(rectX);
          attrs["y"] = util::toString(rectY);
          attrs["width"] = util::toString(boxW);
          attrs["height"] = util::toString(boxH);
          attrs["rx"] = util::toString(boxR);
          attrs["ry"] = util::toString(boxR);
          attrs["fill"] = "#" + fillColor;
          _w.openTag("rect", attrs);
        }
        _w.closeTag();

        {
          std::map<std::string, std::string> attrs;
          attrs["class"] = "line-label";
          attrs["font-weight"] = "bold";
          attrs["font-family"] = "TT Norms Pro";
          attrs["text-anchor"] = "middle";
          attrs["dominant-baseline"] = "middle";
          attrs["alignment-baseline"] = "middle";
          attrs["font-size"] = util::toString(fontSize);
          attrs["fill"] = textColor;
          attrs["x"] = util::toString(rectX + boxW / 2);
          attrs["y"] = util::toString(rectY + padTop + fontSize / 2);
          _w.openTag("text", attrs);
        }

        _w.writeText(label);
        _w.closeTag();
        idx++;
      }
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
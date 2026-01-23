// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdint.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ostream>
#include <vector>

#include "shared/linegraph/Line.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/label/Labeller.h"
#include "transitmap/output/ImageCodec.h"
#include "transitmap/output/MBTilesReader.h"
#include "transitmap/output/SvgRenderer.h"
#include "util/Base64.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using shared::linegraph::Line;
using shared::linegraph::LineNode;
using shared::rendergraph::InnerGeom;
using shared::rendergraph::RenderGraph;
using transitmapper::config::Config;
using transitmapper::label::Labeller;
using transitmapper::output::InnerClique;
using transitmapper::output::MBTilesMetadata;
using transitmapper::output::MBTilesReader;
using transitmapper::output::RgbaImage;
using transitmapper::output::RenderParams;
using transitmapper::output::SvgRenderer;
using util::geo::DPoint;
using util::geo::DPolygon;
using util::geo::LinePoint;
using util::geo::LinePointCmp;
using util::geo::Polygon;
using util::geo::PolyLine;
using util::DEBUG;

namespace {
struct PaperSpec {
  int widthMm = 0;
  int heightMm = 0;
  PaperSpec() {}
  PaperSpec(int w, int h) : widthMm(w), heightMm(h) {}
};

PaperSpec getPaperSpec(const std::string& paper) {
  if (paper == "A4") return PaperSpec(210, 297);
  if (paper == "A4L") return PaperSpec(297, 210);
  if (paper == "A3") return PaperSpec(297, 420);
  if (paper == "A3L") return PaperSpec(420, 297);
  return PaperSpec(297, 210);
}

int computeCanvasHeight(int canvasWidth, const std::string& paper) {
  PaperSpec spec = getPaperSpec(paper);
  if (spec.widthMm == 0) return canvasWidth;
  double ratio =
      static_cast<double>(spec.heightMm) / static_cast<double>(spec.widthMm);
  return static_cast<int>(std::round(canvasWidth * ratio));
}

util::geo::DBox padBox(const util::geo::DBox& box, double padX, double padY) {
  util::geo::DPoint ll(box.getLowerLeft().getX() - padX,
                       box.getLowerLeft().getY() - padY);
  util::geo::DPoint ur(box.getUpperRight().getX() + padX,
                       box.getUpperRight().getY() + padY);
  return util::geo::DBox(ll, ur);
}

util::geo::DBox padBoxSides(const util::geo::DBox& box, double padLeft,
                            double padRight, double padBottom,
                            double padTop) {
  util::geo::DPoint ll(box.getLowerLeft().getX() - padLeft,
                       box.getLowerLeft().getY() - padBottom);
  util::geo::DPoint ur(box.getUpperRight().getX() + padRight,
                       box.getUpperRight().getY() + padTop);
  return util::geo::DBox(ll, ur);
}

util::geo::DBox expandBoxToAspect(const util::geo::DBox& box, double aspect) {
  if (aspect <= 0) return box;
  double width = box.getUpperRight().getX() - box.getLowerLeft().getX();
  double height = box.getUpperRight().getY() - box.getLowerLeft().getY();
  if (width <= 0 || height <= 0) return box;

  double current = height / width;
  if (fabs(current - aspect) < 1e-9) return box;

  util::geo::DPoint ll = box.getLowerLeft();
  util::geo::DPoint ur = box.getUpperRight();
  if (current < aspect) {
    double newH = width * aspect;
    double delta = (newH - height) / 2.0;
    ll.setY(ll.getY() - delta);
    ur.setY(ur.getY() + delta);
  } else {
    double newW = height / aspect;
    double delta = (newW - width) / 2.0;
    ll.setX(ll.getX() - delta);
    ur.setX(ur.getX() + delta);
  }

  return util::geo::DBox(ll, ur);
}

double flippedY(const RenderParams& rparams, double y) {
  return (2.0 * rparams.yOff + rparams.height) - y;
}

const double kEarthRadius = 6378137.0;
const double kOriginShift = 2.0 * M_PI * kEarthRadius / 2.0;
const double kInitialResolution = 2.0 * M_PI * kEarthRadius / 256.0;
const int kTileSize = 256;

double resolutionAtZoom(int z) {
  return kInitialResolution / static_cast<double>(1 << z);
}

struct TileRange {
  int minX = 0;
  int maxX = 0;
  int minY = 0;
  int maxY = 0;
  size_t tileCount() const {
    if (maxX < minX || maxY < minY) return 0;
    return static_cast<size_t>(maxX - minX + 1) *
           static_cast<size_t>(maxY - minY + 1);
  }
};

TileRange tileRangeForBBox(const util::geo::DBox& bbox, int zoom) {
  TileRange range;
  if (zoom < 0) return range;
  double res = resolutionAtZoom(zoom);

  double minPx = (bbox.getLowerLeft().getX() + kOriginShift) / res;
  double maxPx = (bbox.getUpperRight().getX() + kOriginShift) / res;
  double minPy = (kOriginShift - bbox.getUpperRight().getY()) / res;
  double maxPy = (kOriginShift - bbox.getLowerLeft().getY()) / res;

  int maxIndex = (1 << zoom) - 1;

  range.minX = std::max(0, static_cast<int>(std::floor(minPx / kTileSize)));
  range.maxX = std::min(maxIndex,
                        static_cast<int>(std::floor((maxPx - 1) / kTileSize)));
  range.minY = std::max(0, static_cast<int>(std::floor(minPy / kTileSize)));
  range.maxY = std::min(maxIndex,
                        static_cast<int>(std::floor((maxPy - 1) / kTileSize)));

  return range;
}

int nextLowerZoom(int current, const std::vector<int>& allowed, int minZoom) {
  if (allowed.empty()) return std::max(minZoom, current - 1);
  for (int i = static_cast<int>(allowed.size()) - 1; i >= 0; --i) {
    if (allowed[i] < current) return allowed[i];
  }
  return current;
}

int pickZoomAuto(const util::geo::DBox& bbox, int canvasW, int canvasH,
                 double oversample, int minZoom, int maxZoom, int maxTiles,
                 const std::vector<int>& zoomOverrides, size_t* tileCount) {
  double width = bbox.getUpperRight().getX() - bbox.getLowerLeft().getX();
  double height = bbox.getUpperRight().getY() - bbox.getLowerLeft().getY();
  if (canvasW <= 0 || canvasH <= 0 || width <= 0 || height <= 0) {
    if (tileCount) *tileCount = 0;
    return minZoom;
  }

  double resX = width / static_cast<double>(canvasW);
  double resY = height / static_cast<double>(canvasH);
  double res = std::max(resX, resY) * oversample;

  double zIdeal = std::log(kInitialResolution / res) / std::log(2.0);
  int z = static_cast<int>(std::floor(zIdeal));
  if (z < minZoom) z = minZoom;
  if (z > maxZoom) z = maxZoom;

  std::vector<int> allowed = zoomOverrides;
  if (!allowed.empty()) {
    std::vector<int> filtered;
    for (int level : allowed) {
      if (level >= minZoom && level <= maxZoom) filtered.push_back(level);
    }
    allowed.swap(filtered);
  }

  if (!allowed.empty()) {
    int chosen = allowed.front();
    for (int level : allowed) {
      if (level <= z) chosen = level;
      if (level > z) break;
    }
    z = chosen;
  }

  TileRange range = tileRangeForBBox(bbox, z);
  size_t count = range.tileCount();
  while (count > static_cast<size_t>(maxTiles) && z > minZoom) {
    int next = nextLowerZoom(z, allowed, minZoom);
    if (next == z) break;
    z = next;
    range = tileRangeForBBox(bbox, z);
    count = range.tileCount();
  }

  if (tileCount) *tileCount = count;
  return z;
}

struct BackgroundResult {
  std::string dataUri;
  util::geo::DBox bbox;
  int zoom = 0;
  size_t tileCount = 0;
};

bool buildMbtilesBackground(const Config* cfg,
                            const util::geo::DBox& bbox,
                            BackgroundResult* out) {
  if (!cfg || cfg->mbtilesPath.empty() || !out) return false;

  MBTilesReader reader(cfg->mbtilesPath);
  if (!reader.isOpen()) return false;

  const MBTilesMetadata& meta = reader.metadata();
  int canvasH = computeCanvasHeight(cfg->canvasWidth, cfg->paper);

  size_t tileCount = 0;
  int zoom = pickZoomAuto(bbox, cfg->canvasWidth, canvasH, cfg->oversample,
                          meta.minZoom, meta.maxZoom, cfg->maxTiles,
                          cfg->zoomLevels, &tileCount);

  TileRange range = tileRangeForBBox(bbox, zoom);
  if (range.tileCount() == 0) return false;

  const int tilesWide = range.maxX - range.minX + 1;
  const int tilesHigh = range.maxY - range.minY + 1;
  const int mosaicW = tilesWide * kTileSize;
  const int mosaicH = tilesHigh * kTileSize;

  RgbaImage mosaic;
  mosaic.width = mosaicW;
  mosaic.height = mosaicH;
  mosaic.rgba.assign(static_cast<size_t>(mosaicW * mosaicH * 4), 0);

  std::vector<unsigned char> tileBlob;
  for (int ty = range.minY; ty <= range.maxY; ++ty) {
    for (int tx = range.minX; tx <= range.maxX; ++tx) {
      tileBlob.clear();
      if (!reader.getTile(zoom, tx, ty, &tileBlob)) continue;
      RgbaImage tileImage;
      if (!decodeImageToRgba(meta.format, tileBlob.data(), tileBlob.size(),
                             &tileImage)) {
        continue;
      }
      if (tileImage.width <= 0 || tileImage.height <= 0) continue;

      int dstX = (tx - range.minX) * kTileSize;
      int dstY = (ty - range.minY) * kTileSize;
      int copyW = std::min(kTileSize, tileImage.width);
      int copyH = std::min(kTileSize, tileImage.height);

      for (int y = 0; y < copyH; ++y) {
        const unsigned char* srcRow =
            tileImage.rgba.data() + y * tileImage.width * 4;
        unsigned char* dstRow =
            mosaic.rgba.data() +
            (dstY + y) * mosaic.width * 4 + dstX * 4;
        std::memcpy(dstRow, srcRow, copyW * 4);
      }
    }
  }

  double res = resolutionAtZoom(zoom);
  double minPx = (bbox.getLowerLeft().getX() + kOriginShift) / res;
  double maxPx = (bbox.getUpperRight().getX() + kOriginShift) / res;
  double minPy = (kOriginShift - bbox.getUpperRight().getY()) / res;
  double maxPy = (kOriginShift - bbox.getLowerLeft().getY()) / res;

  int minPxI = static_cast<int>(std::floor(minPx));
  int maxPxI = static_cast<int>(std::ceil(maxPx));
  int minPyI = static_cast<int>(std::floor(minPy));
  int maxPyI = static_cast<int>(std::ceil(maxPy));

  int cropX = minPxI - range.minX * kTileSize;
  int cropY = minPyI - range.minY * kTileSize;
  int cropW = maxPxI - minPxI;
  int cropH = maxPyI - minPyI;

  cropX = std::max(0, cropX);
  cropY = std::max(0, cropY);
  cropW = std::min(cropW, mosaic.width - cropX);
  cropH = std::min(cropH, mosaic.height - cropY);

  if (cropW <= 0 || cropH <= 0) return false;

  RgbaImage cropped;
  cropped.width = cropW;
  cropped.height = cropH;
  cropped.rgba.resize(static_cast<size_t>(cropW * cropH * 4));
  for (int y = 0; y < cropH; ++y) {
    const unsigned char* srcRow = mosaic.rgba.data() +
                                  (cropY + y) * mosaic.width * 4 +
                                  cropX * 4;
    unsigned char* dstRow =
        cropped.rgba.data() + y * cropped.width * 4;
    std::memcpy(dstRow, srcRow, cropW * 4);
  }

  // Flip vertically so the image is upright after the SVG Y-flip transform.
  if (cropped.width > 0 && cropped.height > 0) {
    const size_t rowBytes = static_cast<size_t>(cropped.width) * 4;
    for (int y = 0; y < cropped.height / 2; ++y) {
      unsigned char* top = cropped.rgba.data() + y * rowBytes;
      unsigned char* bottom =
          cropped.rgba.data() + (cropped.height - 1 - y) * rowBytes;
      for (size_t i = 0; i < rowBytes; ++i) {
        std::swap(top[i], bottom[i]);
      }
    }
  }

  std::vector<unsigned char> pngBytes;
  if (!encodeRgbaToPng(cropped, &pngBytes)) return false;

  out->dataUri = "data:image/png;base64," + util::base64Encode(pngBytes);
  out->bbox = bbox;
  out->zoom = zoom;
  out->tileCount = tileCount;
  return true;
}
}  // namespace

// _____________________________________________________________________________
SvgRenderer::SvgRenderer(std::ostream* o, const config::Config* cfg)
    : _o(o), _w(o, true), _cfg(cfg) {}

// _____________________________________________________________________________
void SvgRenderer::print(const RenderGraph& outG) {
  std::map<std::string, std::string> params;
  RenderParams rparams;

  auto box = outG.getBBox();

  box = util::geo::pad(
      box, outG.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing));

  Labeller labeller(_cfg);
  if (_cfg->renderLabels) {
    LOGTO(DEBUG, std::cerr) << "Rendering labels...";
    labeller.label(outG, _cfg->dontLabelDeg2);
    box = util::geo::extendBox(labeller.getBBox(), box);
  }

  box = padBoxSides(box, _cfg->outputPaddingLeft, _cfg->outputPaddingRight,
                    _cfg->outputPaddingBottom, _cfg->outputPaddingTop);

  int canvasHeight = computeCanvasHeight(_cfg->canvasWidth, _cfg->paper);
  double targetAspect =
      static_cast<double>(canvasHeight) / static_cast<double>(_cfg->canvasWidth);

  if (!_cfg->mbtilesPath.empty()) {
    double width = box.getUpperRight().getX() - box.getLowerLeft().getX();
    double height = box.getUpperRight().getY() - box.getLowerLeft().getY();
    double padX = width * _cfg->backgroundPadPct;
    double padY = height * _cfg->backgroundPadPct;
    box = padBox(box, padX, padY);
  }

  // Ensure viewBox matches canvas aspect ratio to avoid stretching.
  box = expandBoxToAspect(box, targetAspect);

  rparams.xOff = box.getLowerLeft().getX();
  rparams.yOff = box.getLowerLeft().getY();

  rparams.width = box.getUpperRight().getX() - rparams.xOff;
  rparams.height = box.getUpperRight().getY() - rparams.yOff;

  if (!_cfg->worldFilePath.empty() && _cfg->canvasUnit == "px") {
    std::ofstream file;
    file.open(_cfg->worldFilePath);
    if (file) {
      double res = rparams.width / static_cast<double>(_cfg->canvasWidth);
      file << res << std::endl
           << 0 << std::endl
           << 0 << std::endl
           << -res << std::endl
           << std::fixed << rparams.xOff << std::endl
           << (rparams.yOff + rparams.height) << std::endl;
      file.close();
    }
  }

  auto latLngLL = util::geo::webMercToLatLng<double>(box.getLowerLeft().getX(),
                                                     box.getLowerLeft().getY());
  auto latLngUR = util::geo::webMercToLatLng<double>(
      box.getUpperRight().getX(), box.getUpperRight().getY());

  params["latlng-box"] = std::to_string(latLngLL.getX()) + "," +
                         std::to_string(latLngLL.getY()) + "," +
                         std::to_string(latLngUR.getX()) + "," +
                         std::to_string(latLngUR.getY());

  params["width"] = std::to_string(_cfg->canvasWidth) + _cfg->canvasUnit;
  params["height"] = std::to_string(canvasHeight) + _cfg->canvasUnit;
  params["viewBox"] = std::to_string(rparams.xOff) + " " +
                      std::to_string(rparams.yOff) + " " +
                      std::to_string(rparams.width) + " " +
                      std::to_string(rparams.height);
  params["xmlns"] = "http://www.w3.org/2000/svg";
  params["xmlns:xlink"] = "http://www.w3.org/1999/xlink";

  BackgroundResult bg;
  bool hasBg = false;
  if (!_cfg->mbtilesPath.empty()) {
    hasBg = buildMbtilesBackground(_cfg, box, &bg);
    if (hasBg) {
      LOGTO(INFO, std::cerr) << "MBTiles zoom selected: z=" << bg.zoom
                             << " tiles=" << bg.tileCount;
    } else {
      LOGTO(WARN, std::cerr) << "Failed to build MBTiles background";
    }
  }

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
  for (auto const& m : _markers) {
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

  std::map<std::string, std::string> gParams;
  gParams["transform"] =
      "translate(0 " +
      std::to_string(2.0 * rparams.yOff + rparams.height) +
      ") scale(1 -1)";
  _w.openTag("g", gParams);

  if (hasBg) {
    std::map<std::string, std::string> imgParams;
    imgParams["x"] = std::to_string(bg.bbox.getLowerLeft().getX());
    imgParams["y"] = std::to_string(bg.bbox.getLowerLeft().getY());
    imgParams["width"] =
        std::to_string(bg.bbox.getUpperRight().getX() -
                       bg.bbox.getLowerLeft().getX());
    imgParams["height"] =
        std::to_string(bg.bbox.getUpperRight().getY() -
                       bg.bbox.getLowerLeft().getY());
    imgParams["preserveAspectRatio"] = "xMidYMid slice";
    imgParams["opacity"] = std::to_string(_cfg->backgroundOpacity);
    imgParams["xlink:href"] = bg.dataUri;
    _w.openTag("image", imgParams);
    _w.closeTag();
  }

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

  _w.closeTag();

  LOGTO(DEBUG, std::cerr) << "Writing labels...";
  if (_cfg->renderLabels) {
    renderLineLabels(labeller, rparams);
    renderStationLabels(labeller, rparams);
  }

  _w.closeTags();
}

// _____________________________________________________________________________
void SvgRenderer::outputNodes(const RenderGraph& outG,
                              const RenderParams& rparams) {
  _w.openTag("g");
  for (auto n : outG.getNds()) {
    if (_cfg->renderStations && n->pl().stops().size() > 0 &&
        n->pl().fronts().size() > 0) {
      bool mergedRep = n->pl().stopmergeIsMergedRep() &&
                       n->pl().stopmergeMemberCount() >= 2;
      double legacyRad = (_cfg->lineWidth * 0.5) + _cfg->outlineWidth;
      double baseRad =
          (_cfg->stationRadius > 0.0 ? _cfg->stationRadius : legacyRad);
      baseRad *= _cfg->stationRadiusMult;
      baseRad = std::max(1.0, baseRad);
      DPoint center = *n->pl().getGeom();
      std::map<std::string, std::string> ringParams;
      ringParams["stroke"] = "black";
      ringParams["fill"] = "none";
      ringParams["stroke-width"] = util::toString(_cfg->lineWidth / 2.5);

      std::map<std::string, std::string> innerParams;
      innerParams["stroke"] = "black";
      innerParams["fill"] = "white";
      innerParams["stroke-width"] = util::toString(_cfg->lineWidth / 3.0);

      _w.openTag("g");
      if (mergedRep) {
        switch (_cfg->mergedStopView) {
          case Config::MergedStopView::NONE:
            printCircle(center, baseRad, innerParams, rparams);
            break;
          case Config::MergedStopView::DOUBLE_RING:
            printCircle(center, baseRad * 1.8, ringParams, rparams);
            printCircle(center, baseRad * 1.5, ringParams, rparams);
            printCircle(center, baseRad * 1.2, innerParams, rparams);
            break;
          case Config::MergedStopView::MERGE_COUNT:
            printCircle(center, baseRad, innerParams, rparams);
            break;
          case Config::MergedStopView::SIMPLE_CROSS: {
            printCircle(center, baseRad, innerParams, rparams);
            double crossLen = baseRad * _cfg->mergedCrossScale;
            crossLen = std::max(baseRad * 0.8, std::min(crossLen, baseRad * 1.5));
            double crossWidth =
                std::max(1.0, _cfg->outlineWidth * _cfg->mergedCrossStrokeMult);
            std::stringstream crossStyle;
            crossStyle << "fill:none;stroke:black;stroke-linecap:round;"
                       << "stroke-width:" << crossWidth;
            PolyLine<double> diag1(
                DPoint(center.getX() - crossLen, center.getY() - crossLen),
                DPoint(center.getX() + crossLen, center.getY() + crossLen));
            PolyLine<double> diag2(
                DPoint(center.getX() - crossLen, center.getY() + crossLen),
                DPoint(center.getX() + crossLen, center.getY() - crossLen));
            printLine(diag1, crossStyle.str(), rparams);
            printLine(diag2, crossStyle.str(), rparams);
            break;
          }
        }

        if (_cfg->mergedStopView == Config::MergedStopView::DOUBLE_RING ||
            _cfg->mergedStopView == Config::MergedStopView::MERGE_COUNT) {
          std::string badge =
              "×" + util::toString(n->pl().stopmergeMemberCount());
          std::map<std::string, std::string> badgeParams;
          badgeParams["x"] = util::toString(center.getX() + baseRad * 1.5);
          badgeParams["y"] = util::toString(center.getY() - baseRad * 1.2);
          badgeParams["font-size"] =
              util::toString(_cfg->stationLabelSize * 0.35);
          badgeParams["font-weight"] = "bold";
          badgeParams["fill"] = "black";
          _w.openTag("text", badgeParams);
          _w.writeText(badge);
          _w.closeTag();
        }
      } else {
        printCircle(center, baseRad, innerParams, rparams);
      }
      _w.closeTag();
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::renderNodeFronts(const RenderGraph& outG,
                                   const RenderParams& rparams) {
  _w.openTag("g");
  for (auto n : outG.getNds()) {
    std::string color = n->pl().stops().size() > 0 ? "red" : "black";
    for (auto& f : n->pl().fronts()) {
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
void SvgRenderer::outputEdges(const RenderGraph& outG,
                              const RenderParams& rparams) {
  struct cmp {
    bool operator()(const LineNode* lhs, const LineNode* rhs) const {
      return lhs->getAdjList().size() > rhs->getAdjList().size() ||
             (lhs->getAdjList().size() == rhs->getAdjList().size() &&
              RenderGraph::getConnCardinality(lhs) >
                  RenderGraph::getConnCardinality(rhs)) ||
             (lhs->getAdjList().size() == rhs->getAdjList().size() &&
              lhs > rhs);
    }
  };

  struct cmpEdge {
    bool operator()(const shared::linegraph::LineEdge* lhs,
                    const shared::linegraph::LineEdge* rhs) const {
      return lhs->pl().getLines().size() < rhs->pl().getLines().size() ||
             (lhs->pl().getLines().size() == rhs->pl().getLines().size() &&
              lhs < rhs);
    }
  };

  std::set<const LineNode*, cmp> nodesOrdered;
  std::set<const shared::linegraph::LineEdge*, cmpEdge> edgesOrdered;
  for (auto nd : outG.getNds()) nodesOrdered.insert(nd);

  std::set<const shared::linegraph::LineEdge*> rendered;

  for (const auto n : nodesOrdered) {
    edgesOrdered.insert(n->getAdjList().begin(), n->getAdjList().end());

    for (const auto* e : edgesOrdered) {
      if (rendered.insert(e).second) renderEdgeTripGeom(outG, e, rparams);
    }
  }
}

// _____________________________________________________________________________
void SvgRenderer::renderNodeConnections(const RenderGraph& outG,
                                        const LineNode* n,
                                        const RenderParams& rparams) {
  UNUSED(rparams);
  auto geoms = outG.innerGeoms(n, _cfg->innerGeometryPrecision);

  for (auto& clique : getInnerCliques(n, geoms, 9999)) renderClique(clique, n);
}

// _____________________________________________________________________________
std::multiset<InnerClique> SvgRenderer::getInnerCliques(
    const shared::linegraph::LineNode* n, std::vector<InnerGeom> pool,
    size_t level) const {
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
size_t SvgRenderer::getNextPartner(const InnerClique& forClique,
                                   const std::vector<InnerGeom>& pool,
                                   size_t level) const {
  for (size_t i = 0; i < pool.size(); i++) {
    const auto& ic = pool[i];
    for (auto& ciq : forClique.geoms) {
      if (isNextTo(ic, ciq) || (level > 1 && hasSameOrigin(ic, ciq))) {
        return i;
      }
    }
  }

  return pool.size();
}

// _____________________________________________________________________________
bool SvgRenderer::isNextTo(const InnerGeom& a, const InnerGeom& b) const {
  double THRESHOLD = 0.5 * M_PI + 0.1;

  if (!a.from.edge) return false;
  if (!b.from.edge) return false;
  if (!a.to.edge) return false;
  if (!b.to.edge) return false;

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
bool SvgRenderer::hasSameOrigin(const InnerGeom& a, const InnerGeom& b) const {
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
void SvgRenderer::renderClique(const InnerClique& cc, const LineNode* n) {
  _innerDelegates.push_back(
      std::map<uintptr_t, std::vector<OutlinePrintPair>>());
  std::multiset<InnerClique> renderCliques = getInnerCliques(n, cc.geoms, 0);
  for (const auto& c : renderCliques) {
    // the longest geom will be the ref geom
    InnerGeom ref = c.geoms[0];
    for (size_t i = 1; i < c.geoms.size(); i++) {
      if (c.geoms[i].geom.getLength() > ref.geom.getLength()) ref = c.geoms[i];
    }

    for (size_t i = 0; i < c.geoms.size(); i++) {
      PolyLine<double> pl = c.geoms[i].geom;

      if (ref.geom.getLength() >
          (_cfg->lineWidth + 2 * _cfg->outlineWidth + _cfg->lineSpacing) * 4) {
        double off =
            -(_cfg->lineWidth + _cfg->lineSpacing + 2 * _cfg->outlineWidth) *
            (static_cast<int>(c.geoms[i].slotFrom) -
             static_cast<int>(ref.slotFrom));

        if (ref.from.edge->getTo() == n) off = -off;

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
                          << (_cfg->lineWidth + _cfg->outlineWidth);
      Params paramsOutlineCropped;
      paramsOutlineCropped["style"] = styleOutlineCropped.str();
      paramsOutlineCropped["class"] += " inner-geom-outline";
      paramsOutlineCropped["class"] +=
          " " + getLineClass(c.geoms[i].from.line->id());

      std::stringstream styleStr;
      styleStr << "fill:none;stroke:#" << c.geoms[i].from.line->color();

      styleStr << ";stroke-linecap:round;stroke-opacity:1;stroke-width:"
               << _cfg->lineWidth;
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
                                 const Line& line, const std::string& css,
                                 const std::string& oCss) {
  renderLinePart(p, width, line, css, oCss, "");
}

// _____________________________________________________________________________
void SvgRenderer::renderLinePart(const PolyLine<double> p, double width,
                                 const Line& line, const std::string& css,
                                 const std::string& oCss,
                                 const std::string& endMarker) {
  std::stringstream styleOutline;
  styleOutline << "fill:none;stroke:#000000;stroke-linecap:round;stroke-width:"
               << (width + _cfg->outlineWidth) << ";"
               << oCss;
  Params paramsOutline;
  paramsOutline["style"] = styleOutline.str();
  paramsOutline["class"] = "transit-edge-outline " + getLineClass(line.id());

  std::stringstream styleStr;
  styleStr << "fill:none;stroke:#" << line.color() << ";" << css;

  if (!endMarker.empty()) {
    styleStr << ";marker-end:url(#" << endMarker << ")";
  }

  styleStr << ";stroke-linecap:round;stroke-opacity:1;stroke-width:"
           << width;
  Params params;
  params["style"] = styleStr.str();
  params["class"] = "transit-edge " + getLineClass(line.id());

  _delegates[0].insert(_delegates[0].begin(),
                       OutlinePrintPair(PrintDelegate(params, p),
                                        PrintDelegate(paramsOutline, p)));
}

// _____________________________________________________________________________
void SvgRenderer::renderEdgeTripGeom(const RenderGraph& outG,
                                     const shared::linegraph::LineEdge* e,
                                     const RenderParams& rparams) {
  UNUSED(rparams);
  const shared::linegraph::NodeFront* nfTo = e->getTo()->pl().frontFor(e);
  const shared::linegraph::NodeFront* nfFrom = e->getFrom()->pl().frontFor(e);

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
    const auto& lo = e->pl().lineOccAtPos(i);

    const Line* line = lo.line;
    PolyLine<double> p = center;

    if (p.getLength() < 0.01) continue;

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

    if (_cfg->renderDirMarkers && lo.direction != 0 &&
        center.getLength() > arrowLength * 3) {
      std::stringstream markerName;
      markerName << e << ":" << line << ":" << i;

      std::string markerPathMale = getMarkerPathMale(lineW);
      EndMarker emm(markerName.str() + "_m", "white", markerPathMale, lineW,
                    lineW);

      _markers.push_back(emm);

      PolyLine<double> firstPart = p.getSegmentAtDist(0, p.getLength() / 2);
      PolyLine<double> secondPart =
          p.getSegmentAtDist(p.getLength() / 2, p.getLength());

      if (lo.direction == e->getTo()) {
        renderLinePart(firstPart, lineW, *line, css, oCss,
                       markerName.str() + "_m");
        renderLinePart(secondPart.reversed(), lineW, *line, css, oCss);
      } else {
        renderLinePart(secondPart.reversed(), lineW, *line, css, oCss,
                       markerName.str() + "_m");
        renderLinePart(firstPart, lineW, *line, css, oCss);
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
void SvgRenderer::renderDelegates(const RenderGraph& outG,
                                  const RenderParams& rparams) {
  UNUSED(outG);
  for (auto& a : _delegates) {
    _w.openTag("g");
    for (auto& pd : a.second) {
      if (_cfg->outlineWidth > 0) {
        printLine(pd.back.second, pd.back.first, rparams);
      }
      printLine(pd.front.second, pd.front.first, rparams);
    }
    _w.closeTag();
  }

  for (auto& a : _innerDelegates) {
    _w.openTag("g");
    for (auto& b : a) {
      for (auto& pd : b.second) {
        if (_cfg->outlineWidth > 0) {
          printLine(pd.back.second, pd.back.first, rparams);
        }
      }
      for (auto& pd : b.second) {
        printLine(pd.front.second, pd.front.first, rparams);
      }
    }
    _w.closeTag();
  }
}

// _____________________________________________________________________________
void SvgRenderer::printPoint(const DPoint& p, const std::string& style,
                             const RenderParams& rparams) {
  UNUSED(rparams);
  std::map<std::string, std::string> params;
  params["cx"] = std::to_string(p.getX());
  params["cy"] = std::to_string(p.getY());
  params["r"] = "2";
  params["style"] = style;
  _w.openTag("circle", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::printLine(const PolyLine<double>& l, const std::string& style,
                            const RenderParams& rparams) {
  std::map<std::string, std::string> params;
  params["style"] = style;
  printLine(l, params, rparams);
}

// _____________________________________________________________________________
void SvgRenderer::printLine(const PolyLine<double>& l,
                            const std::map<std::string, std::string>& ps,
                            const RenderParams& rparams) {
  UNUSED(rparams);
  std::map<std::string, std::string> params = ps;
  std::stringstream points;
  points << std::setprecision(15);

  for (auto& p : l.getLine()) {
    points << " " << p.getX() << "," << p.getY();
  }

  params["points"] = points.str();

  _w.openTag("polyline", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::printPolygon(const Polygon<double>& g,
                               const std::map<std::string, std::string>& ps,
                               const RenderParams& rparams) {
  UNUSED(rparams);
  std::map<std::string, std::string> params = ps;
  std::stringstream points;
  points << std::setprecision(15);

  for (auto& p : g.getOuter()) {
    points << " " << p.getX() << "," << p.getY();
  }

  params["points"] = points.str();
  params["class"] = "station-poly";

  _w.openTag("polygon", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::printCircle(const DPoint& center, double rad,
                              const std::string& style,
                              const RenderParams& rparams) {
  std::map<std::string, std::string> params;
  params["style"] = style;
  printCircle(center, rad, params, rparams);
}

// _____________________________________________________________________________
void SvgRenderer::printCircle(const DPoint& center, double rad,
                              const std::map<std::string, std::string>& ps,
                              const RenderParams& rparams) {
  UNUSED(rparams);
  std::map<std::string, std::string> params = ps;
  std::stringstream points;
  points << std::setprecision(15);

  params["cx"] = std::to_string(center.getX());
  params["cy"] = std::to_string(center.getY());
  params["r"] = std::to_string(rad);

  _w.openTag("circle", params);
  _w.closeTag();
}

// _____________________________________________________________________________
size_t InnerClique::getNumBranchesIn(
    const shared::linegraph::LineEdge* edg) const {
  std::set<size_t> slots;
  size_t ret = 0;
  for (const auto& ig : geoms) {
    if (ig.from.edge == edg && !slots.insert(ig.slotFrom).second) ret++;
    if (ig.to.edge == edg && !slots.insert(ig.slotTo).second) ret++;
  }

  return ret;
}

// _____________________________________________________________________________
void SvgRenderer::renderStationLabels(const Labeller& labeller,
                                      const RenderParams& rparams) {
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
    points << std::setprecision(15);
    std::map<std::string, std::string> pathPars;

    points << "M" << textPath.front().getX() << " "
           << flippedY(rparams, textPath.front().getY());

    for (auto& p : textPath.getLine()) {
      points << " L" << p.getX() << " " << flippedY(rparams, p.getY());
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
    params["font-family"] = "Ubuntu Condensed";
    params["dy"] = shift;
    params["font-size"] = util::toString(label.fontSize);

    _w.openTag("text", params);
    _w.openTag("textPath", {{"dy", shift},
                            {"xlink:href", "#" + idStr},
                            {"startOffset", startOffset},
                            {"text-anchor", textAnchor}});

    if (label.labelLines.size() > 1 &&
        _cfg->showMergedStopMembers == "multiline") {
      for (size_t i = 0; i < label.labelLines.size(); ++i) {
        _w.openTag("tspan",
                   {{"x", "0"},
                    {"dy", i == 0 ? "0" : util::toString(label.fontSize)}});
        _w.writeText(label.labelLines[i]);
        _w.closeTag();
      }
    } else {
      _w.writeText(label.labelText.empty() ? label.s.name : label.labelText);
    }
    _w.closeTag();
    _w.closeTag();
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::renderLineLabels(const Labeller& labeller,
                                   const RenderParams& rparams) {
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
    points << std::setprecision(15);
    std::map<std::string, std::string> pathPars;

    points << "M" << textPath.front().getX() << " "
           << flippedY(rparams, textPath.front().getY());

    for (auto& p : textPath.getLine()) {
      points << " L" << p.getX() << " " << flippedY(rparams, p.getY());
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
    params["font-family"] = "Ubuntu";
    params["dy"] = shift;
    params["font-size"] = util::toString(label.fontSize);

    _w.openTag("text", params);
    _w.openTag("textPath", {{"dy", shift},
                            {"xlink:href", "#" + idStr},
                            {"text-anchor", "middle"},
                            {"startOffset", "50%"}});

    double dy = 0;
    for (auto line : label.lines) {
      _w.openTag("tspan",
                 {{"fill", "#" + line->color()}, {"dx", util::toString(dy)}});
      dy = label.fontSize / 3;
      _w.writeText(line->label());
      _w.closeTag();
    }
    _w.closeTag();
    _w.closeTag();
  }
  _w.closeTag();
}

// _____________________________________________________________________________
double InnerClique::getZWeight() const {
  // more weight = more to the bottom

  double BRANCH_WEIGHT = 4;

  double ret = 0;

  ret = geoms.size();  // baseline: threads with more lines to the bottom,
                       // because they are easier to follow

  for (const auto& nf : n->pl().fronts()) {
    ret -= getNumBranchesIn(nf.edge) * BRANCH_WEIGHT;
  }

  return ret;
}

// _____________________________________________________________________________
std::string SvgRenderer::getLineClass(const std::string& id) const {
  auto i = lineClassIds.find(id);
  if (i != lineClassIds.end()) return "line-" + std::to_string(i->second);

  lineClassIds[id] = ++lineClassId;
  return "line-" + std::to_string(lineClassId);
}

// _____________________________________________________________________________
bool InnerClique::operator<(const InnerClique& rhs) const {
  // more weight = more to the bottom
  return getZWeight() > rhs.getZWeight();
}

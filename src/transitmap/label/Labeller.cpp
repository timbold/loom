// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#ifdef LOOM_HAVE_FREETYPE
#ifdef Max
#pragma push_macro("Max")
#undef Max
#endif
#ifdef Min
#pragma push_macro("Min")
#undef Min
#endif
#include <ft2build.h>
#include FT_FREETYPE_H
#include "transitmap/label/tt_norms_pro_regular.h"
#ifdef Min
#pragma pop_macro("Min")
#endif
#ifdef Max
#pragma pop_macro("Max")
#endif
#endif

#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/label/Labeller.h"
#include "util/geo/Geo.h"

using shared::rendergraph::RenderGraph;
using transitmapper::label::Labeller;
using transitmapper::label::LineLabel;
using transitmapper::label::Overlaps;
using transitmapper::label::StationLabel;

using util::geo::MultiLine;
using util::geo::PolyLine;

namespace {

// Penalty for placing terminus labels at non horizontal/vertical angles.
constexpr double kTerminusAnglePen = 3.0;

// Decode UTF-8 string into Unicode code points.
std::vector<char32_t> decodeUtf8(const std::string &s) {
  std::vector<char32_t> cps;
  for (size_t i = 0; i < s.size();) {
    unsigned char c = static_cast<unsigned char>(s[i]);
    char32_t cp = 0;
    if (c < 0x80) {
      cp = c;
      i += 1;
    } else if ((c >> 5) == 0x6 && i + 1 < s.size()) {
      cp = ((c & 0x1F) << 6) | (static_cast<unsigned char>(s[i + 1]) & 0x3F);
      i += 2;
    } else if ((c >> 4) == 0xE && i + 2 < s.size()) {
      cp = ((c & 0x0F) << 12) |
           ((static_cast<unsigned char>(s[i + 1]) & 0x3F) << 6) |
           (static_cast<unsigned char>(s[i + 2]) & 0x3F);
      i += 3;
    } else if ((c >> 3) == 0x1E && i + 3 < s.size()) {
      cp = ((c & 0x07) << 18) |
           ((static_cast<unsigned char>(s[i + 1]) & 0x3F) << 12) |
           ((static_cast<unsigned char>(s[i + 2]) & 0x3F) << 6) |
           (static_cast<unsigned char>(s[i + 3]) & 0x3F);
      i += 4;
    } else {
      // Invalid UTF-8 byte, skip.
      i += 1;
      continue;
    }
    cps.push_back(cp);
  }
  return cps;
}

// Encode Unicode code points back to a UTF-8 string.
std::string encodeUtf8(const std::vector<char32_t> &cps, size_t start,
                       size_t end) {
  std::string out;
  for (size_t i = start; i < end; ++i) {
    char32_t cp = cps[i];
    if (cp <= 0x7F) {
      out.push_back(static_cast<char>(cp));
    } else if (cp <= 0x7FF) {
      out.push_back(static_cast<char>(0xC0 | ((cp >> 6) & 0x1F)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    } else if (cp <= 0xFFFF) {
      out.push_back(static_cast<char>(0xE0 | ((cp >> 12) & 0x0F)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    } else {
      out.push_back(static_cast<char>(0xF0 | ((cp >> 18) & 0x07)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
      out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
    }
  }
  return out;
}

// Unicode-aware whitespace check based on Unicode whitespace characters.
bool isUnicodeWhitespace(char32_t c) {
  return c == 0x0009 || c == 0x000A || c == 0x000B || c == 0x000C ||
         c == 0x000D || c == 0x0020 || c == 0x0085 || c == 0x00A0 ||
         c == 0x1680 || (c >= 0x2000 && c <= 0x200A) || c == 0x2028 ||
         c == 0x2029 || c == 0x202F || c == 0x205F || c == 0x3000;
}

std::string trimCopy(const std::string &s) {
  std::vector<char32_t> cps = decodeUtf8(s);
  size_t start = 0;
  while (start < cps.size() && isUnicodeWhitespace(cps[start])) {
    ++start;
  }
  size_t end = cps.size();
  while (end > start && isUnicodeWhitespace(cps[end - 1])) {
    --end;
  }
  return encodeUtf8(cps, start, end);
}

double getTextWidthFT(const std::string &text, double fontSize,
                      double resolution) {

#ifdef LOOM_HAVE_FREETYPE
  static FT_Library library = nullptr;
  static FT_Face face = nullptr;
  static bool initialized = false;
  if (!initialized) {
    if (FT_Init_FreeType(&library) ||
        FT_New_Memory_Face(library, tt_norms_pro_regular_otf,
                           tt_norms_pro_regular_otf_len, 0, &face)) {
      return (text.size() + 1) * fontSize / 2.1;
    }
    initialized = true;
  }

  if (!face) {
    return (text.size() + 1) * fontSize / 2.1;
  }

  FT_Set_Pixel_Sizes(face, 0,
                     static_cast<FT_UInt>(std::round(fontSize * resolution)));

  double width = 0.0;
  FT_UInt prevIdx = 0;
  for (size_t i = 0; i < text.size();) {
    unsigned char c = static_cast<unsigned char>(text[i]);
    FT_ULong cp = 0;
    size_t extra = 0;
    if (c < 0x80) {
      cp = c;
      extra = 0;
    } else if ((c & 0xE0) == 0xC0 && i + 1 < text.size()) {
      cp = c & 0x1F;
      extra = 1;
    } else if ((c & 0xF0) == 0xE0 && i + 2 < text.size()) {
      cp = c & 0x0F;
      extra = 2;
    } else if ((c & 0xF8) == 0xF0 && i + 3 < text.size()) {
      cp = c & 0x07;
      extra = 3;
    } else {
      ++i;
      continue;
    }

    bool invalid = false;
    for (size_t j = 1; j <= extra; ++j) {
      unsigned char cc = static_cast<unsigned char>(text[i + j]);
      if ((cc & 0xC0) != 0x80) {
        invalid = true;
        break;
      }
      cp = (cp << 6) | (cc & 0x3F);
    }
    if (invalid) {
      ++i;
      continue;
    }
    i += extra + 1;

    FT_UInt glyphIdx = FT_Get_Char_Index(face, cp);
    if (FT_Load_Glyph(face, glyphIdx, FT_LOAD_DEFAULT))
      continue;
    if (prevIdx) {
      FT_Vector delta;
      if (!FT_Get_Kerning(face, prevIdx, glyphIdx, FT_KERNING_DEFAULT,
                          &delta)) {
        width += static_cast<double>(delta.x) / 64.0;
      }
    }
    width += static_cast<double>(face->glyph->advance.x) / 64.0;
    prevIdx = glyphIdx;
  }
  if (prevIdx) {
    width -= static_cast<double>(face->glyph->advance.x -
                                 (face->glyph->metrics.horiBearingX +
                                  face->glyph->metrics.width)) /
             64.0;
  }

  return width / resolution;
#else
  (void)resolution;
  return (text.size() + 1) * fontSize / 2.1;
#endif
}
} // namespace

// _____________________________________________________________________________
Labeller::Labeller(const config::Config *cfg) : _cfg(cfg) {}

// _____________________________________________________________________________
void Labeller::label(const RenderGraph &g, bool notDeg2) {
  labelStations(g, notDeg2);
  labelLines(g);
}

// _____________________________________________________________________________
util::geo::MultiLine<double>
Labeller::getStationLblBand(const shared::linegraph::LineNode *n,
                            double fontSize, uint8_t offset,
                            const RenderGraph &g) {
  // TODO: the hull padding should be the same as in the renderer
  auto statHull = g.getStopGeoms(n, _cfg->tightStations, 4);

  double rad = util::geo::getEnclosingRadius(*n->pl().getGeom(), statHull);

  // measure the label width using FreeType
  std::string lbl = trimCopy(n->pl().stops().front().name);
  double textWidth = getTextWidthFT(lbl, fontSize, _cfg->outputResolution);
  double spaceWidth = getTextWidthFT("___", fontSize, _cfg->outputResolution);
  double offsetW = _cfg->lineSpacing + _cfg->lineWidth;
  double labelW = offsetW + textWidth + spaceWidth;

  util::geo::MultiLine<double> band;

  // TODO: should also be determined based on the font
  double h = fontSize * 0.75;

  util::geo::Line<double> geomBaseLine, geomMiddle, geomTop, capLeft, capRight;
  geomBaseLine.push_back({n->pl().getGeom()->getX() + rad + offsetW,
                          n->pl().getGeom()->getY() - (offset * h / 2)});
  geomBaseLine.push_back({n->pl().getGeom()->getX() + rad + labelW,
                          n->pl().getGeom()->getY() - (offset * h / 2)});

  geomMiddle.push_back({n->pl().getGeom()->getX() + rad + offsetW,
                        n->pl().getGeom()->getY() + h / 2 - (offset * h / 2)});
  geomMiddle.push_back({n->pl().getGeom()->getX() + rad + labelW,
                        n->pl().getGeom()->getY() + h / 2 - (offset * h / 2)});

  geomTop.push_back({n->pl().getGeom()->getX() + rad + offsetW,
                     n->pl().getGeom()->getY() + h - (offset * h / 2)});
  geomTop.push_back({n->pl().getGeom()->getX() + rad + labelW,
                     n->pl().getGeom()->getY() + h - (offset * h / 2)});

  capLeft = util::geo::PolyLine<double>(geomMiddle)
                .getOrthoLineAtDist(util::geo::len(geomMiddle), h)
                .getLine();
  capRight = util::geo::PolyLine<double>(geomMiddle)
                 .getOrthoLineAtDist(0, h)
                 .getLine();

  band.push_back(geomBaseLine);
  band.push_back(geomMiddle);
  band.push_back(geomTop);
  band.push_back(capLeft);
  band.push_back(capRight);

  return band;
}

// _____________________________________________________________________________
void Labeller::labelStations(const RenderGraph &g, bool notdeg2) {
  // Partition nodes so that termini are processed first. This allows the
  // overlap logic to account for already placed terminus labels when placing
  // regular station labels.
  std::vector<const shared::linegraph::LineNode *> termini;
  std::vector<const shared::linegraph::LineNode *> others;
  for (auto n : g.getNds()) {
    if (n->pl().stops().size() == 0 || (notdeg2 && n->getDeg() == 2))
      continue;
    if (g.isTerminus(n)) {
      termini.push_back(n);
    } else {
      others.push_back(n);
    }
  }

  std::sort(termini.begin(), termini.end(), statNdCmp);
  std::sort(others.begin(), others.end(), statNdCmp);

  std::vector<const shared::linegraph::LineNode *> orderedNds;
  orderedNds.reserve(termini.size() + others.size());
  orderedNds.insert(orderedNds.end(), termini.begin(), termini.end());
  orderedNds.insert(orderedNds.end(), others.begin(), others.end());
  auto mapBox = g.getBBox();

  for (auto n : orderedNds) {
    double fontSize = _cfg->stationLabelSize;
    bool isTerminus = g.isTerminus(n);
    if (_cfg->highlightTerminals && isTerminus) {
      fontSize += 10;
    }
    if (_cfg->fontSvgMax >= 0 &&
        fontSize * _cfg->outputResolution > _cfg->fontSvgMax) {
      fontSize = _cfg->fontSvgMax / _cfg->outputResolution;
    }
    int prefDeg = 0;
    if (n->pl().stops().size()) {
      const auto &sp = n->pl().stops().front().pos;
      const auto *cp = n->pl().getGeom();
      double dx = sp.getX() - cp->getX();
      double dy = sp.getY() - cp->getY();
      if (std::abs(dx) > 1e-9 || std::abs(dy) > 1e-9) {
        double ang = std::atan2(dy, dx) * 180.0 / M_PI;
        prefDeg = static_cast<int>(std::round(ang / 30.0));
        prefDeg = (prefDeg % 12 + 12) % 12;
      }
    }

    auto station = n->pl().stops().front();
    station.name = trimCopy(station.name);

    std::vector<StationLabel> cands;

    for (uint8_t offset = 0; offset < 3; offset++) {
      for (size_t deg = 0; deg < 12; deg++) {
        // generate candidate bands using the (possibly boosted) font size
        auto band = getStationLblBand(n, fontSize, offset, g);
        band = util::geo::rotate(band, 30 * deg, *n->pl().getGeom());

        auto box = util::geo::getBoundingBox(band);
        double diag = util::geo::dist(box.getLowerLeft(), box.getUpperRight());
        double searchRad =
            g.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing) +
            std::max(_cfg->stationLabelSize, diag);

        auto overlaps = getOverlaps(band, n, g, searchRad);
        if (overlaps.lineOverlaps + overlaps.statLabelOverlaps +
                overlaps.statOverlaps + overlaps.landmarkOverlaps >
            0)
          continue;

        // measure local crowding to discourage labels in dense regions
        auto neighEdges = g.getNeighborEdges(band[0], searchRad);
        std::set<const shared::linegraph::LineNode*> neighNodes;
        for (auto e : neighEdges) {
          neighNodes.insert(e->getFrom());
          neighNodes.insert(e->getTo());
        }
        double clusterPen =
            static_cast<double>(neighEdges.size() + neighNodes.size());

        bool outside = box.getLowerLeft().getX() < mapBox.getLowerLeft().getX() ||
                        box.getLowerLeft().getY() < mapBox.getLowerLeft().getY() ||
                        box.getUpperRight().getX() > mapBox.getUpperRight().getX() ||
                        box.getUpperRight().getY() > mapBox.getUpperRight().getY();
        double outsidePen = outside ? -5.0 : 0.0;

        size_t diff = (deg + 12 - prefDeg) % 12;
        if (diff > 6)
          diff = 12 - diff;
        double sidePen = static_cast<double>(diff) * 5.0;
        double termPen = isTerminus && (deg % 6 != 0) && (deg % 6 != 3)
                             ? kTerminusAnglePen
                             : 0;
        cands.emplace_back(PolyLine<double>(band[0]), band, fontSize,
                           isTerminus, deg, offset, overlaps, sidePen + termPen,
                           _cfg->stationLineOverlapPenalty, clusterPen,
                           outsidePen, station);
      }
    }

    std::sort(cands.begin(), cands.end());
    if (cands.size() == 0)
      continue;
    auto cand = cands.front();
    // Recompute band and geometry in case the font size was adjusted.
    cand.band = util::geo::rotate(
        getStationLblBand(n, cand.fontSize, static_cast<uint8_t>(cand.pos), g),
        30 * cand.deg, *n->pl().getGeom());
    cand.geom = PolyLine<double>(cand.band[0]);

    auto* nn = const_cast<shared::linegraph::LineNode*>(n);
    if (!nn->pl().stops().empty()) {
      nn->pl().stops()[0].labelDeg = cand.deg;
      cand.s.labelDeg = cand.deg;
    }

    _stationLabels.push_back(cand);
    _statLblIdx.add(cand.band, _stationLabels.size() - 1);
  }
}

// _____________________________________________________________________________
Overlaps Labeller::getOverlaps(const util::geo::MultiLine<double> &band,
                               const shared::linegraph::LineNode *forNd,
                               const RenderGraph &g, double radius) const {
  std::set<const shared::linegraph::LineEdge *> proced;

  Overlaps ret{0, 0, 0, 0, 0, 0};

  std::set<const shared::linegraph::LineNode *> procedNds{forNd};

  for (auto line : band) {
    auto neighs = g.getNeighborEdges(line, radius);
    for (auto neigh : neighs) {
      if (proced.count(neigh))
        continue;

      if (util::geo::dist(*neigh->pl().getGeom(), band) <
          g.getTotalWidth(neigh) / 2) {
        ret.lineOverlaps++;
      }
      proced.insert(neigh);

      std::vector<const shared::linegraph::LineNode *> nds = {neigh->getTo(),
                                                              neigh->getFrom()};

      for (auto nd : nds) {
        if (nd->pl().stops().size() && !procedNds.count(nd)) {
          procedNds.insert(nd);

          auto statHull = g.getStopGeoms(nd, _cfg->tightStations, 4);

          double rad =
              util::geo::getEnclosingRadius(*nd->pl().getGeom(), statHull);

          if (util::geo::dist(*nd->pl().getGeom(), band) <
              rad + (_cfg->lineWidth + _cfg->lineSpacing) / 2) {
            ret.statOverlaps++;
          }
        }
      }
    }
  }

  std::set<size_t> labelNeighs;
  _statLblIdx.get(band, radius, &labelNeighs);

  for (auto id : labelNeighs) {
    auto labelNeigh = _stationLabels[id];
    if (util::geo::dist(labelNeigh.band, band) < 1) {
      if (labelNeigh.bold)
        ret.termLabelOverlaps++;
      else
        ret.statLabelOverlaps++;
    }
  }

  std::set<size_t> landmarkNeighs;
  _landmarkIdx.get(band, radius, &landmarkNeighs);
  ret.landmarkOverlaps += landmarkNeighs.size();

  return ret;
}

// _____________________________________________________________________________
void Labeller::labelLines(const RenderGraph &g) {
  LineLblIdx labelIdx = LineLblIdx();
  for (auto n : g.getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n)
        continue;
      double geomLen = util::geo::len(*e->pl().getGeom());

      // estimate label width using precise text measurement
      double fontSize = _cfg->lineLabelSize;
      double spacing = fontSize / 3.0;
      double labelW = 0.0;
      bool first = true;
      for (auto lo : e->pl().getLines()) {
        if (!first) {
          labelW += spacing;
        }
        labelW +=
            getTextWidthFT(lo.line->label(), fontSize, _cfg->outputResolution);
        first = false;
      }

      // try out positions
      double step = fontSize;

      std::vector<LineLabel> cands;

      for (int dir = -1; dir < 2; dir += 2) {
        double start = 0;
        while (start + labelW <= geomLen) {
          PolyLine<double> cand(util::geo::segment(
              *e->pl().getGeom(), start / geomLen, (start + labelW) / geomLen));
          if (cand.getLength() < labelW) {
            start += step;
            continue;
          }
          if (cand.getLength() < 5)
            break;
          cand.offsetPerp(dir * (g.getTotalWidth(e) / 2 +
                                 (_cfg->lineSpacing + _cfg->lineWidth)));

          // Reject candidates that bend too much.
          auto line = cand.getLine();
          double chord = util::geo::dist(line.front(), line.back());
          if (chord > 0 &&
              cand.getLength() / chord > _cfg->lineLabelLengthRatio) {
            start += step;
            continue;
          }
          double maxBend = 0.0;
          if (line.size() >= 3) {
            for (size_t i = 1; i + 1 < line.size(); ++i) {
              auto a = line[i - 1];
              auto b = line[i];
              auto c = line[i + 1];
              double v1x = a.getX() - b.getX();
              double v1y = a.getY() - b.getY();
              double v2x = c.getX() - b.getX();
              double v2y = c.getY() - b.getY();
              double len1 = std::hypot(v1x, v1y);
              double len2 = std::hypot(v2x, v2y);
              if (len1 == 0 || len2 == 0)
                continue;
              double cosTheta = (v1x * v2x + v1y * v2y) / (len1 * len2);
              cosTheta = std::max(-1.0, std::min(1.0, cosTheta));
              double theta = std::acos(cosTheta);
              double bend = M_PI - theta;
              if (bend > maxBend)
                maxBend = bend;
            }
          }
          if (maxBend > _cfg->lineLabelBendAngle) {
            start += step;
            continue;
          }

          bool block = false;

          for (auto neigh : g.getNeighborEdges(
                   cand.getLine(),
                   g.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing) +
                       fontSize * 4)) {
            if (neigh == e)
              continue;
            if (util::geo::dist(cand.getLine(), *neigh->pl().getGeom()) <
                (g.getTotalWidth(neigh) / 2) + (fontSize)) {
              block = true;
              break;
            }
          }

          std::set<size_t> labelNeighs;
          _statLblIdx.get(MultiLine<double>{cand.getLine()},
                          g.getMaxLineNum() *
                              (_cfg->lineWidth + _cfg->lineSpacing),
                          &labelNeighs);

          for (auto neighId : labelNeighs) {
            auto neigh = _stationLabels[neighId];
            if (util::geo::dist(cand.getLine(), neigh.band) < (fontSize)) {
              block = true;
              break;
            }
          }

          if (dir < 0)
            cand.reverse();

          std::set<size_t> lineLabelNeighs;
          labelIdx.get(cand.getLine(),
                       20 * (_cfg->lineWidth + _cfg->lineSpacing),
                       &lineLabelNeighs);

          std::vector<const shared::linegraph::Line *> lines;
          for (auto lo : e->pl().getLines()) {
            lines.push_back(lo.line);
          }

          for (auto neighLabelId : lineLabelNeighs) {
            auto neighLabel = _lineLabels[neighLabelId];
            if (neighLabel.lines == lines &&
                util::geo::dist(cand.getLine(), neighLabel.geom.getLine()) <
                    20 * (_cfg->lineWidth + _cfg->lineSpacing)) {
              block = true;
              break;
            }
          }

          std::set<size_t> landmarkNeighs;
          _landmarkIdx.get(MultiLine<double>{cand.getLine()}, fontSize,
                           &landmarkNeighs);
          if (!landmarkNeighs.empty())
            block = true;

          if (!block)
            cands.push_back({cand, fabs((geomLen / 2) - (start + (labelW / 2))),
                             fontSize, lines});
          start += step;
        }
      }

      std::sort(cands.begin(), cands.end());
      if (cands.size() == 0)
        continue;
      _lineLabels.push_back(cands.front());
      labelIdx.add(cands.front().geom.getLine(), _lineLabels.size() - 1);
    }
  }
}

// _____________________________________________________________________________
const std::vector<LineLabel> &Labeller::getLineLabels() const {
  return _lineLabels;
}

// _____________________________________________________________________________
const std::vector<StationLabel> &Labeller::getStationLabels() const {
  return _stationLabels;
}

// _____________________________________________________________________________
std::vector<size_t> Labeller::getStationLabelDegrees() const {
  std::vector<size_t> ret;
  ret.reserve(_stationLabels.size());
  for (const auto &sl : _stationLabels) ret.push_back(sl.deg);
  return ret;
}

// _____________________________________________________________________________
util::geo::Box<double> Labeller::getBBox() const {
  util::geo::Box<double> ret;

  for (auto lbl : _lineLabels)
    ret = util::geo::extendBox(lbl.geom.getLine(), ret);
  for (auto lbl : _stationLabels)
    ret = util::geo::extendBox(lbl.band, ret);

  return ret;
}

// _____________________________________________________________________________
bool Labeller::collidesWithLabels(const util::geo::Box<double> &box) const {
  std::set<size_t> overlaps;
  _landmarkIdx.get(box, 0, &overlaps);
  if (!overlaps.empty())
    return true;
  _statLblIdx.get(box, 0, &overlaps);
  return !overlaps.empty();
}

// _____________________________________________________________________________
bool Labeller::addLandmark(const util::geo::Box<double> &box) {
  if (collidesWithLabels(box))
    return false;
  _landmarkIdx.add(box, _landmarks.size());
  _landmarks.push_back(box);
  return true;
}

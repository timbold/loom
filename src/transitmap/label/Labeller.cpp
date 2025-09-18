// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
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
#include "util/String.h"
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

// Number of candidate angles for station labels and their step size in degrees.
constexpr size_t kStationAngleSteps = 24;
constexpr double kStationAngleDeg = 15.0;

double pointToBoxDistance(const util::geo::DPoint &p,
                          const util::geo::Box<double> &box) {
  double minX = box.getLowerLeft().getX();
  double maxX = box.getUpperRight().getX();
  double minY = box.getLowerLeft().getY();
  double maxY = box.getUpperRight().getY();
  if (minX > maxX || minY > maxY) {
    return std::numeric_limits<double>::infinity();
  }

  double dx = 0.0;
  if (p.getX() < minX) {
    dx = minX - p.getX();
  } else if (p.getX() > maxX) {
    dx = p.getX() - maxX;
  }

  double dy = 0.0;
  if (p.getY() < minY) {
    dy = minY - p.getY();
  } else if (p.getY() > maxY) {
    dy = p.getY() - maxY;
  }

  if (dx == 0.0 && dy == 0.0) {
    return 0.0;
  }
  return std::sqrt(dx * dx + dy * dy);
}

std::pair<double, double> getLandmarkSizeMapUnits(
    const shared::rendergraph::Landmark &lm,
    const transitmapper::config::Config *cfg) {
  if (!cfg || cfg->outputResolution <= 0.0) {
    return {0.0, 0.0};
  }

  double maxWidthPx = cfg->stationLabelSize * cfg->outputResolution * 0.6 * 10.0;
  double widthPx = lm.size * cfg->outputResolution;
  double heightPx = widthPx;

  if (!lm.label.empty()) {
    double labelHeightPx = lm.fontSize;
    double labelWidthPx = util::toWStr(lm.label).size() * (labelHeightPx * 0.6);
    if (labelWidthPx > maxWidthPx && labelWidthPx > 0) {
      double factor = maxWidthPx / labelWidthPx;
      labelWidthPx = maxWidthPx;
      labelHeightPx *= factor;
    }
    widthPx = labelWidthPx;
    heightPx = labelHeightPx;
  } else {
    if (widthPx > maxWidthPx && widthPx > 0) {
      double factor = maxWidthPx / widthPx;
      widthPx = maxWidthPx;
      heightPx *= factor;
    }
  }

  return {widthPx / cfg->outputResolution, heightPx / cfg->outputResolution};
}

util::geo::Box<double> makeLandmarkBox(const shared::rendergraph::Landmark &lm,
                                       const transitmapper::config::Config *cfg) {
  auto dims = getLandmarkSizeMapUnits(lm, cfg);
  double halfW = dims.first / 2.0;
  double halfH = dims.second / 2.0;
  return util::geo::Box<double>(
      util::geo::DPoint(lm.coord.getX() - halfW, lm.coord.getY() - halfH),
      util::geo::DPoint(lm.coord.getX() + halfW, lm.coord.getY() + halfH));
}

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

struct NeighborSideInfo {
  const shared::linegraph::LineNode *neighbor = nullptr;
  int sign = 0;
  size_t multiplicity = 0;
};

struct StationLabelCandidate {
  StationLabel label;
  bool opposite;
  std::vector<NeighborSideInfo> neighborSides;
};

struct StationNodeCandidates {
  const shared::linegraph::LineNode *node = nullptr;
  std::vector<StationLabelCandidate> candidates;
};

struct StationCandidateConflict {
  size_t otherNodeIdx;
  size_t otherCandidateIdx;
  double penalty;
};

std::vector<size_t> optimizeStationLabelAssignments(
    const transitmapper::config::Config *cfg,
    const std::vector<StationNodeCandidates> &nodeCandidates,
    const std::vector<const shared::linegraph::LineNode *> &orderedNds) {
  size_t nodeCount = nodeCandidates.size();
  std::vector<size_t> assignment(nodeCount, 0);
  if (nodeCount == 0) {
    return assignment;
  }

  std::unordered_map<const shared::linegraph::LineNode *, size_t> nodeIndex;
  nodeIndex.reserve(nodeCount);
  for (size_t i = 0; i < nodeCount; ++i) {
    nodeIndex[nodeCandidates[i].node] = i;
  }

  std::unordered_map<const shared::linegraph::LineNode *, size_t> orderIndex;
  orderIndex.reserve(orderedNds.size());
  for (size_t i = 0; i < orderedNds.size(); ++i) {
    orderIndex[orderedNds[i]] = i;
  }

  std::vector<std::vector<double>> baseCosts(nodeCount);
  size_t totalCandidates = 0;
  for (size_t i = 0; i < nodeCount; ++i) {
    const auto &cands = nodeCandidates[i].candidates;
    baseCosts[i].reserve(cands.size());
    for (const auto &cand : cands) {
      baseCosts[i].push_back(cand.label.getPen());
    }
    totalCandidates += cands.size();
  }

  for (size_t i = 0; i < nodeCount; ++i) {
    size_t bestIdx = 0;
    double bestCost = std::numeric_limits<double>::infinity();
    for (size_t j = 0; j < baseCosts[i].size(); ++j) {
      if (baseCosts[i][j] < bestCost) {
        bestCost = baseCosts[i][j];
        bestIdx = j;
      }
    }
    assignment[i] = bestIdx;
  }

  util::geo::RTree<size_t, util::geo::MultiLine, double> candidateIdx;
  std::vector<std::pair<size_t, size_t>> idToCandidate;
  idToCandidate.reserve(totalCandidates);
  for (size_t i = 0; i < nodeCount; ++i) {
    for (size_t j = 0; j < nodeCandidates[i].candidates.size(); ++j) {
      size_t id = idToCandidate.size();
      idToCandidate.emplace_back(i, j);
      candidateIdx.add(nodeCandidates[i].candidates[j].label.band, id);
    }
  }

  std::vector<std::vector<std::vector<StationCandidateConflict>>> conflicts(nodeCount);
  for (size_t i = 0; i < nodeCount; ++i) {
    conflicts[i].resize(nodeCandidates[i].candidates.size());
  }

  for (size_t id = 0; id < idToCandidate.size(); ++id) {
    auto ref = idToCandidate[id];
    size_t nodeIdx = ref.first;
    size_t candIdx = ref.second;
    const auto &cand = nodeCandidates[nodeIdx].candidates[candIdx];
    std::set<size_t> neighborIds;
    candidateIdx.get<util::geo::MultiLine>(cand.label.band, 0.0, &neighborIds);
    for (auto neighborId : neighborIds) {
      if (neighborId <= id) continue;
      auto otherRef = idToCandidate[neighborId];
      size_t otherNodeIdx = otherRef.first;
      size_t otherCandIdx = otherRef.second;
      if (otherNodeIdx == nodeIdx) continue;
      const auto &otherCand = nodeCandidates[otherNodeIdx].candidates[otherCandIdx];
      if (util::geo::dist(cand.label.band, otherCand.label.band) >= 1.0) continue;
      size_t orderA = orderIndex[nodeCandidates[nodeIdx].node];
      size_t orderB = orderIndex[nodeCandidates[otherNodeIdx].node];
      size_t earlierNodeIdx = nodeIdx;
      size_t earlierCandIdx = candIdx;
      if (orderB < orderA || (orderB == orderA && otherNodeIdx < nodeIdx)) {
        earlierNodeIdx = otherNodeIdx;
        earlierCandIdx = otherCandIdx;
      }
      double pairPenalty =
          nodeCandidates[earlierNodeIdx].candidates[earlierCandIdx].label.bold
              ? 10.0
              : 20.0;
      conflicts[nodeIdx][candIdx].push_back({otherNodeIdx, otherCandIdx, pairPenalty});
      conflicts[otherNodeIdx][otherCandIdx].push_back({nodeIdx, candIdx, pairPenalty});
    }
  }

  double sameSidePenalty = cfg ? cfg->sameSidePenalty : 0.0;
  if (sameSidePenalty > 0.0) {
    for (size_t i = 0; i < nodeCount; ++i) {
      for (size_t ci = 0; ci < nodeCandidates[i].candidates.size(); ++ci) {
        const auto &cand = nodeCandidates[i].candidates[ci];
        for (const auto &nb : cand.neighborSides) {
          auto itNode = nodeIndex.find(nb.neighbor);
          if (itNode == nodeIndex.end()) continue;
          size_t otherIdx = itNode->second;
          if (otherIdx == i) continue;
          if (otherIdx < i) continue;
          double penaltyVal =
              sameSidePenalty * static_cast<double>(nb.multiplicity);
          if (penaltyVal <= 0.0) continue;
          for (size_t cj = 0; cj < nodeCandidates[otherIdx].candidates.size(); ++cj) {
            const auto &otherCand = nodeCandidates[otherIdx].candidates[cj];
            for (const auto &otherNb : otherCand.neighborSides) {
              if (otherNb.neighbor != nodeCandidates[i].node) continue;
              if (nb.sign * otherNb.sign < 0) {
                conflicts[i][ci].push_back({otherIdx, cj, penaltyVal});
                conflicts[otherIdx][cj].push_back({i, ci, penaltyVal});
              }
              break;
            }
          }
        }
      }
    }
  }

  std::vector<size_t> order(nodeCount);
  std::iota(order.begin(), order.end(), 0);
  std::mt19937 rng(42);

  for (size_t iter = 0; iter < 10; ++iter) {
    bool improved = false;
    std::shuffle(order.begin(), order.end(), rng);
    for (size_t idx : order) {
      size_t current = assignment[idx];
      double bestDelta = 0.0;
      size_t bestCand = current;
      for (size_t candIdx = 0; candIdx < nodeCandidates[idx].candidates.size(); ++candIdx) {
        if (candIdx == current) continue;
        double delta = baseCosts[idx][candIdx] - baseCosts[idx][current];
        for (const auto &conf : conflicts[idx][current]) {
          if (assignment[conf.otherNodeIdx] == conf.otherCandidateIdx) {
            delta -= conf.penalty;
          }
        }
        for (const auto &conf : conflicts[idx][candIdx]) {
          if (assignment[conf.otherNodeIdx] == conf.otherCandidateIdx) {
            delta += conf.penalty;
          }
        }
        if (delta < bestDelta - 1e-9) {
          bestDelta = delta;
          bestCand = candIdx;
        }
      }
      if (bestCand != current) {
        assignment[idx] = bestCand;
        improved = true;
      }
    }
    if (!improved) break;
  }

  return assignment;
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

Labeller::StationCrowdContext Labeller::computeStationFarCrowd(
    const util::geo::MultiLine<double> &band,
    const shared::linegraph::LineNode *stationNode, double searchRadius,
    const RenderGraph &g) const {
  StationCrowdContext ctx;
  if (band.empty()) {
    return ctx;
  }

  ctx.neighborEdges = g.getNeighborEdges(band[0], searchRadius);
  for (auto edge : ctx.neighborEdges) {
    if (!edge) continue;
    ctx.neighborNodes.insert(edge->getFrom());
    ctx.neighborNodes.insert(edge->getTo());
  }

  if (_cfg->stationLabelFarCrowdRadius <= 0 || band.size() <= 1 ||
      band[1].empty() || !stationNode) {
    return ctx;
  }

  auto stationPos = *stationNode->pl().getGeom();
  auto farPoint = band[1].front();
  double maxDist = util::geo::dist(farPoint, stationPos);
  for (const auto &p : band[1]) {
    double dist = util::geo::dist(p, stationPos);
    if (dist > maxDist) {
      maxDist = dist;
      farPoint = p;
    }
  }

  double radius = _cfg->stationLabelFarCrowdRadius;
  int farCrowdCount = 0;

  for (auto edge : ctx.neighborEdges) {
    if (!edge || edge->getFrom() == stationNode || edge->getTo() == stationNode)
      continue;
    double width = g.getTotalWidth(edge) / 2.0;
    double dist = util::geo::dist(farPoint, *edge->pl().getGeom());
    if (dist <= radius + width) {
      ++farCrowdCount;
    }
  }

  std::set<size_t> nearbyLabels;
  _statLblIdx.get(farPoint, radius, &nearbyLabels);
  for (auto labelIdx : nearbyLabels) {
    if (labelIdx >= _stationLabels.size()) continue;
    const auto &label = _stationLabels[labelIdx];
    double dist = util::geo::dist(farPoint, label.band);
    if (dist <= radius) {
      ++farCrowdCount;
    }
  }

  std::set<const shared::linegraph::LineNode *> processedStations;
  for (auto node : ctx.neighborNodes) {
    if (!node || node == stationNode) continue;
    if (!processedStations.insert(node).second) continue;
    if (node->pl().stops().empty()) continue;
    auto hulls = g.getStopGeoms(node, _cfg->tightStations, 4);
    for (const auto &hull : hulls) {
      if (util::geo::dist(farPoint, hull) <= radius) {
        ++farCrowdCount;
        break;
      }
    }
  }

  const auto &landmarks = g.getLandmarks();
  for (const auto &lm : landmarks) {
    auto box = makeLandmarkBox(lm, _cfg);
    double dist = pointToBoxDistance(farPoint, box);
    if (dist <= radius) {
      ++farCrowdCount;
    }
  }

  if (farCrowdCount > 0) {
    ctx.farCrowdPen =
        static_cast<double>(farCrowdCount) * _cfg->stationLabelFarCrowdPenalty;
  }

  return ctx;
}

// _____________________________________________________________________________
void Labeller::labelStations(const RenderGraph &g, bool notdeg2) {
  _stationLabels.clear();
  _statLblNodes.clear();
  _statLblIdx = StatLblIdx();
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

  std::unordered_map<const shared::linegraph::LineNode*,
                     std::vector<std::pair<const shared::linegraph::LineNode*, int>>>
      sidePrefs;
  std::unordered_map<const shared::linegraph::LineNode*, size_t> labelIndex;
  for (auto n : orderedNds) {
    if (!n || n->pl().stops().empty()) continue;
    auto *mutableNode = const_cast<shared::linegraph::LineNode *>(n);
    mutableNode->pl().stops()[0].labelDeg =
        std::numeric_limits<size_t>::max();
  }

  std::vector<StationNodeCandidates> nodeCandidates;
  nodeCandidates.reserve(orderedNds.size());

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
        prefDeg = static_cast<int>(std::round(ang / kStationAngleDeg));
        prefDeg = (prefDeg % static_cast<int>(kStationAngleSteps) +
                   static_cast<int>(kStationAngleSteps)) %
                  static_cast<int>(kStationAngleSteps);
      }
    }

    auto station = n->pl().stops().front();
    station.name = trimCopy(station.name);

    std::vector<StationLabelCandidate> cands;

    for (uint8_t offset = 0; offset < 3; offset++) {
      for (size_t deg = 0; deg < kStationAngleSteps; deg++) {
        // generate candidate bands using the (possibly boosted) font size
        auto band = getStationLblBand(n, fontSize, offset, g);
        band =
            util::geo::rotate(band, kStationAngleDeg * deg, *n->pl().getGeom());

        auto box = util::geo::getBoundingBox(band);
        double diag = util::geo::dist(box.getLowerLeft(), box.getUpperRight());
        double searchRad =
            g.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing) +
            std::max(_cfg->stationLabelSize, diag);

        auto overlaps = getOverlaps(band, n, g, searchRad);

        auto crowd = computeStationFarCrowd(band, n, searchRad, g);
        const auto &neighEdges = crowd.neighborEdges;
        const auto &neighNodes = crowd.neighborNodes;
        double farCrowdPen = crowd.farCrowdPen;
        double area = (box.getUpperRight().getX() - box.getLowerLeft().getX()) *
                      (box.getUpperRight().getY() - box.getLowerLeft().getY());
        double neighborCount =
            static_cast<double>(neighEdges.size() + neighNodes.size());
        double clusterPen =
            area > 0.0 ? neighborCount / area : neighborCount;

        bool outside =
            box.getLowerLeft().getX() < mapBox.getLowerLeft().getX() ||
            box.getLowerLeft().getY() < mapBox.getLowerLeft().getY() ||
            box.getUpperRight().getX() > mapBox.getUpperRight().getX() ||
            box.getUpperRight().getY() > mapBox.getUpperRight().getY();

        size_t diff = (deg + kStationAngleSteps - prefDeg) % kStationAngleSteps;
        if (diff > kStationAngleSteps / 2)
          diff = kStationAngleSteps - diff;
        double sidePen =
            static_cast<double>(diff) * _cfg->sidePenaltyWeight;
        double termPen = isTerminus && (deg % (kStationAngleSteps / 4) != 0)
                             ? kTerminusAnglePen
                             : 0;

        double candAng = deg * kStationAngleDeg * M_PI / 180.0;
        double candVecX = std::cos(candAng);
        double candVecY = std::sin(candAng);
        double sameSidePen = 0.0;
        std::unordered_map<const shared::linegraph::LineNode *, NeighborSideInfo>
            neighborInfo;
        for (auto e : n->getAdjList()) {
          auto neigh = e->getFrom() == n ? e->getTo() : e->getFrom();
          if (neigh->pl().stops().empty())
            continue;
          double edgeVecX =
              neigh->pl().getGeom()->getX() - n->pl().getGeom()->getX();
          double edgeVecY =
              neigh->pl().getGeom()->getY() - n->pl().getGeom()->getY();
          double candSide = edgeVecX * candVecY - edgeVecY * candVecX;
          auto &info = neighborInfo[neigh];
          info.neighbor = neigh;
          info.sign = candSide >= 0 ? 1 : -1;
          info.multiplicity += 1;
          size_t neighDeg = neigh->pl().stops()[0].labelDeg;
          if (neighDeg == std::numeric_limits<size_t>::max())
            continue;
          double neighAng = neighDeg * kStationAngleDeg * M_PI / 180.0;
          double neighVecX = std::cos(neighAng);
          double neighVecY = std::sin(neighAng);
          double neighSide = edgeVecX * neighVecY - edgeVecY * neighVecX;
          if (candSide * neighSide < 0)
            sameSidePen += _cfg->sameSidePenalty;
        }

        bool opposite = sameSidePen > 0.0;
        StationLabelCandidate cand{
            StationLabel(PolyLine<double>(band[0]), band, fontSize, isTerminus,
                         deg, offset, overlaps,
                         sidePen + termPen + sameSidePen,
                         _cfg->stationLineOverlapPenalty, clusterPen,
                         farCrowdPen, outside,
                         _cfg->clusterPenScale, _cfg->outsidePenalty,
                         &_cfg->orientationPenalties, station),
            opposite,
            {}};
        cand.neighborSides.reserve(neighborInfo.size());
        for (auto &entry : neighborInfo) {
          cand.neighborSides.push_back(entry.second);
        }
        cands.push_back(std::move(cand));
      }
    }

    bool hasSame =
        std::any_of(cands.begin(), cands.end(), [](const StationLabelCandidate &c) {
          return !c.opposite;
        });
    std::vector<StationLabelCandidate> filtered;
    filtered.reserve(cands.size());
    for (const auto &c : cands) {
      if (!hasSame || !c.opposite) filtered.push_back(c);
    }

    std::sort(filtered.begin(), filtered.end(),
              [](const StationLabelCandidate &a,
                 const StationLabelCandidate &b) { return a.label < b.label; });
    if (filtered.empty()) continue;

    StationNodeCandidates entry;
    entry.node = n;
    entry.candidates = std::move(filtered);
    nodeCandidates.push_back(std::move(entry));
  }

  auto assignment =
      optimizeStationLabelAssignments(_cfg, nodeCandidates, orderedNds);

  for (size_t idx = 0; idx < nodeCandidates.size(); ++idx) {
    const auto *node = nodeCandidates[idx].node;
    if (!node) continue;
    if (idx >= assignment.size()) continue;
    size_t candIdx = assignment[idx];
    if (candIdx >= nodeCandidates[idx].candidates.size()) continue;

    StationLabel cand = nodeCandidates[idx].candidates[candIdx].label;
    cand.band = util::geo::rotate(
        getStationLblBand(node, cand.fontSize, static_cast<uint8_t>(cand.pos), g),
        kStationAngleDeg * cand.deg, *node->pl().getGeom());
    cand.geom = PolyLine<double>(cand.band[0]);

    _stationLabels.push_back(cand);
    _statLblNodes.push_back(node);
    size_t storedIdx = _stationLabels.size() - 1;
    _statLblIdx.add(cand.band, storedIdx);
    labelIndex[node] = storedIdx;

    auto *mutableNode = const_cast<shared::linegraph::LineNode *>(node);
    if (!mutableNode->pl().stops().empty()) {
      mutableNode->pl().stops()[0].labelDeg = cand.deg;
      _stationLabels.back().s.labelDeg = cand.deg;
    }

    for (const auto &pref : nodeCandidates[idx].candidates[candIdx].neighborSides) {
      if (!pref.neighbor) continue;
      sidePrefs[pref.neighbor].push_back({node, pref.sign});
    }
  }

  for (size_t idx = 0; idx < _stationLabels.size(); ++idx) {
    auto &placed = _stationLabels[idx];
    const auto *n = _statLblNodes[idx];
    if (!n) continue;

    std::vector<size_t> crowding;
    _statLblIdx.get(placed.band, 0, &crowding);
    if (crowding.size() <= 1) {
      continue;
    }

    _statLblIdx.remove(idx);

    bool isTerminus = g.isTerminus(n);
    int prefDeg = 0;
    if (n->pl().stops().size()) {
      const auto &sp = n->pl().stops().front().pos;
      const auto *cp = n->pl().getGeom();
      double dx = sp.getX() - cp->getX();
      double dy = sp.getY() - cp->getY();
      if (std::abs(dx) > 1e-9 || std::abs(dy) > 1e-9) {
        double ang = std::atan2(dy, dx) * 180.0 / M_PI;
        prefDeg = static_cast<int>(std::round(ang / kStationAngleDeg));
        prefDeg = (prefDeg % static_cast<int>(kStationAngleSteps) +
                   static_cast<int>(kStationAngleSteps)) %
                  static_cast<int>(kStationAngleSteps);
      }
    }

    size_t flippedDeg =
        (placed.deg + kStationAngleSteps / 2) % kStationAngleSteps;
    auto flippedBand = util::geo::rotate(
        getStationLblBand(n, placed.fontSize, static_cast<uint8_t>(placed.pos), g),
        kStationAngleDeg * flippedDeg, *n->pl().getGeom());

    auto box = util::geo::getBoundingBox(flippedBand);
    double diag = util::geo::dist(box.getLowerLeft(), box.getUpperRight());
    double searchRad =
        g.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing) +
        std::max(_cfg->stationLabelSize, diag);

    auto overlaps = getOverlaps(flippedBand, n, g, searchRad);

    auto crowd = computeStationFarCrowd(flippedBand, n, searchRad, g);
    const auto &neighEdges = crowd.neighborEdges;
    const auto &neighNodes = crowd.neighborNodes;
    double farCrowdPen = crowd.farCrowdPen;
    double area =
        (box.getUpperRight().getX() - box.getLowerLeft().getX()) *
        (box.getUpperRight().getY() - box.getLowerLeft().getY());
    double neighborCount =
        static_cast<double>(neighEdges.size() + neighNodes.size());
    double clusterPen = area > 0.0 ? neighborCount / area : neighborCount;

    bool outside =
        box.getLowerLeft().getX() < mapBox.getLowerLeft().getX() ||
        box.getLowerLeft().getY() < mapBox.getLowerLeft().getY() ||
        box.getUpperRight().getX() > mapBox.getUpperRight().getX() ||
        box.getUpperRight().getY() > mapBox.getUpperRight().getY();

    size_t diff =
        (flippedDeg + kStationAngleSteps - static_cast<size_t>(prefDeg)) %
        kStationAngleSteps;
    if (diff > kStationAngleSteps / 2)
      diff = kStationAngleSteps - diff;
    double sidePen = static_cast<double>(diff) * _cfg->sidePenaltyWeight;
    double termPen =
        isTerminus && (flippedDeg % (kStationAngleSteps / 4) != 0)
            ? kTerminusAnglePen
            : 0;

    double candAng = flippedDeg * kStationAngleDeg * M_PI / 180.0;
    double candVecX = std::cos(candAng);
    double candVecY = std::sin(candAng);
    double sameSidePen = 0.0;
    for (auto e : n->getAdjList()) {
      auto neigh = e->getFrom() == n ? e->getTo() : e->getFrom();
      if (neigh->pl().stops().empty()) continue;
      size_t neighDeg = neigh->pl().stops()[0].labelDeg;
      if (neighDeg == std::numeric_limits<size_t>::max()) continue;
      double edgeVecX =
          neigh->pl().getGeom()->getX() - n->pl().getGeom()->getX();
      double edgeVecY =
          neigh->pl().getGeom()->getY() - n->pl().getGeom()->getY();
      double neighAng = neighDeg * kStationAngleDeg * M_PI / 180.0;
      double neighVecX = std::cos(neighAng);
      double neighVecY = std::sin(neighAng);
      double candSide = edgeVecX * candVecY - edgeVecY * candVecX;
      double neighSide = edgeVecX * neighVecY - edgeVecY * neighVecX;
      if (candSide * neighSide < 0)
        sameSidePen +=
            _cfg->sameSidePenalty * _cfg->crowdingSameSideScale;
    }
    auto prefIt = sidePrefs.find(n);
    if (prefIt != sidePrefs.end()) {
      for (auto &pref : prefIt->second) {
        const auto *prefNeigh = pref.first;
        int desired = pref.second;
        double edgeVecX =
            prefNeigh->pl().getGeom()->getX() - n->pl().getGeom()->getX();
        double edgeVecY =
            prefNeigh->pl().getGeom()->getY() - n->pl().getGeom()->getY();
        double candSide = edgeVecX * candVecY - edgeVecY * candVecX;
        if (candSide * desired < 0)
          sameSidePen +=
              _cfg->sameSidePenalty * _cfg->crowdingSameSideScale;
      }
    }

    StationLabel flipped(
        PolyLine<double>(flippedBand[0]), flippedBand, placed.fontSize,
        placed.bold, flippedDeg, placed.pos, overlaps,
        sidePen + termPen + sameSidePen, _cfg->stationLineOverlapPenalty,
        clusterPen, farCrowdPen, outside, _cfg->clusterPenScale,
        _cfg->outsidePenalty,
        &_cfg->orientationPenalties, placed.s);

    if (flipped.getPen() < placed.getPen()) {
      placed = flipped;
      auto *mutableNode2 = const_cast<shared::linegraph::LineNode *>(n);
      if (!mutableNode2->pl().stops().empty()) {
        mutableNode2->pl().stops()[0].labelDeg = flippedDeg;
        placed.s.labelDeg = flippedDeg;
      }
      for (auto e : n->getAdjList()) {
        auto neigh = e->getFrom() == n ? e->getTo() : e->getFrom();
        if (neigh->pl().stops().empty()) continue;
        double edgeVecX =
            neigh->pl().getGeom()->getX() - n->pl().getGeom()->getX();
        double edgeVecY =
            neigh->pl().getGeom()->getY() - n->pl().getGeom()->getY();
        double candSide = edgeVecX * candVecY - edgeVecY * candVecX;
        int sign = candSide >= 0 ? 1 : -1;
        auto prefIt2 = sidePrefs.find(neigh);
        if (prefIt2 == sidePrefs.end()) continue;
        for (auto &pref : prefIt2->second) {
          if (pref.first == n) {
            pref.second = sign;
            break;
          }
        }
      }
    }

    _statLblIdx.add(_stationLabels[idx].band, idx);
  }

  for (auto &entry : nodeCandidates) {
    entry.candidates.clear();
  }

  for (auto n : orderedNds) {
    auto it = labelIndex.find(n);
    if (it == labelIndex.end()) continue;
    size_t idx = it->second;
    auto &lbl = _stationLabels[idx];
    int same = 0;
    int opp = 0;
    for (auto e : n->getAdjList()) {
      auto neigh = e->getFrom() == n ? e->getTo() : e->getFrom();
      auto nit = labelIndex.find(neigh);
      if (nit == labelIndex.end()) continue;
      auto &nlbl = _stationLabels[nit->second];
      double edgeVecX =
          neigh->pl().getGeom()->getX() - n->pl().getGeom()->getX();
      double edgeVecY =
          neigh->pl().getGeom()->getY() - n->pl().getGeom()->getY();
      double lblAng = lbl.deg * kStationAngleDeg * M_PI / 180.0;
      double lblVecX = std::cos(lblAng);
      double lblVecY = std::sin(lblAng);
      double lblSide = edgeVecX * lblVecY - edgeVecY * lblVecX;
      double neighAng = nlbl.deg * kStationAngleDeg * M_PI / 180.0;
      double neighVecX = std::cos(neighAng);
      double neighVecY = std::sin(neighAng);
      double neighSide = edgeVecX * neighVecY - edgeVecY * neighVecX;
      if (lblSide * neighSide < 0)
        opp++;
      else
        same++;
    }
    if (opp > same) {
      _statLblIdx.remove(idx);
      lbl.deg = (lbl.deg + kStationAngleSteps / 2) % kStationAngleSteps;
      lbl.band = util::geo::rotate(
          getStationLblBand(n, lbl.fontSize, static_cast<uint8_t>(lbl.pos), g),
          kStationAngleDeg * lbl.deg, *n->pl().getGeom());
      lbl.geom = PolyLine<double>(lbl.band[0]);
      auto *nn = const_cast<shared::linegraph::LineNode *>(n);
      if (!nn->pl().stops().empty()) {
        nn->pl().stops()[0].labelDeg = lbl.deg;
        lbl.s.labelDeg = lbl.deg;
      }
      _statLblIdx.add(lbl.band, idx);
    }
  }

  for (int i = 0; i < _cfg->repositionLabel; ++i) {
    repositionStationLabels(g);
  }
}

// _____________________________________________________________________________
void Labeller::repositionStationLabels(const RenderGraph &g) {
  auto mapBox = g.getBBox();
  for (size_t idx = 0; idx < _stationLabels.size(); ++idx) {
    auto &placed = _stationLabels[idx];
    const auto *n = _statLblNodes[idx];
    if (!n) continue;

    _statLblIdx.remove(idx);

    int prefDegInt = 0;
    if (n->pl().stops().size()) {
      const auto &sp = n->pl().stops().front().pos;
      const auto *cp = n->pl().getGeom();
      double dx = sp.getX() - cp->getX();
      double dy = sp.getY() - cp->getY();
      if (std::abs(dx) > 1e-9 || std::abs(dy) > 1e-9) {
        double ang = std::atan2(dy, dx) * 180.0 / M_PI;
        prefDegInt = static_cast<int>(std::round(ang / kStationAngleDeg));
        int steps = static_cast<int>(kStationAngleSteps);
        prefDegInt %= steps;
        if (prefDegInt < 0) {
          prefDegInt += steps;
        }
      }
    }

    size_t prefDeg = static_cast<size_t>(prefDegInt);

    StationLabel best = placed;
    double bestScore = placed.getPen();
    bool isTerminus = g.isTerminus(n);

    for (uint8_t pos = 0; pos < 3; ++pos) {
      for (size_t flip = 0; flip < 2; ++flip) {
        size_t deg = placed.deg;
        if (flip == 1)
          deg = (placed.deg + kStationAngleSteps / 2) % kStationAngleSteps;

        auto band = util::geo::rotate(
            getStationLblBand(n, placed.fontSize, pos, g),
            kStationAngleDeg * deg, *n->pl().getGeom());

        auto box = util::geo::getBoundingBox(band);
        double diag = util::geo::dist(box.getLowerLeft(), box.getUpperRight());
        double searchRad = g.getMaxLineNum() *
                               (_cfg->lineWidth + _cfg->lineSpacing) +
                           std::max(_cfg->stationLabelSize, diag);

        auto overlaps = getOverlaps(band, n, g, searchRad);

        auto crowd = computeStationFarCrowd(band, n, searchRad, g);
        const auto &neighEdges = crowd.neighborEdges;
        const auto &neighNodes = crowd.neighborNodes;
        double farCrowdPen = crowd.farCrowdPen;
        double area = (box.getUpperRight().getX() - box.getLowerLeft().getX()) *
                      (box.getUpperRight().getY() - box.getLowerLeft().getY());
        double neighborCount =
            static_cast<double>(neighEdges.size() + neighNodes.size());
        double clusterPen =
            area > 0.0 ? neighborCount / area : neighborCount;

        bool outside = box.getLowerLeft().getX() < mapBox.getLowerLeft().getX() ||
                        box.getLowerLeft().getY() < mapBox.getLowerLeft().getY() ||
                        box.getUpperRight().getX() > mapBox.getUpperRight().getX() ||
                        box.getUpperRight().getY() > mapBox.getUpperRight().getY();

        size_t diff =
            (deg + kStationAngleSteps - prefDeg) % kStationAngleSteps;
        if (diff > kStationAngleSteps / 2)
          diff = kStationAngleSteps - diff;
        double sidePen =
            static_cast<double>(diff) * _cfg->sidePenaltyWeight;

        double candAng = deg * kStationAngleDeg * M_PI / 180.0;
        double candVecX = std::cos(candAng);
        double candVecY = std::sin(candAng);
        double sameSidePen = 0.0;
        for (auto e : n->getAdjList()) {
          auto neigh = e->getFrom() == n ? e->getTo() : e->getFrom();
          if (neigh->pl().stops().empty())
            continue;
          size_t neighDeg = neigh->pl().stops()[0].labelDeg;
          if (neighDeg == std::numeric_limits<size_t>::max())
            continue;
          double edgeVecX =
              neigh->pl().getGeom()->getX() - n->pl().getGeom()->getX();
          double edgeVecY =
              neigh->pl().getGeom()->getY() - n->pl().getGeom()->getY();
          double neighAng = neighDeg * kStationAngleDeg * M_PI / 180.0;
          double neighVecX = std::cos(neighAng);
          double neighVecY = std::sin(neighAng);
          double candSide = edgeVecX * candVecY - edgeVecY * candVecX;
          double neighSide = edgeVecX * neighVecY - edgeVecY * neighVecX;
          if (candSide * neighSide < 0)
            sameSidePen +=
                _cfg->sameSidePenalty * _cfg->crowdingSameSideScale;
        }

        double termPen =
            isTerminus && (deg % (kStationAngleSteps / 4) != 0)
                ? kTerminusAnglePen
                : 0;

        StationLabel cand(
            PolyLine<double>(band[0]), band, placed.fontSize, placed.bold, deg,
            pos, overlaps, sidePen + termPen + sameSidePen,
            _cfg->stationLineOverlapPenalty, clusterPen, farCrowdPen, outside,
            _cfg->clusterPenScale, _cfg->outsidePenalty,
            &_cfg->orientationPenalties, placed.s);

        double score = cand.getPen();
        if (score < bestScore) {
          bestScore = score;
          best = cand;
        }
      }
    }

    placed = best;
    auto *nn = const_cast<shared::linegraph::LineNode *>(n);
    if (!nn->pl().stops().empty()) {
      nn->pl().stops()[0].labelDeg = best.deg;
      placed.s.labelDeg = best.deg;
    }
    _statLblIdx.add(placed.band, idx);
  }
}

// _____________________________________________________________________________
Overlaps Labeller::getOverlaps(const util::geo::MultiLine<double> &band,
                               const shared::linegraph::LineNode *forNd,
                               const RenderGraph &g, double radius) const {
  std::set<const shared::linegraph::LineEdge *> proced;

  Overlaps ret{0, 0, 0, 0, 0};

  std::unordered_set<const shared::linegraph::Line *> overlappedLines;
  std::set<const shared::linegraph::LineNode *> procedNds{forNd};

  for (auto line : band) {
    auto neighs = g.getNeighborEdges(line, radius);
    for (auto neigh : neighs) {
      if (proced.count(neigh))
        continue;

      if (util::geo::dist(*neigh->pl().getGeom(), band) <
          g.getTotalWidth(neigh) / 2) {
        if (_cfg->stationLineOverlapPerLine) {
          for (const auto &lineOcc : neigh->pl().getLines()) {
            overlappedLines.insert(lineOcc.line);
          }
        } else {
          ret.lineOverlaps++;
        }
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

  if (_cfg->stationLineOverlapPerLine) {
    ret.lineOverlaps += overlappedLines.size();
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

          Overlaps overlaps{0, 0, 0, 0, 0};

          for (auto neigh : g.getNeighborEdges(
                   cand.getLine(),
                   g.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing) +
                       fontSize * 4)) {
            if (neigh == e)
              continue;
            if (util::geo::dist(cand.getLine(), *neigh->pl().getGeom()) <
                (g.getTotalWidth(neigh) / 2) + (fontSize)) {
              overlaps.lineOverlaps++;
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
              overlaps.statLabelOverlaps++;
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
              overlaps.lineLabelOverlaps++;
            }
          }

          cands.push_back({cand, fabs((geomLen / 2) - (start + (labelW / 2))),
                           fontSize, lines, overlaps});
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
  for (const auto &sl : _stationLabels)
    ret.push_back(sl.deg);
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

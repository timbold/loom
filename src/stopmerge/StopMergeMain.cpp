// Copyright 2016
// University of Freiburg - Chair of Algorithms and Data Structures
// Author: Patrick Brosi

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "3rdparty/json.hpp"
#include "shared/linegraph/LineGraph.h"
#include "stopmerge/config/ConfigReader.h"
#include "stopmerge/config/StopMergeConfig.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"
#include "util/geo/RTree.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using shared::linegraph::LineEdge;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::linegraph::Station;
using util::geo::DPoint;
using util::geo::Point;
using util::geo::PolyLine;
using util::geo::LinePoint;
using util::DEBUG;
using util::ERROR;
using util::WARN;

namespace {

struct UnionFind {
  explicit UnionFind(size_t n) : parent(n), size(n, 1) {
    for (size_t i = 0; i < n; ++i) parent[i] = i;
  }

  size_t find(size_t x) {
    if (parent[x] == x) return x;
    parent[x] = find(parent[x]);
    return parent[x];
  }

  void unite(size_t a, size_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return;
    if (size[a] < size[b]) std::swap(a, b);
    parent[b] = a;
    size[a] += size[b];
  }

  size_t compSize(size_t x) { return size[find(x)]; }

  std::vector<size_t> parent;
  std::vector<size_t> size;
};

struct EdgeInfo {
  LineEdge* edge = nullptr;
  size_t index = 0;
  double length = 0.0;
  std::set<std::string> lines;
  LineNode* from = nullptr;
  LineNode* to = nullptr;
};

struct StationInfo {
  LineNode* node = nullptr;
  std::string id;
  std::string name;
  DPoint pos;
  std::set<std::string> lines;

  LineEdge* edge = nullptr;
  size_t edgeIndex = std::numeric_limits<size_t>::max();
  double projDist = std::numeric_limits<double>::infinity();
  double chainageM = 0.0;
  DPoint projPoint;
  double bearingDeg = 0.0;

  LineNode* nearNode = nullptr;

  size_t corridorGroup = std::numeric_limits<size_t>::max();
  size_t superGroup = std::numeric_limits<size_t>::max();

  double axisPos = 0.0;
  bool assigned = false;
};

struct CorridorGroup {
  size_t id = 0;
  std::vector<size_t> edgeIndices;
  double length = 0.0;
  double axisBearingDeg = 0.0;
  std::set<std::string> lines;
  std::vector<DPoint> samples;
  util::geo::Box<double> bbox;
  DPoint centroid;
};

struct Overrides {
  std::vector<std::vector<std::string>> forceMerge;
  std::unordered_set<std::string> neverMerge;
  std::unordered_map<std::string, std::string> rename;
  bool empty() const {
    return forceMerge.empty() && neverMerge.empty() && rename.empty();
  }
};

struct MergeCluster {
  std::vector<size_t> members;  // indices into stations
  std::string decisionReason;
  std::string modeUsed;
  size_t groupId = std::numeric_limits<size_t>::max();
};

static double degToRad(double deg) { return deg * M_PI / 180.0; }
static double radToDeg(double rad) { return rad * 180.0 / M_PI; }

static double angleDiffDeg(double a, double b) {
  double diff = fabs(a - b);
  while (diff > 360.0) diff -= 360.0;
  if (diff > 180.0) diff = 360.0 - diff;
  return diff;
}

static double angleDiffUndirectedDeg(double a, double b) {
  double diff = angleDiffDeg(a, b);
  if (diff > 90.0) diff = 180.0 - diff;
  return diff;
}

static double bearingDeg(const DPoint& a, const DPoint& b) {
  return radToDeg(util::geo::angBetween(a, b));
}

static DPoint edgeDirPoint(LineEdge* e, LineNode* n) {
  const auto& pl = e->pl().getPolyline();
  if (pl.getLine().size() > 2) {
    if (e->getTo() == n) {
      return pl.getLine()[pl.getLine().size() - 2];
    }
    return pl.getLine()[1];
  }
  return *e->getOtherNd(n)->pl().getGeom();
}

static double edgeAngleAtNode(LineEdge* e, LineNode* n) {
  const auto* np = n->pl().getGeom();
  if (!np) return 0.0;
  DPoint other = edgeDirPoint(e, n);
  return bearingDeg(*np, other);
}

static std::set<std::string> edgeLineIds(const LineEdge* e) {
  std::set<std::string> ret;
  for (const auto& occ : e->pl().getLines()) {
    ret.insert(occ.line->id());
  }
  return ret;
}

static std::set<std::string> stationLineIds(const LineNode* n) {
  std::set<std::string> ret;
  for (const auto* l : LineGraph::servedLines(n)) ret.insert(l->id());
  return ret;
}

static bool setIntersects(const std::set<std::string>& a,
                          const std::set<std::string>& b) {
  if (a.empty() || b.empty()) return false;
  auto itA = a.begin();
  auto itB = b.begin();
  while (itA != a.end() && itB != b.end()) {
    if (*itA == *itB) return true;
    if (*itA < *itB)
      ++itA;
    else
      ++itB;
  }
  return false;
}

static size_t setIntersectionCount(const std::set<std::string>& a,
                                   const std::set<std::string>& b) {
  size_t count = 0;
  auto itA = a.begin();
  auto itB = b.begin();
  while (itA != a.end() && itB != b.end()) {
    if (*itA == *itB) {
      ++count;
      ++itA;
      ++itB;
      continue;
    }
    if (*itA < *itB)
      ++itA;
    else
      ++itB;
  }
  return count;
}

static std::string coordKey(const DPoint& p) {
  std::ostringstream out;
  out.setf(std::ios::fixed);
  out.precision(6);
  out << p.getX() << "|" << p.getY();
  return out.str();
}

static std::string stationStableId(const StationInfo& s) {
  if (!s.id.empty()) return s.id;
  return std::string("coord:") + coordKey(s.pos);
}

static bool parseBaseSuffix(const std::string& name, std::string* base,
                            std::string* suffix) {
  size_t firstSlash = name.find('/');
  if (firstSlash == std::string::npos) return false;
  size_t lastSlash = name.find_last_of('/');
  if (lastSlash <= firstSlash) return false;
  if (lastSlash != name.size() - 1) return false;
  std::string rawBase = util::trim(name.substr(0, firstSlash));
  std::string rawSuffix = util::trim(name.substr(firstSlash + 1,
                                                 lastSlash - firstSlash - 1));
  if (rawBase.empty() || rawSuffix.empty()) return false;
  *base = rawBase;
  *suffix = rawSuffix;
  return true;
}

static std::vector<std::string> splitSuffixTokens(const std::string& suffix) {
  std::vector<std::string> ret;
  std::stringstream ss(suffix);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token = util::trim(token);
    if (!token.empty()) ret.push_back(token);
  }
  return ret;
}

static double median(std::vector<double>* vals) {
  if (vals->empty()) return 0.0;
  std::sort(vals->begin(), vals->end());
  size_t mid = vals->size() / 2;
  if (vals->size() % 2 == 1) return (*vals)[mid];
  return 0.5 * ((*vals)[mid - 1] + (*vals)[mid]);
}

static std::vector<DPoint> samplePolyline(const PolyLine<double>& pl,
                                          double spacing, size_t cap) {
  std::vector<DPoint> pts;
  double len = pl.getLength();
  if (len <= 0.0) {
    pts.push_back(pl.front());
    return pts;
  }
  for (double d = 0.0; d <= len; d += spacing) {
    pts.push_back(pl.getPointAtDist(d).p);
  }
  if (util::geo::dist(pts.back(), pl.back()) > 0.1) pts.push_back(pl.back());

  if (pts.size() > cap) {
    size_t step = static_cast<size_t>(
        std::ceil(static_cast<double>(pts.size()) / static_cast<double>(cap)));
    std::vector<DPoint> reduced;
    for (size_t i = 0; i < pts.size(); i += step) reduced.push_back(pts[i]);
    if (util::geo::dist(reduced.back(), pts.back()) > 0.1) reduced.push_back(pts.back());
    pts.swap(reduced);
  }
  return pts;
}

static util::geo::Box<double> bboxFromPoints(const std::vector<DPoint>& pts) {
  if (pts.empty()) return util::geo::Box<double>();
  util::geo::Box<double> box(pts.front(), pts.front());
  for (const auto& p : pts) box = util::geo::extendBox(p, box);
  return box;
}

static double boxMinDist(const util::geo::Box<double>& a,
                         const util::geo::Box<double>& b) {
  double dx = 0.0;
  if (a.getUpperRight().getX() < b.getLowerLeft().getX())
    dx = b.getLowerLeft().getX() - a.getUpperRight().getX();
  else if (b.getUpperRight().getX() < a.getLowerLeft().getX())
    dx = a.getLowerLeft().getX() - b.getUpperRight().getX();

  double dy = 0.0;
  if (a.getUpperRight().getY() < b.getLowerLeft().getY())
    dy = b.getLowerLeft().getY() - a.getUpperRight().getY();
  else if (b.getUpperRight().getY() < a.getLowerLeft().getY())
    dy = a.getLowerLeft().getY() - b.getUpperRight().getY();

  return std::sqrt(dx * dx + dy * dy);
}

static util::json::Val toUtilJson(const nlohmann::json& j) {
  using util::json::Val;
  if (j.is_object()) {
    util::json::Dict d;
    for (auto it = j.begin(); it != j.end(); ++it) {
      d[it.key()] = toUtilJson(it.value());
    }
    return Val(d);
  }
  if (j.is_array()) {
    util::json::Array arr;
    for (const auto& v : j) arr.push_back(toUtilJson(v));
    return Val(arr);
  }
  if (j.is_string()) return Val(j.get<std::string>());
  if (j.is_boolean()) return Val(j.get<bool>());
  if (j.is_number_integer()) return Val(static_cast<int>(j.get<int64_t>()));
  if (j.is_number_unsigned()) return Val(static_cast<uint64_t>(j.get<uint64_t>()));
  if (j.is_number_float()) return Val(j.get<double>());
  return Val(util::json::Null());
}

static nlohmann::json configSnapshot(const stopmerge::config::StopMergeConfig& cfg) {
  nlohmann::json j;
  j["merge_stops"] = cfg.mergeStops;
  j["merge_stops_parallel_corridors"] = cfg.parallelCorridors;
  j["merge_stops_hub_fallback"] = cfg.mergeStopsHubFallback;
  j["merge_stops_snap_dist_m"] = cfg.mergeStopsSnapDistM;
  j["merge_stops_radius_m"] = cfg.mergeStopsRadiusM;
  j["merge_stops_chainage_m"] = cfg.mergeStopsChainageM;
  j["merge_stops_node_zone_m"] = cfg.mergeStopsNodeZoneM;
  j["merge_stops_corridor_angle_deg"] = cfg.mergeStopsCorridorAngleDeg;
  j["merge_stops_max_cluster_size"] = cfg.mergeStopsMaxClusterSize;
  j["parallel_pair_max_dist_m"] = cfg.parallelPairMaxDistM;
  j["parallel_pair_max_angle_deg"] = cfg.parallelPairMaxAngleDeg;
  j["parallel_pair_min_overlap_ratio"] = cfg.parallelPairMinOverlapRatio;
  j["parallel_pair_require_lines"] = cfg.parallelPairRequireLines;
  j["parallel_pair_min_group_length_m"] = cfg.parallelPairMinGroupLengthM;
  j["parallel_pair_max_median_dist_m"] = cfg.parallelPairMaxMedianDistM;
  j["hub_radius_m"] = cfg.hubRadiusM;
  j["merge_stops_max_local_density"] = cfg.mergeStopsMaxLocalDensity;
  j["hub_auto_require_base_name"] = cfg.hubAutoRequireBaseName;
  return j;
}

static Overrides loadOverrides(const std::string& path) {
  Overrides ov;
  if (path.empty()) return ov;
  std::ifstream in(path);
  if (!in) {
    LOG(ERROR) << "Could not open overrides file " << path;
    return ov;
  }
  nlohmann::json j;
  in >> j;
  if (j.contains("force_merge") && j["force_merge"].is_array()) {
    for (const auto& arr : j["force_merge"]) {
      if (!arr.is_array()) continue;
      std::vector<std::string> ids;
      for (const auto& id : arr) {
        if (id.is_string()) ids.push_back(id.get<std::string>());
      }
      if (!ids.empty()) ov.forceMerge.push_back(ids);
    }
  }
  if (j.contains("never_merge") && j["never_merge"].is_array()) {
    for (const auto& id : j["never_merge"]) {
      if (id.is_string()) ov.neverMerge.insert(id.get<std::string>());
    }
  }
  if (j.contains("rename") && j["rename"].is_object()) {
    for (auto it = j["rename"].begin(); it != j["rename"].end(); ++it) {
      if (it.value().is_string()) ov.rename[it.key()] = it.value().get<std::string>();
    }
  }
  return ov;
}

}  // namespace

int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  stopmerge::config::StopMergeConfig cfg;
  stopmerge::config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  if (cfg.mergeStops == "never") {
    std::cout << std::cin.rdbuf();
    return 0;
  }

  LineGraph lg;
  lg.readFromJson(&std::cin);

  // collect edges
  std::vector<EdgeInfo> edges;
  edges.reserve(lg.numEdgs());
  std::unordered_map<LineEdge*, size_t> edgeIndex;

  for (auto n : lg.getNds()) {
    for (auto e : n->getAdjListOut()) {
      if (e->getFrom() != n) continue;
      EdgeInfo info;
      info.edge = e;
      info.index = edges.size();
      info.length = e->pl().getPolyline().getLength();
      info.lines = edgeLineIds(e);
      info.from = e->getFrom();
      info.to = e->getTo();
      edgeIndex[e] = info.index;
      edges.push_back(info);
    }
  }

  // collect stations
  std::vector<StationInfo> stations;
  stations.reserve(lg.numNds());
  std::unordered_map<std::string, std::string> memberIdToLabel;
  std::unordered_map<LineNode*, std::string> nodeToStableId;
  for (auto n : lg.getNds()) {
    if (n->pl().stops().empty()) continue;
    const auto& st = n->pl().stops().front();
    StationInfo si;
    si.node = n;
    si.id = st.id;
    si.name = st.name;
    si.pos = *n->pl().getGeom();
    si.lines = stationLineIds(n);
    std::string stableId = stationStableId(si);
    memberIdToLabel[stableId] = si.name;
    nodeToStableId[n] = stableId;
    stations.push_back(si);
  }

  if (stations.empty()) {
    util::geo::output::GeoGraphJsonOutput out;
    util::geo::output::GeoJsonOutput jsonOut(std::cout);
    out.printLatLng(lg, &jsonOut);
    jsonOut.flush();
    return 0;
  }

  // assign stations to closest edge
  const auto* edgeGrid = lg.getEdgGrid();
  for (auto& st : stations) {
    std::vector<LineEdge*> candidates;
    edgeGrid->get(*st.node->pl().getGeom(), cfg.mergeStopsSnapDistM, &candidates);
    double bestScore = std::numeric_limits<double>::infinity();
    LineEdge* bestEdge = nullptr;
    LinePoint<double> bestLp;

    for (auto* e : candidates) {
      const auto& pl = e->pl().getPolyline();
      LinePoint<double> lp = pl.projectOn(*st.node->pl().getGeom());
      double dist = util::geo::dist(lp.p, *st.node->pl().getGeom());
      double score = dist;
      if (!st.lines.empty()) {
        auto eLines = edgeLineIds(e);
        if (setIntersects(st.lines, eLines)) score -= 1e-6;
      }

      if (score < bestScore - 1e-9 ||
          (fabs(score - bestScore) <= 1e-9 &&
           (!bestEdge || util::toString(e) < util::toString(bestEdge)))) {
        bestScore = score;
        bestEdge = e;
        bestLp = lp;
      }
    }

    if (bestEdge && bestScore <= cfg.mergeStopsSnapDistM) {
      st.edge = bestEdge;
      st.edgeIndex = edgeIndex[bestEdge];
      st.projDist = bestScore;
      st.projPoint = bestLp.p;
      const auto& pl = bestEdge->pl().getPolyline();
      double len = pl.getLength();
      st.chainageM = bestLp.totalPos * len;

      size_t segIdx = bestLp.lastIndex;
      const auto& line = pl.getLine();
      size_t nextIdx = std::min(segIdx + 1, line.size() - 1);
      st.bearingDeg = bearingDeg(line[segIdx], line[nextIdx]);
      st.assigned = true;

      double distFrom = util::geo::dist(st.projPoint,
                                       *bestEdge->getFrom()->pl().getGeom());
      double distTo = util::geo::dist(st.projPoint,
                                     *bestEdge->getTo()->pl().getGeom());
      if (distFrom <= cfg.mergeStopsNodeZoneM &&
          distFrom <= distTo) {
        st.nearNode = bestEdge->getFrom();
      } else if (distTo <= cfg.mergeStopsNodeZoneM) {
        st.nearNode = bestEdge->getTo();
      }
    }
  }

  // corridor grouping
  UnionFind edgeUf(edges.size());
  for (auto n : lg.getNds()) {
    const auto& adj = n->getAdjList();
    for (size_t i = 0; i < adj.size(); ++i) {
      for (size_t j = i + 1; j < adj.size(); ++j) {
        LineEdge* a = adj[i];
        LineEdge* b = adj[j];
        if (!edgeIndex.count(a) || !edgeIndex.count(b)) continue;
        double angA = edgeAngleAtNode(a, n);
        double angB = edgeAngleAtNode(b, n);
        double diff = angleDiffUndirectedDeg(angA, angB);
        if (diff > cfg.mergeStopsCorridorAngleDeg) continue;
        auto aLines = edgeLineIds(a);
        auto bLines = edgeLineIds(b);
        if (!aLines.empty() && !bLines.empty() && !setIntersects(aLines, bLines))
          continue;
        edgeUf.unite(edgeIndex[a], edgeIndex[b]);
      }
    }
  }

  // build corridor groups
  std::unordered_map<size_t, size_t> rootToGroup;
  std::vector<CorridorGroup> corridorGroups;
  for (const auto& e : edges) {
    size_t root = edgeUf.find(e.index);
    if (!rootToGroup.count(root)) {
      size_t id = corridorGroups.size();
      rootToGroup[root] = id;
      corridorGroups.emplace_back();
      corridorGroups.back().id = id;
    }
    CorridorGroup& g = corridorGroups[rootToGroup[root]];
    g.edgeIndices.push_back(e.index);
    g.length += e.length;
    g.lines.insert(e.lines.begin(), e.lines.end());
  }

  // compute corridor axes and samples
  for (auto& g : corridorGroups) {
    std::vector<double> bearings;
    bearings.reserve(g.edgeIndices.size() * 4);
    std::vector<DPoint> samples;

    for (size_t idx : g.edgeIndices) {
      const auto& pl = edges[idx].edge->pl().getPolyline();
      const auto& line = pl.getLine();
      for (size_t i = 1; i < line.size(); ++i) {
        double ang = bearingDeg(line[i - 1], line[i]);
        if (ang < 0) ang += 360.0;
        if (ang >= 180.0) ang -= 180.0;
        bearings.push_back(ang);
      }
      auto pts = samplePolyline(pl, 10.0, 200);
      samples.insert(samples.end(), pts.begin(), pts.end());
    }

    if (!bearings.empty()) {
      std::sort(bearings.begin(), bearings.end());
      g.axisBearingDeg = bearings[bearings.size() / 2];
    }

    if (samples.empty()) {
      g.centroid = DPoint(0, 0);
    } else {
      double sx = 0, sy = 0;
      for (const auto& p : samples) {
        sx += p.getX();
        sy += p.getY();
      }
      g.centroid = DPoint(sx / samples.size(), sy / samples.size());
    }

    if (samples.size() > 500) {
      size_t step = static_cast<size_t>(
          std::ceil(static_cast<double>(samples.size()) / 500.0));
      std::vector<DPoint> reduced;
      for (size_t i = 0; i < samples.size(); i += step) reduced.push_back(samples[i]);
      samples.swap(reduced);
    }

    g.samples = samples;
    g.bbox = bboxFromPoints(samples);
  }

  // assign corridor group to stations
  for (auto& st : stations) {
    if (!st.assigned) continue;
    size_t root = edgeUf.find(st.edgeIndex);
    if (!rootToGroup.count(root)) continue;
    st.corridorGroup = rootToGroup[root];
  }

  // parallel corridor pairing
  UnionFind corridorUf(corridorGroups.size());
  std::ofstream pairDebug;
  if (!cfg.parallelPairDebugCsv.empty()) {
    pairDebug.open(cfg.parallelPairDebugCsv);
    pairDebug << "groupA_id,groupB_id,lenA_m,lenB_m,min_dist_m,median_dist_m,"
                 "angle_diff_deg,overlap_ratio,shared_lines_count,decision,reason\n";
  }

  if (cfg.parallelCorridors != "off") {
    for (size_t i = 0; i < corridorGroups.size(); ++i) {
      for (size_t j = i + 1; j < corridorGroups.size(); ++j) {
        const auto& a = corridorGroups[i];
        const auto& b = corridorGroups[j];
        std::string decision = "skipped";
        std::string reason = "";

        if (a.length < cfg.parallelPairMinGroupLengthM ||
            b.length < cfg.parallelPairMinGroupLengthM) {
          reason = "length";
        } else if (angleDiffUndirectedDeg(a.axisBearingDeg, b.axisBearingDeg) >
                   cfg.parallelPairMaxAngleDeg) {
          reason = "angle";
        } else if (boxMinDist(a.bbox, b.bbox) > cfg.parallelPairMaxDistM) {
          reason = "bbox_dist";
        } else {
          const auto& small = (a.samples.size() <= b.samples.size()) ? a : b;
          const auto& large = (a.samples.size() <= b.samples.size()) ? b : a;
          std::vector<double> dists;
          dists.reserve(small.samples.size());
          double minDist = std::numeric_limits<double>::infinity();
          for (const auto& p : small.samples) {
            double best = std::numeric_limits<double>::infinity();
            for (const auto& q : large.samples) {
              double d = util::geo::dist(p, q);
              if (d < best) best = d;
            }
            dists.push_back(best);
            if (best < minDist) minDist = best;
          }

          double medDist = median(&dists);
          if (minDist > cfg.parallelPairMaxDistM) {
            reason = "min_dist";
          } else if (medDist > cfg.parallelPairMaxMedianDistM) {
            reason = "median_dist";
          } else {
            double axisDeg = (a.length >= b.length) ? a.axisBearingDeg
                                                    : b.axisBearingDeg;
            double axisRad = degToRad(axisDeg);
            DPoint axis(std::cos(axisRad), std::sin(axisRad));

            auto projInterval = [&](const std::vector<DPoint>& pts,
                                    double* minp, double* maxp) {
              *minp = std::numeric_limits<double>::infinity();
              *maxp = -std::numeric_limits<double>::infinity();
              for (const auto& p : pts) {
                double proj = p.getX() * axis.getX() + p.getY() * axis.getY();
                if (proj < *minp) *minp = proj;
                if (proj > *maxp) *maxp = proj;
              }
            };

            double minA, maxA, minB, maxB;
            projInterval(a.samples, &minA, &maxA);
            projInterval(b.samples, &minB, &maxB);
            double overlap = std::max(0.0, std::min(maxA, maxB) - std::max(minA, minB));
            double lenA = maxA - minA;
            double lenB = maxB - minB;
            double overlapRatio = (std::min(lenA, lenB) > 0.0)
                                      ? overlap / std::min(lenA, lenB)
                                      : 0.0;

            size_t sharedLines = setIntersectionCount(a.lines, b.lines);
            bool requireLines = cfg.parallelPairRequireLines ||
                                cfg.parallelCorridors == "auto";

            if (overlapRatio < cfg.parallelPairMinOverlapRatio) {
              reason = "overlap";
            } else if (requireLines && sharedLines == 0) {
              reason = "lines";
            } else {
              corridorUf.unite(i, j);
              decision = "paired";
              reason = "ok";
            }

            if (pairDebug.is_open()) {
              pairDebug << i << "," << j << "," << a.length << "," << b.length
                        << "," << minDist << "," << medDist << ","
                        << angleDiffUndirectedDeg(a.axisBearingDeg,
                                                  b.axisBearingDeg)
                        << "," << overlapRatio << "," << sharedLines << ","
                        << decision << "," << reason << "\n";
            }
            continue;
          }
        }

        if (pairDebug.is_open()) {
          size_t sharedLines = setIntersectionCount(a.lines, b.lines);
          pairDebug << i << "," << j << "," << a.length << "," << b.length
                    << "," << 0 << "," << 0 << ","
                    << angleDiffUndirectedDeg(a.axisBearingDeg, b.axisBearingDeg)
                    << "," << 0 << "," << sharedLines << "," << decision
                    << "," << reason << "\n";
        }
      }
    }
  }

  if (pairDebug.is_open()) pairDebug.close();

  // build super corridor groups
  std::unordered_map<size_t, size_t> superRootToGroup;
  std::vector<CorridorGroup> superGroups;
  for (size_t i = 0; i < corridorGroups.size(); ++i) {
    size_t root = corridorUf.find(i);
    if (!superRootToGroup.count(root)) {
      size_t id = superGroups.size();
      superRootToGroup[root] = id;
      superGroups.emplace_back();
      superGroups.back().id = id;
    }
    auto& g = superGroups[superRootToGroup[root]];
    g.length += corridorGroups[i].length;
    g.lines.insert(corridorGroups[i].lines.begin(),
                   corridorGroups[i].lines.end());
    g.edgeIndices.insert(g.edgeIndices.end(),
                         corridorGroups[i].edgeIndices.begin(),
                         corridorGroups[i].edgeIndices.end());
    g.samples.insert(g.samples.end(),
                     corridorGroups[i].samples.begin(),
                     corridorGroups[i].samples.end());
  }

  for (auto& g : superGroups) {
    if (g.samples.empty()) continue;
    if (g.samples.size() > 500) {
      size_t step = static_cast<size_t>(
          std::ceil(static_cast<double>(g.samples.size()) / 500.0));
      std::vector<DPoint> reduced;
      for (size_t i = 0; i < g.samples.size(); i += step) reduced.push_back(g.samples[i]);
      g.samples.swap(reduced);
    }
    double sx = 0, sy = 0;
    for (const auto& p : g.samples) {
      sx += p.getX();
      sy += p.getY();
    }
    g.centroid = DPoint(sx / g.samples.size(), sy / g.samples.size());
    g.bbox = bboxFromPoints(g.samples);

    std::vector<double> bearings;
    for (const auto& p : g.samples) {
      double ang = bearingDeg(g.centroid, p);
      if (ang < 0) ang += 360.0;
      if (ang >= 180.0) ang -= 180.0;
      bearings.push_back(ang);
    }
    if (!bearings.empty()) {
      std::sort(bearings.begin(), bearings.end());
      g.axisBearingDeg = bearings[bearings.size() / 2];
    }
  }

  // assign super groups to stations
  for (auto& st : stations) {
    if (!st.assigned) continue;
    size_t root = corridorUf.find(st.corridorGroup);
    if (!superRootToGroup.count(root)) continue;
    st.superGroup = superRootToGroup[root];
  }

  // compute axis positions per group
  auto computeAxisPos = [](const CorridorGroup& g, StationInfo* st) {
    double axisRad = degToRad(g.axisBearingDeg);
    DPoint axis(std::cos(axisRad), std::sin(axisRad));
    DPoint rel = DPoint(st->projPoint.getX() - g.centroid.getX(),
                        st->projPoint.getY() - g.centroid.getY());
    st->axisPos = rel.getX() * axis.getX() + rel.getY() * axis.getY();
  };

  if (cfg.parallelCorridors == "off") {
    for (auto& st : stations) {
      if (st.assigned && st.corridorGroup < corridorGroups.size())
        computeAxisPos(corridorGroups[st.corridorGroup], &st);
    }
  } else {
    for (auto& st : stations) {
      if (st.assigned && st.superGroup < superGroups.size())
        computeAxisPos(superGroups[st.superGroup], &st);
    }
  }

  // load overrides
  Overrides overrides = loadOverrides(cfg.mergeStopsOverrides);

  std::unordered_map<std::string, size_t> idToStation;
  for (size_t i = 0; i < stations.size(); ++i) {
    if (stations[i].id.empty()) continue;
    if (idToStation.count(stations[i].id)) {
      LOG(WARN) << "Duplicate station_id " << stations[i].id
                << " in overrides map; keeping first.";
      continue;
    }
    idToStation[stations[i].id] = i;
  }

  std::vector<bool> locked(stations.size(), false);
  for (const auto& id : overrides.neverMerge) {
    auto it = idToStation.find(id);
    if (it != idToStation.end()) locked[it->second] = true;
  }

  // force merge clusters
  std::vector<MergeCluster> finalClusters;
  std::vector<bool> assignedStation(stations.size(), false);

  for (const auto& group : overrides.forceMerge) {
    MergeCluster cluster;
    cluster.modeUsed = "override";
    cluster.decisionReason = "force_merge";
    for (const auto& id : group) {
      auto it = idToStation.find(id);
      if (it == idToStation.end()) continue;
      size_t idx = it->second;
      if (!stations[idx].assigned) {
        LOG(WARN) << "Station " << id
                  << " is not assigned to a corridor; skipping forced merge.";
        continue;
      }
      if (locked[idx]) {
        LOG(WARN) << "Station " << id
                  << " is in never_merge; skipping forced merge.";
        continue;
      }
      cluster.members.push_back(idx);
    }
    std::sort(cluster.members.begin(), cluster.members.end());
    cluster.members.erase(std::unique(cluster.members.begin(),
                                      cluster.members.end()),
                          cluster.members.end());
    if (cluster.members.size() >= 2) {
      for (auto idx : cluster.members) assignedStation[idx] = true;
      finalClusters.push_back(cluster);
    }
  }

  // build mergeable clusters by group
  std::unordered_map<size_t, std::vector<size_t>> groupToStations;
  for (size_t i = 0; i < stations.size(); ++i) {
    if (!stations[i].assigned) continue;
    if (locked[i]) continue;
    if (assignedStation[i]) continue;
    size_t gid = (cfg.parallelCorridors == "off") ? stations[i].corridorGroup
                                                   : stations[i].superGroup;
    if (gid == std::numeric_limits<size_t>::max()) continue;
    groupToStations[gid].push_back(i);
  }

  auto stationGeoDist = [&](size_t a, size_t b) {
    return util::geo::dist(stations[a].pos, stations[b].pos);
  };

  auto hasSharedBase = [&](const std::vector<size_t>& members, std::string* base) {
    std::unordered_map<std::string, size_t> baseCount;
    for (auto idx : members) {
      std::string b, s;
      if (!parseBaseSuffix(stations[idx].name, &b, &s)) continue;
      baseCount[b]++;
    }
    size_t bestCount = 0;
    std::string best;
    for (const auto& kv : baseCount) {
      if (kv.second > bestCount) {
        bestCount = kv.second;
        best = kv.first;
      }
    }
    if (bestCount > members.size() / 2) {
      if (base) *base = best;
      return true;
    }
    return false;
  };

  auto sharedLines = [&](const std::vector<size_t>& members) {
    if (members.empty()) return false;
    std::set<std::string> inter = stations[members.front()].lines;
    for (size_t i = 1; i < members.size(); ++i) {
      std::set<std::string> next;
      for (const auto& l : inter) {
        if (stations[members[i]].lines.count(l)) next.insert(l);
      }
      inter.swap(next);
      if (inter.empty()) return false;
    }
    return !inter.empty();
  };

  auto tightNodeZone = [&](const std::vector<size_t>& members) {
    if (members.size() > 4) return false;
    LineNode* base = nullptr;
    for (auto idx : members) {
      if (!stations[idx].nearNode) return false;
      if (!base) base = stations[idx].nearNode;
      if (base != stations[idx].nearNode) return false;
    }
    return base != nullptr;
  };

  for (const auto& kv : groupToStations) {
    const auto& stIdxs = kv.second;
    if (stIdxs.size() < 2) continue;

    UnionFind stUf(stIdxs.size());
    for (size_t i = 0; i < stIdxs.size(); ++i) {
      for (size_t j = i + 1; j < stIdxs.size(); ++j) {
        size_t a = stIdxs[i];
        size_t b = stIdxs[j];
        if (stationGeoDist(a, b) > cfg.mergeStopsRadiusM) continue;
        bool nearSameNode = (stations[a].nearNode &&
                             stations[a].nearNode == stations[b].nearNode);
        bool chainageClose =
            fabs(stations[a].axisPos - stations[b].axisPos) <=
            cfg.mergeStopsChainageM;
        if (nearSameNode || chainageClose) stUf.unite(i, j);
      }
    }

    std::unordered_map<size_t, std::vector<size_t>> comps;
    for (size_t i = 0; i < stIdxs.size(); ++i) {
      comps[stUf.find(i)].push_back(stIdxs[i]);
    }

    for (auto& comp : comps) {
      auto members = comp.second;
      if (members.size() < 2) continue;

      MergeCluster cluster;
      cluster.members = members;
      cluster.groupId = kv.first;
      cluster.modeUsed = cfg.mergeStops;

      bool sizeOk = (cfg.mergeStops != "auto") ||
                    (members.size() <= cfg.mergeStopsMaxClusterSize);
      if (!sizeOk) {
        continue;
      }

      std::string baseName;
      bool baseShared = hasSharedBase(members, &baseName);
      bool linesShared = sharedLines(members);
      bool tightZone = tightNodeZone(members);

      if (cfg.mergeStops == "auto") {
        if (!(baseShared || linesShared || tightZone)) {
          continue;
        }
      }

      // hub protection
      double sx = 0, sy = 0;
      for (auto idx : members) {
        sx += stations[idx].pos.getX();
        sy += stations[idx].pos.getY();
      }
      DPoint centroid(sx / members.size(), sy / members.size());
      size_t localDensity = 0;
      for (const auto& st : stations) {
        if (util::geo::dist(st.pos, centroid) <= cfg.hubRadiusM) localDensity++;
      }

      if (localDensity > cfg.mergeStopsMaxLocalDensity) {
        bool allow = false;
        if (cfg.hubAutoRequireBaseName) {
          allow = baseShared || linesShared;
        } else {
          allow = linesShared || tightZone;
        }
        if (!allow) {
          continue;
        }
      }

      cluster.decisionReason = "auto_rules";
      finalClusters.push_back(cluster);
      for (auto idx : members) assignedStation[idx] = true;
    }
  }

  // hub fallback (cross-corridor) merging
  std::ofstream rejectDebug;
  if (!cfg.mergeStopsDebugRejectsCsv.empty()) {
    rejectDebug.open(cfg.mergeStopsDebugRejectsCsv);
    rejectDebug << "stopA,stopB,dist_m,same_corridor_group,"
                   "same_supercorridor_group,node_zone_ok,chainage_ok,"
                   "density,decision,reject_reason\n";
  }

  bool hubFallbackEnabled =
      (cfg.mergeStopsHubFallback != "off") && cfg.mergeStops != "never";
  if (hubFallbackEnabled) {
    std::vector<size_t> candidates;
    candidates.reserve(stations.size());
    for (size_t i = 0; i < stations.size(); ++i) {
      if (!stations[i].assigned) continue;
      if (locked[i]) continue;
      if (assignedStation[i]) continue;
      candidates.push_back(i);
    }

    util::geo::RTree<size_t, util::geo::Point, double> stopGrid;
    for (auto idx : candidates) {
      stopGrid.add(stations[idx].pos, idx);
    }

    std::unordered_map<size_t, std::string> baseNameFor;
    for (auto idx : candidates) {
      std::string base, suffix;
      if (parseBaseSuffix(stations[idx].name, &base, &suffix)) {
        baseNameFor[idx] = base;
      }
    }

    UnionFind hubUf(candidates.size());
    std::unordered_map<size_t, size_t> idxToLocal;
    for (size_t i = 0; i < candidates.size(); ++i) {
      idxToLocal[candidates[i]] = i;
    }

    for (size_t i = 0; i < candidates.size(); ++i) {
      size_t a = candidates[i];
      if (!baseNameFor.count(a)) continue;
      std::vector<size_t> neighs;
      stopGrid.get(stations[a].pos, cfg.mergeStopsRadiusM, &neighs);
      for (auto b : neighs) {
        if (b <= a) continue;
        if (!baseNameFor.count(b)) continue;
        if (baseNameFor[a] != baseNameFor[b]) continue;

        double distAB = util::geo::dist(stations[a].pos, stations[b].pos);
        if (distAB > cfg.mergeStopsRadiusM) continue;

        bool sameCorridor = stations[a].corridorGroup == stations[b].corridorGroup;
        bool sameSuper = stations[a].superGroup == stations[b].superGroup;
        bool nodeZoneOk = (stations[a].nearNode &&
                           stations[a].nearNode == stations[b].nearNode);
        bool chainageOk =
            fabs(stations[a].axisPos - stations[b].axisPos) <=
            cfg.mergeStopsChainageM;

        DPoint centroid((stations[a].pos.getX() + stations[b].pos.getX()) / 2.0,
                        (stations[a].pos.getY() + stations[b].pos.getY()) / 2.0);
        size_t density = 0;
        for (const auto& st : stations) {
          if (util::geo::dist(st.pos, centroid) <= cfg.hubRadiusM) density++;
        }

        std::string decision = "skipped";
        std::string reason = "";

        if (density > cfg.mergeStopsMaxLocalDensity) {
          reason = "density";
        } else {
          size_t la = idxToLocal[a];
          size_t lb = idxToLocal[b];
          size_t sizeA = hubUf.compSize(la);
          size_t sizeB = hubUf.compSize(lb);
          if (sizeA + sizeB > cfg.mergeStopsMaxClusterSize) {
            reason = "max_cluster";
          } else {
            hubUf.unite(la, lb);
            decision = "merged";
            reason = "hub_fallback";
          }
        }

        if (rejectDebug.is_open()) {
          rejectDebug << stationStableId(stations[a]) << ","
                      << stationStableId(stations[b]) << ","
                      << distAB << "," << (sameCorridor ? "true" : "false")
                      << "," << (sameSuper ? "true" : "false") << ","
                      << (nodeZoneOk ? "true" : "false") << ","
                      << (chainageOk ? "true" : "false") << ","
                      << density << "," << decision << "," << reason << "\n";
        }
      }
    }

    std::unordered_map<size_t, std::vector<size_t>> comps;
    for (size_t i = 0; i < candidates.size(); ++i) {
      comps[hubUf.find(i)].push_back(candidates[i]);
    }

    for (auto& comp : comps) {
      auto members = comp.second;
      if (members.size() < 2) continue;
      MergeCluster cluster;
      cluster.members = members;
      cluster.modeUsed = "hub_fallback";
      cluster.decisionReason = "hub_fallback";
      finalClusters.push_back(cluster);
      for (auto idx : members) assignedStation[idx] = true;
    }
  }

  if (rejectDebug.is_open()) rejectDebug.close();

  // build mapping and merge
  std::unordered_map<std::string, std::string> memberToMerged;
  std::unordered_map<std::string, std::vector<std::string>> mergedMembers;
  std::unordered_map<std::string, std::vector<std::string>> mergedMemberLabels;
  std::unordered_map<std::string, std::string> mergedReasons;
  nlohmann::json mergeMapJson;
  mergeMapJson["config"] = configSnapshot(cfg);
  mergeMapJson["merged_stops"] = nlohmann::json::object();
  mergeMapJson["member_to_merged"] = nlohmann::json::object();

  std::ofstream mergeDebug;
  if (!cfg.mergeStopsDebugCsv.empty()) {
    mergeDebug.open(cfg.mergeStopsDebugCsv);
    mergeDebug << "merged_id,merged_name,members_ids,members_names,"
                  "centroid_lat,centroid_lon,group_id,cluster_size,mode_used,"
                  "decision_reason\n";
  }

  for (auto& cluster : finalClusters) {
    std::vector<size_t> members = cluster.members;
    std::sort(members.begin(), members.end());

    // choose merged id
    std::string mergedId;
    std::string bestKey;
    size_t baseIdx = members.front();
    for (auto idx : members) {
      std::string key = stationStableId(stations[idx]);
      if (mergedId.empty() || key < bestKey) {
        bestKey = key;
        mergedId = key;
        baseIdx = idx;
      }
    }

    // build merged name
    std::vector<std::string> memberNames;
    memberNames.reserve(members.size());
    for (auto idx : members) memberNames.push_back(stations[idx].name);

    std::string baseName;
    bool baseShared = hasSharedBase(members, &baseName);
    std::string mergedName;
    if (baseShared) {
      std::vector<std::string> suffixes;
      std::unordered_set<std::string> seen;
      for (auto idx : members) {
        std::string b, s;
        if (!parseBaseSuffix(stations[idx].name, &b, &s)) continue;
        if (b != baseName) continue;
        for (const auto& tok : splitSuffixTokens(s)) {
          if (seen.insert(tok).second) suffixes.push_back(tok);
        }
      }
      std::ostringstream out;
      out << baseName << " /";
      for (size_t i = 0; i < suffixes.size(); ++i) {
        if (i) out << ", ";
        out << suffixes[i];
      }
      out << "/";
      mergedName = out.str();
    } else {
      std::unordered_set<std::string> seen;
      std::ostringstream out;
      bool first = true;
      for (auto idx : members) {
        const auto& nm = stations[idx].name;
        if (seen.insert(nm).second) {
          if (!first) out << ", ";
          out << nm;
          first = false;
        }
      }
      mergedName = out.str();
    }

    if (overrides.rename.count(mergedId)) mergedName = overrides.rename[mergedId];

    double sx = 0, sy = 0;
    for (auto idx : members) {
      sx += stations[idx].pos.getX();
      sy += stations[idx].pos.getY();
    }
    DPoint centroid(sx / members.size(), sy / members.size());

    // merge nodes into base
    LineNode* base = stations[baseIdx].node;
    for (auto idx : members) {
      if (idx == baseIdx) continue;
      LineNode* other = stations[idx].node;
      base = lg.mergeNds(other, base);
    }

    base->pl().setGeom(centroid);
    base->pl().clearStops();
    base->pl().addStop(Station(mergedId, mergedName, centroid));

    for (auto idx : members) {
      memberToMerged[stationStableId(stations[idx])] = mergedId;
      mergeMapJson["member_to_merged"][stationStableId(stations[idx])] = mergedId;
    }

    nlohmann::json clusterInfo;
    std::vector<std::string> memberIds;
    for (auto idx : members) memberIds.push_back(stationStableId(stations[idx]));
    std::sort(memberIds.begin(), memberIds.end());
    memberIds.erase(std::unique(memberIds.begin(), memberIds.end()),
                    memberIds.end());
    clusterInfo["members"] = memberIds;
    clusterInfo["merged_name"] = mergedName;
    auto ll = util::geo::webMercToLatLng<double>(centroid.getX(), centroid.getY());
    clusterInfo["centroid"] = {ll.getX(), ll.getY()};
    clusterInfo["group_id"] = cluster.groupId;
    clusterInfo["reasons"] = cluster.decisionReason;
    mergeMapJson["merged_stops"][mergedId] = clusterInfo;

    mergedMembers[mergedId] = memberIds;
    mergedReasons[mergedId] = cluster.decisionReason;
    std::vector<std::string> labels;
    labels.reserve(memberIds.size());
    for (const auto& mid : memberIds) {
      auto it = memberIdToLabel.find(mid);
      if (it != memberIdToLabel.end()) labels.push_back(it->second);
    }
    mergedMemberLabels[mergedId] = labels;

    if (mergeDebug.is_open()) {
      std::ostringstream ids;
      std::ostringstream names;
      for (size_t i = 0; i < members.size(); ++i) {
        if (i) {
          ids << "|";
          names << "|";
        }
        ids << stationStableId(stations[members[i]]);
        names << stations[members[i]].name;
      }
      mergeDebug << mergedId << "," << mergedName << "," << ids.str() << ","
                 << names.str() << "," << ll.getY() << "," << ll.getX() << ","
                 << cluster.groupId << "," << members.size() << ","
                 << cluster.modeUsed << "," << cluster.decisionReason << "\n";
    }
  }

  if (mergeDebug.is_open()) mergeDebug.close();

  if (!cfg.stopMergeMapJson.empty()) {
    std::ofstream out(cfg.stopMergeMapJson);
    out << mergeMapJson.dump(2);
  }

  // update graph properties with mapping
  nlohmann::json graphProps = lg.getGraphProps();
  graphProps["stopmerge"]["member_to_merged"] = mergeMapJson["member_to_merged"];
  graphProps["stopmerge"]["config"] = mergeMapJson["config"];

  // write stopmerge metadata into station nodes
  for (auto n : lg.getNds()) {
    if (n->pl().stops().empty()) continue;
    const auto& st = n->pl().stops().front();
    std::string stationId = st.id;
    if (stationId.empty()) {
      auto it = nodeToStableId.find(n);
      if (it != nodeToStableId.end()) stationId = it->second;
    }
    std::string repId = stationId;
    auto it = memberToMerged.find(stationId);
    if (it != memberToMerged.end()) repId = it->second;

    auto mit = mergedMembers.find(repId);
    if (mit != mergedMembers.end() && mit->second.size() >= 2 &&
        stationId == repId) {
      const auto& members = mit->second;
      const auto& labels = mergedMemberLabels[repId];
      const auto& reason = mergedReasons[repId];
      n->pl().setStopmergeMeta(true, members.size(), members, labels, reason);
    } else {
      n->pl().setStopmergeMeta(false, 1, {}, {}, "");
    }
  }

  util::geo::output::GeoGraphJsonOutput out;
  util::geo::output::GeoJsonOutput jsonOut(std::cout, toUtilJson(graphProps));
  out.printLatLng(lg, &jsonOut);
  jsonOut.flush();

  return 0;
}

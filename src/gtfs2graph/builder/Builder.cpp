// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_set>
#include <vector>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "gtfs2graph/builder/Builder.h"
#include "gtfs2graph/graph/BuildGraph.h"
#include "gtfs2graph/graph/EdgePL.h"
#include "gtfs2graph/graph/NodePL.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/log/Log.h"

using namespace gtfs2graph;
using namespace graph;

using util::geo::Box;
using util::geo::DBox;
using util::geo::DPoint;
using util::geo::extendBox;
using util::geo::Grid;
using util::geo::Point;
using util::geo::PolyLine;
using util::geo::SharedSegments;

using graph::Edge;
using graph::Node;

using ad::cppgtfs::gtfs::Feed;
using ad::cppgtfs::gtfs::Shape;
using ad::cppgtfs::gtfs::Stop;
using ad::cppgtfs::gtfs::StopTime;
using ad::cppgtfs::gtfs::Trip;

using util::DEBUG;

namespace {

struct TripCloseStats {
  size_t stopCount = 0;
  size_t uniqueStopCount = 0;
  size_t uniqueNodeCount = 0;
  Node* firstNode = 0;
  Node* lastNode = 0;
  Edge* firstEdge = 0;
  Edge* lastEdge = 0;
  const Stop* firstStop = 0;
  const Stop* lastStop = 0;
  double minLat = 0.0;
  double minLon = 0.0;
  double maxLat = 0.0;
  double maxLon = 0.0;
  bool bboxInit = false;
  std::vector<double> speedSamples;
};

bool timeToSeconds(const ad::cppgtfs::gtfs::Time& t, int* out) {
  if (t.empty()) return false;
  *out = t.seconds();
  return true;
}

double median(std::vector<double> values) {
  if (values.empty()) return 0.0;
  std::sort(values.begin(), values.end());
  size_t mid = values.size() / 2;
  if (values.size() % 2 == 0) {
    return 0.5 * (values[mid - 1] + values[mid]);
  }
  return values[mid];
}

}  // namespace

// _____________________________________________________________________________
Builder::Builder(const config::Config* cfg) : _cfg(cfg) {}

// _____________________________________________________________________________
void Builder::consume(const Feed& f, BuildGraph* g) {
  DBox graphBox(getProjP(f.getMinLat(), f.getMinLon()),
                getProjP(f.getMaxLat(), f.getMaxLon()));

  NodeGrid ngrid(2000, 2000, graphBox);

  size_t i = 0;

  for (auto t = f.getTrips().begin(); t != f.getTrips().end(); ++t) {
    i++;
    // ignore trips with only one stop
    if (t->second->getStopTimes().size() < 2) continue;
    if (!_cfg->useMots.count(t->second->getRoute()->getType())) continue;

    auto st = t->second->getStopTimes().begin();
    TripCloseStats stats;
    stats.stopCount = t->second->getStopTimes().size();
    std::unordered_set<const Stop*> uniqueStops;
    std::unordered_set<Node*> uniqueNodes;

    auto prev = *st;
    const Edge* prevEdge = 0;
    Node* firstNode = addStop(prev.getStop(), g, &ngrid);
    uniqueStops.insert(prev.getStop());
    uniqueNodes.insert(firstNode);
    stats.firstNode = firstNode;
    stats.lastNode = firstNode;
    stats.firstStop = prev.getStop();
    stats.lastStop = prev.getStop();
    stats.minLat = prev.getStop()->getLat();
    stats.maxLat = prev.getStop()->getLat();
    stats.minLon = prev.getStop()->getLng();
    stats.maxLon = prev.getStop()->getLng();
    stats.bboxInit = true;
    ++st;

    if (i % 100 == 0)
      LOGTO(DEBUG, std::cerr) << "@ trip " << i << "/" << f.getTrips().size();

    for (; st != t->second->getStopTimes().end(); ++st) {
      const auto& cur = *st;

      Node* fromNode = getNodeByStop(g, prev.getStop());
      Node* toNode = addStop(cur.getStop(), g, &ngrid);
      stats.lastNode = toNode;
      stats.lastStop = cur.getStop();
      uniqueStops.insert(cur.getStop());
      uniqueNodes.insert(toNode);
      if (!stats.bboxInit) {
        stats.minLat = cur.getStop()->getLat();
        stats.maxLat = cur.getStop()->getLat();
        stats.minLon = cur.getStop()->getLng();
        stats.maxLon = cur.getStop()->getLng();
        stats.bboxInit = true;
      } else {
        const double curLat = static_cast<double>(cur.getStop()->getLat());
        const double curLon = static_cast<double>(cur.getStop()->getLng());
        stats.minLat = std::min(stats.minLat, curLat);
        stats.maxLat = std::max(stats.maxLat, curLat);
        stats.minLon = std::min(stats.minLon, curLon);
        stats.maxLon = std::max(stats.maxLon, curLon);
      }

      int prevSec = 0;
      int curSec = 0;
      bool prevOk = timeToSeconds(prev.getDepartureTime(), &prevSec);
      if (!prevOk) prevOk = timeToSeconds(prev.getArrivalTime(), &prevSec);
      bool curOk = timeToSeconds(cur.getArrivalTime(), &curSec);
      if (!curOk) curOk = timeToSeconds(cur.getDepartureTime(), &curSec);
      if (prevOk && curOk) {
        int dt = curSec - prevSec;
        if (dt > 0) {
          double dist =
              util::geo::haversine(prev.getStop()->getLat(),
                                   prev.getStop()->getLng(),
                                   cur.getStop()->getLat(),
                                   cur.getStop()->getLng());
          if (dist > 0.0) stats.speedSamples.push_back(dist / dt);
        }
      }

      // TODO: we should also allow this, for round-trips
      if (fromNode == toNode) {
        prev = cur;
        continue;
      }

      Edge* exE = g->getEdg(fromNode, toNode);

      if (!exE) {
        exE = g->addEdg(fromNode, toNode, EdgePL());
        exE->pl().setEdge(exE);
      }

      Node* directionNode = toNode;

      std::pair<bool, PolyLine<double>> edgeGeom;
      edgeGeom = getSubPolyLine(prev.getStop(), cur.getStop(), t->second,
                                prev.getShapeDistanceTravelled(),
                                cur.getShapeDistanceTravelled());

      if (prevEdge) {
        fromNode->pl().connOccurs(t->second->getRoute(), prevEdge, exE);
      }

      exE->pl().addTrip(t->second, edgeGeom.second, directionNode);

      if (!stats.firstEdge) stats.firstEdge = exE;
      stats.lastEdge = exE;

      prev = cur;
      prevEdge = exE;
    }

    stats.uniqueStopCount = uniqueStops.size();
    stats.uniqueNodeCount = uniqueNodes.size();

    if (_cfg->closeCircularRoutesMode !=
        config::CloseCircularRoutesMode::Never) {
      bool shouldClose = true;

      if (!stats.firstNode || !stats.lastNode || !stats.firstStop ||
          !stats.lastStop) {
        shouldClose = false;
      }

      if (shouldClose && stats.uniqueNodeCount < 2) shouldClose = false;
      if (shouldClose && stats.firstNode == stats.lastNode) shouldClose = false;

      double distClose = 0.0;
      if (shouldClose) {
        distClose = util::geo::haversine(stats.firstStop->getLat(),
                                         stats.firstStop->getLng(),
                                         stats.lastStop->getLat(),
                                         stats.lastStop->getLng());
        if (_cfg->circularMaxCloseDistanceM > 0 &&
            distClose > _cfg->circularMaxCloseDistanceM) {
          LOG(DEBUG) << "Skipping circular closure for trip "
                     << t->second->getId() << ": distance " << distClose
                     << " > max " << _cfg->circularMaxCloseDistanceM;
          shouldClose = false;
        }
      }

      if (shouldClose &&
          _cfg->closeCircularRoutesMode == config::CloseCircularRoutesMode::Auto) {
        if (stats.stopCount < _cfg->circularMinStops) {
          LOG(DEBUG) << "Skipping AUTO circular closure for trip "
                     << t->second->getId() << ": stops " << stats.stopCount
                     << " < min " << _cfg->circularMinStops;
          shouldClose = false;
        }
        if (stats.uniqueStopCount < _cfg->circularMinUniqueStops) {
          LOG(DEBUG) << "Skipping AUTO circular closure for trip "
                     << t->second->getId() << ": unique stops "
                     << stats.uniqueStopCount << " < min "
                     << _cfg->circularMinUniqueStops;
          shouldClose = false;
        }
        double diag = 0.0;
        if (stats.bboxInit) {
          diag = util::geo::haversine(stats.minLat, stats.minLon, stats.maxLat,
                                      stats.maxLon);
        }
        if (!(diag > 0.0)) {
          LOG(DEBUG) << "Skipping AUTO circular closure for trip "
                     << t->second->getId() << ": degenerate bbox";
          shouldClose = false;
        } else if (distClose >
                   _cfg->circularProximityRatio * diag) {
          LOG(DEBUG) << "Skipping AUTO circular closure for trip "
                     << t->second->getId() << ": ratio failed dist="
                     << distClose << " threshold="
                     << _cfg->circularProximityRatio * diag;
          shouldClose = false;
        }
      }

      if (shouldClose) {
        if (!stats.firstEdge || !stats.lastEdge) {
          shouldClose = false;
        }
      }

      if (shouldClose) {
        if (g->getEdg(stats.lastNode, stats.firstNode)) {
          LOG(DEBUG) << "Skipping circular closure for trip "
                     << t->second->getId() << ": edge already exists";
          shouldClose = false;
        }
      }

      if (shouldClose) {
        DPoint lastP = getProjP(stats.lastStop->getLat(),
                                stats.lastStop->getLng());
        DPoint firstP = getProjP(stats.firstStop->getLat(),
                                 stats.firstStop->getLng());
        PolyLine<double> closeGeom(lastP, firstP);

        Edge* closeEdge =
            g->addEdg(stats.lastNode, stats.firstNode, EdgePL());
        closeEdge->pl().setEdge(closeEdge);
        closeEdge->pl().addTrip(t->second, closeGeom, stats.firstNode);

        double speedMps = median(stats.speedSamples);
        if (speedMps <= 0.0) {
          speedMps = _cfg->circularDefaultSpeedKmh / 3.6;
        }
        double dtClose = distClose / speedMps;
        if (dtClose < 1.0) dtClose = 1.0;

        closeEdge->pl().setDistanceMeters(distClose);
        closeEdge->pl().setTimeSeconds(dtClose);
        closeEdge->pl().setSyntheticCircularClose(true);

        stats.lastNode->pl().connOccurs(t->second->getRoute(), stats.lastEdge,
                                        closeEdge);
        stats.firstNode->pl().connOccurs(t->second->getRoute(), closeEdge,
                                         stats.firstEdge);

        double diag = 0.0;
        if (stats.bboxInit) {
          diag = util::geo::haversine(stats.minLat, stats.minLon, stats.maxLat,
                                      stats.maxLon);
        }
        LOG(INFO) << "Closed circular trip " << t->second->getId()
                  << " dist=" << distClose << " diag=" << diag
                  << " ratio=" << _cfg->circularProximityRatio
                  << " speed_mps=" << speedMps << " dt=" << dtClose;
      }
    }
  }
}

// _____________________________________________________________________________
DPoint Builder::getProjP(double lat, double lng) const {
  return util::geo::latLngToWebMerc<double>(lat, lng);
}

// _____________________________________________________________________________
void Builder::simplify(BuildGraph* g) {
  // calculate average number of trip occurences per line
  double avg = 0;
  int c = 0;
  for (auto n : g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      for (auto& etg : *e->pl().getEdgeTripGeoms()) {
        for (auto& r : *etg.getTripsUnordered()) {
          avg += r.trips.size();
          c++;
        }
      }
    }
  }

  avg /= c;

  // try to merge both-direction edges into a single one
  // also prune edges with few trips
  for (auto n : g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      e->pl().simplify(avg * _cfg->pruneThreshold);
    }
  }

  // delete edges without a reference ETG
  std::vector<Edge*> toDel;
  for (auto n : g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (!e->pl().getRefETG()) toDel.push_back(e);
    }
  }

  for (auto e : toDel) g->delEdg(e->getFrom(), e->getTo());
}

// _____________________________________________________________________________
std::pair<bool, PolyLine<double>> Builder::getSubPolyLine(const Stop* a,
                                                          const Stop* b,
                                                          Trip* t, double distA,
                                                          double distB) {
  UNUSED(distA);
  UNUSED(distB);
  DPoint ap = getProjP(a->getLat(), a->getLng());
  DPoint bp = getProjP(b->getLat(), b->getLng());

  if (!t->getShape()) {
    return std::pair<bool, PolyLine<double>>(false, PolyLine<double>(ap, bp));
  }

  auto pl = _polyLines.find(t->getShape());
  if (pl == _polyLines.end()) {
    // generate polyline for this shape
    pl = _polyLines
             .insert(std::pair<Shape*, PolyLine<double>>(t->getShape(),
                                                         PolyLine<double>()))
             .first;

    for (const auto& sp : t->getShape()->getPoints()) {
      pl->second << getProjP(sp.lat, sp.lng);
    }
  }

  PolyLine<double> p;

  p = pl->second.getSegment(ap, bp);

  return std::pair<bool, PolyLine<double>>(true, p);
}

// _____________________________________________________________________________
Node* Builder::addStop(const Stop* curStop, BuildGraph* g, NodeGrid* grid) {
  Node* n = getNodeByStop(g, curStop);
  if (n) return n;

  DPoint p = getProjP(curStop->getLat(), curStop->getLng());

  if (n) {
    n->pl().addStop(curStop);
    _stopNodes[curStop] = n;
  } else {
    n = g->addNd(NodePL(p, curStop));
    n->pl().setNode(n);
    grid->add(n->pl().getPos(), n);
    _stopNodes[curStop] = n;
  }

  return n;
}

// _____________________________________________________________________________
Node* Builder::getNodeByStop(const BuildGraph* g, const gtfs::Stop* s) const {
  if (_stopNodes.find(s) != _stopNodes.end()) return _stopNodes.find(s)->second;

  for (const auto n : g->getNds()) {
    if (n->pl().getStops().find(const_cast<gtfs::Stop*>(s)) !=
        n->pl().getStops().end()) {
      return n;
    }
  }
  return 0;
}

// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <vector>
#include "ad/cppgtfs/gtfs/Trip.h"
#include "gtfs2graph/graph/EdgePL.h"
#include "gtfs2graph/graph/EdgeTripGeom.h"
#include "gtfs2graph/graph/NodePL.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"

using namespace gtfs2graph;
using namespace graph;
using namespace ad::cppgtfs;

using util::geo::PolyLine;

// _____________________________________________________________________________
EdgePL::EdgePL(const Edge* e) : _e(e) {}

// _____________________________________________________________________________
EdgePL::EdgePL() : _e(0) {}

// _____________________________________________________________________________
void EdgePL::setEdge(const Edge* e) { _e = e; }

// _____________________________________________________________________________
bool EdgePL::addTrip(gtfs::Trip* t, PolyLine<double> pl, Node* toNode) {
  assert(toNode == _e->getFrom() || toNode == _e->getTo());
  bool inserted = false;
  for (auto& e : _tripsContained) {
    if (e.getGeom().equals(pl, 10)) {
      e.addTrip(t, toNode, pl);
      inserted = true;
      break;
    }
  }

  if (!inserted) {
    EdgeTripGeom etg(pl, toNode);
    etg.addTrip(t, toNode);
    addEdgeTripGeom(etg);
  }

  return true;
}

// _____________________________________________________________________________
const std::vector<EdgeTripGeom>& EdgePL::getEdgeTripGeoms() const {
  return _tripsContained;
}

// _____________________________________________________________________________
std::vector<EdgeTripGeom>* EdgePL::getEdgeTripGeoms() {
  return &_tripsContained;
}

// _____________________________________________________________________________
const EdgeTripGeom* EdgePL::getRefETG() const {
  if (_tripsContained.size()) return &_tripsContained.front();
  return 0;
}

// _____________________________________________________________________________
EdgeTripGeom* EdgePL::getRefETG() {
  if (_tripsContained.size()) return &_tripsContained.front();
  return 0;
}

// _____________________________________________________________________________
void EdgePL::addEdgeTripGeom(const EdgeTripGeom& e) {
  assert(e.getGeomDir() == _e->getFrom() || e.getGeomDir() == _e->getTo());

  for (const auto& to : e.getTripsUnordered()) {
    assert(to.direction == _e->getFrom() || to.direction == _e->getTo());
  }

  _tripsContained.push_back(e);
  if (e.getGeomDir() != _e->getTo()) {
    const_cast<PolyLine<double>*>(&_tripsContained.back().getGeom())->reverse();
    _tripsContained.back().setGeomDir(_e->getTo());
  }
}

// _____________________________________________________________________________
void EdgePL::simplify(double pruneThreshold) {
  auto eit = _tripsContained.begin();
  while (eit != _tripsContained.end()) {
    auto etg = *eit;
    auto it = etg.getTripsUnordered()->begin();
    while (it != etg.getTripsUnordered()->end()) {
      if (it->trips.size() < pruneThreshold) {
        it = etg.getTripsUnordered()->erase(it);
      } else {
        it++;
      }
    }

    // delete if it runs empty
    if (etg.getTripsUnordered()->size() == 0) {
      eit = _tripsContained.erase(eit);
    } else {
      eit++;
    }
  }


  combineIncludedGeoms();
  averageCombineGeom();
}

// _____________________________________________________________________________
void EdgePL::averageCombineGeom() {
  if (_tripsContained.size() < 2) {
    return;
  }

  std::vector<const PolyLine<double>*> lines;

  for (auto& et : _tripsContained) {
    assert(et.getGeomDir() == _e->getTo());
    lines.push_back(&et.getGeom());
  }

  PolyLine<double> pl = PolyLine<double>::average(lines);

  EdgeTripGeom combined(pl, _e->getTo());

  for (auto& et : _tripsContained) {
    for (auto& r : *et.getTripsUnordered()) {
      for (auto& t : r.trips) {
        combined.addTrip(t, r.direction);
      }
    }
  }

  _tripsContained.clear();
  _tripsContained.push_back(combined);
}

// _____________________________________________________________________________
void EdgePL::combineIncludedGeoms() {
  if (_tripsContained.size() < 2) {
    return;
  }

  for (auto et = _tripsContained.begin(); et != _tripsContained.end();) {
    bool combined = false;
    for (auto& toCheckAgainst : _tripsContained) {
      if (toCheckAgainst.getGeom().getLength() > et->getGeom().getLength() &&
          toCheckAgainst.getGeom().contains(et->getGeom(), 50) &&
          !et->getGeom().contains(toCheckAgainst.getGeom(), 50)) {
        for (auto& r : *et->getTripsUnordered()) {
          for (auto& t : r.trips) {
            toCheckAgainst.addTrip(t, r.direction);
          }
        }
        combined = true;
        break;
      }
    }
    if (combined) {
      // delete the old EdgeTripGeom
      et = _tripsContained.erase(et);
    } else {
      et++;
    }
  }
}

// _____________________________________________________________________________
const util::geo::DLine* EdgePL::getGeom() const {
  if (!getRefETG()) return 0;
  return &getRefETG()->getGeom().getLine();
}

// _____________________________________________________________________________
util::json::Dict EdgePL::getAttrs() const {
  util::json::Dict obj;
  util::json::Array lines;

  if (!getRefETG()) return obj;

  for (auto r : getRefETG()->getTripsUnordered()) {
    util::json::Dict route;
    route["id"] = util::toString(r.route);
    route["label"] = r.route->getShortName();
    route["color"] = r.route->getColorString();

    route["direction"] = util::toString(r.direction);

    lines.push_back(route);
  }
  obj["lines"] = lines;

  return obj;
}

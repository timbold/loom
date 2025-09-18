// Copyright 2016
// Author: Patrick Brosi

#include <fstream>
#include <set>
#include <string>

#include "ad/cppgtfs/Parser.h"
#include "ad/cppgtfs/gtfs/Feed.h"
#include "ad/cppgtfs/gtfs/Route.h"
#include "ad/cppgtfs/gtfs/flat/Route.h"
#include "gtfs2graph/builder/Builder.h"
#include "gtfs2graph/config/GraphBuilderConfig.h"
#include "gtfs2graph/graph/BuildGraph.h"
#include "gtfs2graph/graph/EdgePL.h"
#include "gtfs2graph/graph/NodePL.h"
#include "util/Misc.h"
#include "util/String.h"
#include "util/json/Writer.h"

namespace {

void writeGtfsFile(const std::string& baseDir, const std::string& name,
                   const std::string& content) {
  std::ofstream out(baseDir + "/" + name, std::ios::out | std::ios::trunc);
  TEST(out.good());
  out << content;
  TEST(out.good());
}

bool nodeHasStop(const gtfs2graph::graph::Node* node,
                 const std::string& stopId) {
  for (const auto stop : node->pl().getStops()) {
    if (stop->getId() == stopId) return true;
  }
  return false;
}

}  // namespace

// _____________________________________________________________________________
int main(int argc, char** argv) {
  UNUSED(argc);
  UNUSED(argv);

  char tmpTemplate[] = "/tmp/gtfsXXXXXX";
  char* dirName = mkdtemp(tmpTemplate);
  TEST(dirName != nullptr);
  std::string baseDir(dirName);

  writeGtfsFile(baseDir, "agency.txt",
                "agency_id,agency_name,agency_url,agency_timezone\n"
                "AG,Agency,http://example.com,Europe/Berlin\n");
  writeGtfsFile(baseDir, "routes.txt",
                "route_id,agency_id,route_short_name,route_long_name,route_type\n"
                "R1,AG,1,Sample Route,3\n");
  writeGtfsFile(baseDir, "trips.txt",
                "route_id,service_id,trip_id\n"
                "R1,WK,TRIP1\n");
  writeGtfsFile(baseDir, "stop_times.txt",
                "trip_id,arrival_time,departure_time,stop_id,stop_sequence\n"
                "TRIP1,08:00:00,08:00:00,STOP_A,1\n"
                "TRIP1,08:05:00,08:05:00,STOP_B,2\n"
                "TRIP1,08:10:00,08:10:00,STOP_C,3\n"
                "TRIP1,08:15:00,08:15:00,STOP_B,4\n"
                "TRIP1,08:20:00,08:20:00,STOP_D,5\n");
  writeGtfsFile(baseDir, "stops.txt",
                "stop_id,stop_name,stop_lat,stop_lon\n"
                "STOP_A,A,0.0,0.0\n"
                "STOP_B,B,0.1,0.0\n"
                "STOP_C,C,0.1,0.1\n"
                "STOP_D,D,0.0,0.1\n");
  writeGtfsFile(baseDir, "calendar.txt",
                "service_id,monday,tuesday,wednesday,thursday,friday,saturday,sunday,start_date,end_date\n"
                "WK,1,1,1,1,1,0,0,20220101,20221231\n");

  ad::cppgtfs::gtfs::Feed feed;
  ad::cppgtfs::Parser parser(baseDir);
  parser.parse(&feed);

  gtfs2graph::config::Config cfg;
  cfg.pruneThreshold = 0.0;
  for (auto mot :
       ad::cppgtfs::gtfs::flat::Route::getTypesFromString("bus")) {
    cfg.useMots.insert(mot);
  }

  gtfs2graph::Builder builder(&cfg);
  gtfs2graph::graph::BuildGraph graph;
  builder.consume(feed, &graph);
  builder.simplify(&graph);

  const gtfs2graph::graph::Node* nodeB = nullptr;
  const gtfs2graph::graph::Node* nodeC = nullptr;
  for (auto node : graph.getNds()) {
    if (!nodeB && nodeHasStop(node, "STOP_B")) nodeB = node;
    if (!nodeC && nodeHasStop(node, "STOP_C")) nodeC = node;
  }

  TEST(nodeB != nullptr);
  TEST(nodeC != nullptr);

  const gtfs2graph::graph::Edge* bcEdge = nullptr;
  for (auto node : graph.getNds()) {
    for (auto edge : node->getAdjList()) {
      bool connectsB = nodeHasStop(edge->getFrom(), "STOP_B") ||
                       nodeHasStop(edge->getTo(), "STOP_B");
      bool connectsC = nodeHasStop(edge->getFrom(), "STOP_C") ||
                       nodeHasStop(edge->getTo(), "STOP_C");
      if (connectsB && connectsC) {
        bcEdge = edge;
        break;
      }
    }
    if (bcEdge) break;
  }

  TEST(bcEdge != nullptr);

  auto attrs = bcEdge->pl().getAttrs();
  auto linesIt = attrs.find("lines");
  TEST(linesIt != attrs.end());
  const auto& linesVal = linesIt->second;
  TEST(linesVal.type, ==, util::json::Val::ARRAY);
  TEST(linesVal.arr.size(), ==, 2u);

  std::string expectedRouteId;
  std::set<std::string> directions;
  for (const auto& entry : linesVal.arr) {
    TEST(entry.type, ==, util::json::Val::DICT);
    const auto& dict = entry.dict;
    auto idIt = dict.find("id");
    TEST(idIt != dict.end());
    TEST(idIt->second.type, ==, util::json::Val::STRING);
    if (expectedRouteId.empty()) {
      expectedRouteId = idIt->second.str;
    } else {
      TEST(idIt->second.str, ==, expectedRouteId);
    }

    auto dirIt = dict.find("direction");
    TEST(dirIt != dict.end());
    TEST(dirIt->second.type, ==, util::json::Val::STRING);
    directions.insert(dirIt->second.str);
  }

  TEST(directions.size(), ==, 2u);
  std::string nodeBId = util::toString(nodeB);
  std::string nodeCId = util::toString(nodeC);
  TEST(directions.count(nodeBId), ==, 1u);
  TEST(directions.count(nodeCId), ==, 1u);

  return 0;
}

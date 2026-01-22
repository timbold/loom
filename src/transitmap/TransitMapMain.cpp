// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>

#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include "shared/rendergraph/Penalties.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/ConfigReader.cpp"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/GraphBuilder.h"
#include "transitmap/output/SvgRenderer.h"
#include "util/log/Log.h"

using shared::linegraph::LineGraph;
using shared::rendergraph::RenderGraph;
using transitmapper::graph::GraphBuilder;

using util::DEBUG;
using util::ERROR;

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // initialize randomness
  srand(time(NULL) + rand());

  transitmapper::config::Config cfg;

  transitmapper::config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  T_START(TIMER);

  GraphBuilder b(&cfg);

  LOGTO(DEBUG, std::cerr) << "Reading graph...";

  if (cfg.renderMethod == "svg") {
    RenderGraph g(cfg.lineWidth, cfg.outlineWidth, cfg.lineSpacing);
    if (cfg.fromDot)
      g.readFromDot(&std::cin);
    else
      g.readFromJson(&std::cin);

    if (cfg.randomColors) g.fillMissingColors();

    // snap orphan stations
    g.snapOrphanStations();

    g.contractStrayNds();
    g.smooth(cfg.inputSmoothing);
    b.writeNodeFronts(&g);
    b.expandOverlappinFronts(&g);
    g.createMetaNodes();

    if (true) {
      b.dropOverlappingStations(&g);
      g.contractStrayNds();
      b.expandOverlappinFronts(&g);
      g.createMetaNodes();
    }

    LOGTO(DEBUG, std::cerr) << "Outputting to SVG ...";
    transitmapper::output::SvgRenderer svgOut(&std::cout, &cfg);
    svgOut.print(g);
  } else {
    LOG(ERROR) << "Unknown render method " << cfg.renderMethod;
    exit(1);
  }

  double took = T_STOP(TIMER);

  if (cfg.writeStats) {
    util::json::Writer wr(&std::cout);
    wr.obj();
    wr.keyVal("time", took);
    wr.closeAll();
  }

  return (0);
}

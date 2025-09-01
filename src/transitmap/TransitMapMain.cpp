// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#ifndef _WIN32
#include <unistd.h>
#endif

#include <iostream>
#include <string>
#include <random>

#include "shared/rendergraph/Penalties.h"
#include "shared/rendergraph/RenderGraph.h"
#include "shared/rendergraph/Landmark.h"
#include "transitmap/config/ConfigReader.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/GraphBuilder.h"
#include "transitmap/output/MvtRenderer.h"
#include "transitmap/output/SvgRenderer.h"
#include "transitmap/util/String.h"
#include "util/log/Log.h"

using shared::linegraph::LineGraph;
using shared::rendergraph::RenderGraph;
using transitmapper::graph::GraphBuilder;

// _____________________________________________________________________________
int main(int argc, char **argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  std::random_device rd;
  std::mt19937 rng(rd());
  srand(static_cast<unsigned>(rng()));

  transitmapper::config::Config cfg;

  transitmapper::config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  T_START(TIMER);

  GraphBuilder b(&cfg);

  LOGTO(DEBUG, std::cerr) << "Reading graph...";

  if (cfg.renderMethod == "mvt") {
#ifdef PROTOBUF_FOUND
    LineGraph lg;
    if (cfg.fromDot)
      lg.readFromDot(&std::cin);
    else
      lg.readFromJson(&std::cin);

    if (cfg.randomColors)
      lg.fillMissingColors();

    // snap orphan stations
    lg.snapOrphanStations();

    for (size_t z : cfg.mvtZooms) {
      double lWidth = cfg.lineWidth;
      double lSpacing = cfg.lineSpacing;
      double lOutlineWidth = cfg.outlineWidth;

      lWidth *= 156543.0 / (1 << z);
      lSpacing *= 156543.0 / (1 << z);
      lOutlineWidth *= 156543.0 / (1 << z);

      RenderGraph g(lg, lWidth, lOutlineWidth, lSpacing);

      g.contractStrayNds();
      g.smooth(cfg.inputSmoothing);
      b.writeNodeFronts(&g);
      b.expandOverlappinFronts(&g);

      g.createMetaNodes();

      // avoid overlapping stations
      b.dropOverlappingStations(&g);
      g.contractStrayNds();
      b.expandOverlappinFronts(&g);
      g.createMetaNodes();

      if (!cfg.meStation.empty()) {
        for (auto n : g.getNds()) {
          if (!n->pl().stops().size()) continue;
          const auto &st = n->pl().stops().front();
          if (util::sanitizeStationLabel(st.name) == cfg.meStation) {
            cfg.meLandmark.coord = st.pos;
            cfg.meLandmark.color = cfg.meStationFill;
            cfg.renderMe = true;
            break;
          }
        }
      }

      LOGTO(DEBUG, std::cerr) << "Outputting to MVT ...";
      transitmapper::output::MvtRenderer mvtOut(&cfg, z);
      mvtOut.print(g);
    }
#else
    LOG(ERROR) << "transitmap was not compiled with protocol buffers support, "
                  "cannot use render method "
               << cfg.renderMethod;
    exit(1);
#endif
  } else if (cfg.renderMethod == "svg") {
    RenderGraph g(cfg.lineWidth, cfg.outlineWidth, cfg.lineSpacing);
    if (cfg.fromDot)
      g.readFromDot(&std::cin);
    else
      g.readFromJson(&std::cin);

    if (cfg.randomColors)
      g.fillMissingColors();

    // snap orphan stations
    g.snapOrphanStations();

    g.contractStrayNds();
    g.smooth(cfg.inputSmoothing);
    b.writeNodeFronts(&g);
    b.expandOverlappinFronts(&g);
    g.createMetaNodes();

    b.dropOverlappingStations(&g);
    g.contractStrayNds();
    b.expandOverlappinFronts(&g);
    g.createMetaNodes();

    if (!cfg.meStation.empty()) {
      for (auto n : g.getNds()) {
        if (!n->pl().stops().size()) continue;
        const auto &st = n->pl().stops().front();
        if (util::sanitizeStationLabel(st.name) == cfg.meStation) {
          cfg.meLandmark.coord = st.pos;
          cfg.meLandmark.color = cfg.meStationFill;
          cfg.renderMe = true;
          break;
        }
      }
    }

    // Attach landmarks.
    for (const auto &lmCfg : cfg.landmarks) {
      shared::rendergraph::Landmark lm;
      if (!lmCfg.iconPath.empty()) {
        // Landmark with an icon - keep existing SVG loading behavior.
        lm = lmCfg;
      } else if (!lmCfg.label.empty()) {
        // Landmark with a label, color and size but no icon.
        lm.label = lmCfg.label;
        lm.color = lmCfg.color;
        lm.size = lmCfg.size;
        lm.coord = lmCfg.coord;
      } else {
        continue;
      }
      g.addLandmark(lm);
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

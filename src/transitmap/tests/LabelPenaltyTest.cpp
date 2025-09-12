#include "transitmap/tests/LabelPenaltyTest.h"

#include "shared/linegraph/LineNodePL.h"
#include "transitmap/label/Labeller.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"
#include "util/geo/Line.h"

using transitmapper::config::Config;
using transitmapper::label::StationLabel;
using transitmapper::label::Overlaps;
using util::geo::PolyLine;
using util::geo::MultiLine;
using shared::linegraph::Station;

void LabelPenaltyTest::run() {
  Config cfgA;
  Config cfgB;
  cfgA.sidePenaltyWeight = 1.0;
  cfgB.sidePenaltyWeight = 5.0;
  cfgA.orientationPenalties = {0, 100, 0, 0, 0, 0, 0, 0};
  cfgB.orientationPenalties = cfgA.orientationPenalties;

  PolyLine<double> geom;
  MultiLine<double> band;
  Overlaps ov{0, 0, 0, 0, 0};
  Station st("id", "name", util::geo::DPoint());

  double sideA = 2 * cfgA.sidePenaltyWeight;
  double sideB = 2 * cfgB.sidePenaltyWeight;

  StationLabel lblA(geom, band, 10, false, 0, 0, ov, sideA,
                    cfgA.stationLineOverlapPenalty, 0, 0,
                    &cfgA.orientationPenalties, st);
  StationLabel lblB(geom, band, 10, false, 0, 0, ov, sideB,
                    cfgB.stationLineOverlapPenalty, 0, 0,
                    &cfgB.orientationPenalties, st);
  TEST(lblB.getPen() > lblA.getPen());

  StationLabel lblC(geom, band, 10, false, 1, 0, ov, sideA,
                    cfgA.stationLineOverlapPenalty, 0, 0,
                    &cfgA.orientationPenalties, st);
  TEST(lblC.getPen() > lblA.getPen());
}

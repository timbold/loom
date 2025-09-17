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
                    cfgA.stationLineOverlapPenalty, 0, 0, false,
                    cfgA.clusterPenScale, cfgA.outsidePenalty,
                    &cfgA.orientationPenalties, st);
  StationLabel lblB(geom, band, 10, false, 0, 0, ov, sideB,
                    cfgB.stationLineOverlapPenalty, 0, 0, false,
                    cfgB.clusterPenScale, cfgB.outsidePenalty,
                    &cfgB.orientationPenalties, st);
  TEST(lblB.getPen() > lblA.getPen());

  StationLabel lblC(geom, band, 10, false, 1, 0, ov, sideA,
                    cfgA.stationLineOverlapPenalty, 0, 0, false,
                    cfgA.clusterPenScale, cfgA.outsidePenalty,
                    &cfgA.orientationPenalties, st);
  TEST(lblC.getPen() > lblA.getPen());

  StationLabel lblFar(geom, band, 10, false, 0, 0, ov, sideA,
                      cfgA.stationLineOverlapPenalty, 0,
                      cfgA.stationLabelFarCrowdPenalty, false,
                      cfgA.clusterPenScale, cfgA.outsidePenalty,
                      &cfgA.orientationPenalties, st);
  TEST(lblFar.getPen() > lblA.getPen());

  // crowding penalty scale
  Config cfgC;
  Config cfgD;
  cfgC.clusterPenScale = 1.0;
  cfgD.clusterPenScale = 5.0;
  StationLabel lblD(geom, band, 10, false, 0, 0, ov, 0,
                    cfgC.stationLineOverlapPenalty, 1.0, 0, false,
                    cfgC.clusterPenScale, cfgC.outsidePenalty,
                    &cfgC.orientationPenalties, st);
  StationLabel lblE(geom, band, 10, false, 0, 0, ov, 0,
                    cfgD.stationLineOverlapPenalty, 1.0, 0, false,
                    cfgD.clusterPenScale, cfgD.outsidePenalty,
                    &cfgD.orientationPenalties, st);
  TEST(lblE.getPen() > lblD.getPen());

  // outside penalty bonus/penalty
  Config cfgOutPen;
  Config cfgOutBon;
  cfgOutPen.outsidePenalty = 5.0;
  cfgOutBon.outsidePenalty = -5.0;
  StationLabel lblOutPen(geom, band, 10, false, 0, 0, ov, 0,
                         cfgOutPen.stationLineOverlapPenalty, 0, 0, true,
                         cfgOutPen.clusterPenScale, cfgOutPen.outsidePenalty,
                         &cfgOutPen.orientationPenalties, st);
  StationLabel lblOutBon(geom, band, 10, false, 0, 0, ov, 0,
                         cfgOutBon.stationLineOverlapPenalty, 0, 0, true,
                         cfgOutBon.clusterPenScale, cfgOutBon.outsidePenalty,
                         &cfgOutBon.orientationPenalties, st);
  TEST(lblOutPen.getPen() > lblOutBon.getPen());
}

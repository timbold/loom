// Copyright 2016
// Author: Patrick Brosi

#include "util/Misc.h"
#include "transitmap/tests/SanitizeSvgTest.h"
#include "transitmap/tests/DirMarkerTest.h"
#include "transitmap/tests/DropOverlappingStationsTest.h"
#include "transitmap/tests/ArrowHeadDirectionTest.h"
#include "transitmap/tests/BgMapTest.h"
#include "transitmap/tests/ConfigParseTest.h"
#include "transitmap/tests/LabelPenaltyTest.h"
#include "transitmap/tests/LandmarkProjectionTest.h"
#include "transitmap/tests/LandmarkSizeTest.h"
#include "transitmap/tests/LandmarkDisplacementTest.h"
#include "transitmap/tests/MeBadgeSizingTest.h"
#include "transitmap/tests/StationFarCrowdTest.h"
#include "transitmap/tests/StationLabelOptimizerTest.h"
#include "transitmap/tests/TerminusLabelPlacementTest.h"
#include "transitmap/tests/TerminusReverseTest.h"

// _____________________________________________________________________________
int main(int argc, char** argv) {
  UNUSED(argc);
  UNUSED(argv);
  SanitizeSvgTest st;
  st.run();
  DirMarkerTest dmt;
  dmt.run();
  DropOverlappingStationsTest dost;
  dost.run();
  ArrowHeadDirectionTest ahdt;
  ahdt.run();
  BgMapTest bgmt;
  bgmt.run();
  ConfigParseTest cpt;
  cpt.run();
  LabelPenaltyTest lbt;
  lbt.run();
  LandmarkProjectionTest lpt;
  lpt.run();
  LandmarkSizeTest lst;
  lst.run();
  LandmarkDisplacementTest ldt;
  ldt.run();
  MeBadgeSizingTest mbst;
  mbst.run();
  StationFarCrowdTest sfct;
  sfct.run();
  StationLabelOptimizerTest slot;
  slot.run();
  TerminusLabelPlacementTest tlpt;
  tlpt.run();
  TerminusReverseTest trt;
  trt.run();
  return 0;
}

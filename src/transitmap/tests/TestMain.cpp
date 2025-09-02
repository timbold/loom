// Copyright 2016
// Author: Patrick Brosi

#include "util/Misc.h"
#include "transitmap/tests/SanitizeSvgTest.h"
#include "transitmap/tests/DirMarkerTest.h"

// _____________________________________________________________________________
int main(int argc, char** argv) {
  UNUSED(argc);
  UNUSED(argv);
  SanitizeSvgTest st;
  st.run();
  DirMarkerTest dmt;
  dmt.run();
  return 0;
}

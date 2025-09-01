// Copyright 2016
// Author: Patrick Brosi

#include <string>

#include "transitmap/tests/SanitizeSvgTest.h"
#include "util/Misc.h"

// Forward declaration of the sanitizer
bool sanitizeSvg(std::string &s);

// _____________________________________________________________________________
void SanitizeSvgTest::run() {
  {
    std::string svg =
        "<svg><rect style=\"fill:url('javascript:alert(1)')\"/></svg>";
    bool unsafe = sanitizeSvg(svg);
    TEST(unsafe, ==, true);
    TEST(svg.find("style"), ==, std::string::npos);
  }

  {
    std::string svg =
        "<svg><image xlink:href=\"data:text/html,<script>alert(1)</script>\"/></svg>";
    bool unsafe = sanitizeSvg(svg);
    TEST(unsafe, ==, true);
    TEST(svg.find("data:"), ==, std::string::npos);
    TEST(svg.find("script"), ==, std::string::npos);
  }

  {
    std::string svg =
        "<svg><style>*{display:none}</style><rect /></svg>";
    bool unsafe = sanitizeSvg(svg);
    TEST(unsafe, ==, true);
    TEST(svg.find("<style"), ==, std::string::npos);
  }
}


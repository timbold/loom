// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_LOG_LOG_H_
#define UTIL_LOG_LOG_H_

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace util {
enum Level : char { ERROR = 0, WARN = 1, INFO = 2, DEBUG = 3, VDEBUG = 4 };

#ifndef LOGLEVEL
#define LOGLEVEL 2
#endif
#ifndef UTIL_LOGLVL
#define UTIL_LOGLVL LOGLEVEL
#endif

// compiler will optimize statement away if x > LOGLEVEL
#define LOG(x) if (x <= LOGLEVEL) util::Log<x>().log()
#define LOGTO(x, os) if (x <= LOGLEVEL) util::Log<x>(&os).log()

using std::setfill;
using std::setw;
using std::chrono::system_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::time_point_cast;

const static char* LOGS[] = {"ERROR", "WARN ", "INFO ", "DEBUG", "DEBUG"};

template <char LVL>
class Log {
 public:
  Log() { if (LVL < INFO) os = &std::cerr; else os = &std::cout; }
  explicit Log(std::ostream* s) { os = s; }
  ~Log() { buf << std::endl; (*os) << buf.str(); }
  std::ostream& log() { return ts() << LOGS[(size_t)LVL] << ": "; }

 private:
  std::ostream* os;
  std::ostringstream buf;
  std::ostream& ts() {
    char tl[20];
    auto n = system_clock::now();
    time_t tt = system_clock::to_time_t(n);
    int m = duration_cast<milliseconds>(n-time_point_cast<seconds>(n)).count();
    struct tm t = *localtime(&tt);
    strftime(tl, 20, "%Y-%m-%d %H:%M:%S", &t);
    return buf << "[" << tl << "." << setfill('0') << setw(3) << m << "] ";
  }
};
}

using util::ERROR;
using util::WARN;
using util::INFO;
using util::DEBUG;
using util::VDEBUG;

#endif  // UTIL_LOG_LOG_H_

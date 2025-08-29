// Utility string functions for LOOM.
// Provides basic helpers including station label sanitization.

#ifndef UTIL_STRING_H_
#define UTIL_STRING_H_

#include <algorithm>
#include <codecvt>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

namespace util {

inline std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  std::string item;
  std::stringstream ss(s);
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

inline void replaceAll(std::string &str, const std::string &from,
                       const std::string &to) {
  if (from.empty()) return;
  size_t start = 0;
  while ((start = str.find(from, start)) != std::string::npos) {
    str.replace(start, from.length(), to);
    start += to.length();
  }
}

inline std::string trimCopy(const std::string &s) {
  const std::string WHITESPACE = " \n\r\t";
  auto start = s.find_first_not_of(WHITESPACE);
  if (start == std::string::npos) return "";
  auto end = s.find_last_not_of(WHITESPACE);
  return s.substr(start, end - start + 1);
}

inline std::string toString(double d) {
  std::ostringstream oss;
  oss << d;
  return oss.str();
}

inline std::string sanitizeStationLabel(const std::string &in) {
  std::wstring_convert<std::codecvt_utf8<wchar_t>> conv;
  std::wstring ws = conv.from_bytes(in);
  std::wstring out;
  bool prevUnderscore = false;
  std::locale loc;
  for (wchar_t c : ws) {
    if (std::iswalnum(c, loc)) {
      out.push_back(c);
      prevUnderscore = false;
    } else if (!prevUnderscore) {
      out.push_back(L'_');
      prevUnderscore = true;
    }
  }
  while (!out.empty() && out.front() == L'_') out.erase(out.begin());
  while (!out.empty() && out.back() == L'_') out.pop_back();
  return conv.to_bytes(out);
}

}  // namespace util

#endif  // UTIL_STRING_H_

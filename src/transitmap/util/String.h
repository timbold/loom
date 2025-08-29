#ifndef TRANSITMAP_UTIL_STRING_H_
#define TRANSITMAP_UTIL_STRING_H_

#include <codecvt>
#include <cwctype>
#include <locale>
#include <string>

namespace util {

inline std::string sanitizeStationLabel(const std::string& input) {
  std::wstring_convert<std::codecvt_utf8<wchar_t>> conv;
  std::wstring ws = conv.from_bytes(input);
  std::wstring out;
  bool prevUnderscore = false;
  std::locale loc;
  for (wchar_t c : ws) {
    if (std::isalnum(c, loc)) {
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

#endif  // TRANSITMAP_UTIL_STRING_H_

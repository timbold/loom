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
  std::locale loc("");
  std::wstring out;
  for (wchar_t c : ws) {
    if (std::iswalnum(c, loc)) {
      out.push_back(std::towlower(c, loc));
    }
  }
  return conv.to_bytes(out);
}

}  // namespace util

#endif  // TRANSITMAP_UTIL_STRING_H_

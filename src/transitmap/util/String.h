#ifndef TRANSITMAP_UTIL_STRING_H_
#define TRANSITMAP_UTIL_STRING_H_

#include <string>
#include <unicode/unistr.h>  // icu::UnicodeString
#include <unicode/locid.h>   // icu::Locale
#include <unicode/utf16.h>   // U16_LENGTH

namespace util {

inline std::string sanitizeStationLabel(const std::string& input) {
  // 1) UTF-8 -> ICU string
  icu::UnicodeString u = icu::UnicodeString::fromUTF8(input);

  // 2) Replace runs of disallowed with single underscore
  icu::UnicodeString tmp;
  tmp.remove();
  bool prevUnderscore = false;

  for (int32_t i = 0; i < u.length();) {
    UChar32 c = u.char32At(i);
    i += U16_LENGTH(c);

    const bool ascii_digit      = (c >= 0x30 && c <= 0x39);   // 0-9
    const bool ascii_upper      = (c >= 0x41 && c <= 0x5A);   // A-Z
    const bool ascii_lower      = (c >= 0x61 && c <= 0x7A);   // a-z
    const bool bmp_non_ascii    = (c >= 0x80  && c <= 0xFFFF); // U+0080..U+FFFF
    const bool allowed = ascii_digit || ascii_upper || ascii_lower || bmp_non_ascii;

    if (allowed) {
      tmp.append(c);
      prevUnderscore = false;
    } else {
      if (!prevUnderscore) {
        tmp.append(UChar(0x5F)); // '_'
        prevUnderscore = true;
      }
    }
  }

  // 3) Trim leading/trailing underscores
  int32_t start = 0, end = tmp.length();
  while (start < end && tmp[start] == 0x5F) ++start;
  while (end > start && tmp[end - 1] == 0x5F) --end;
  icu::UnicodeString trimmed = tmp.tempSubStringBetween(start, end);

  // 4) Unicode lowercasing with ROOT locale (locale-independent; matches Python str.lower)
  icu::UnicodeString lowered = trimmed.toLower(icu::Locale::getRoot());

  // 5) ICU -> UTF-8
  std::string result;
  lowered.toUTF8String(result);
  return result;
}

}  // namespace util

#endif  // TRANSITMAP_UTIL_STRING_H_

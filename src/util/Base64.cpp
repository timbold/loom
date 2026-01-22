// Copyright 2026
// Base64 helpers

#include "util/Base64.h"

namespace util {

static const char* BASE64_ALPHABET =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

std::string base64Encode(const unsigned char* data, size_t len) {
  std::string out;
  out.reserve(((len + 2) / 3) * 4);

  size_t i = 0;
  while (i + 2 < len) {
    unsigned int n = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
    out.push_back(BASE64_ALPHABET[(n >> 18) & 63]);
    out.push_back(BASE64_ALPHABET[(n >> 12) & 63]);
    out.push_back(BASE64_ALPHABET[(n >> 6) & 63]);
    out.push_back(BASE64_ALPHABET[n & 63]);
    i += 3;
  }

  if (i < len) {
    unsigned int n = data[i] << 16;
    if (i + 1 < len) n |= data[i + 1] << 8;

    out.push_back(BASE64_ALPHABET[(n >> 18) & 63]);
    out.push_back(BASE64_ALPHABET[(n >> 12) & 63]);

    if (i + 1 < len) {
      out.push_back(BASE64_ALPHABET[(n >> 6) & 63]);
    } else {
      out.push_back('=');
    }

    out.push_back('=');
  }

  return out;
}

std::string base64Encode(const std::vector<unsigned char>& data) {
  if (data.empty()) return "";
  return base64Encode(data.data(), data.size());
}

}  // namespace util

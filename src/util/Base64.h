// Copyright 2026
// Base64 helpers

#ifndef UTIL_BASE64_H_
#define UTIL_BASE64_H_

#include <cstddef>
#include <string>
#include <vector>

namespace util {

std::string base64Encode(const unsigned char* data, size_t len);
std::string base64Encode(const std::vector<unsigned char>& data);

}  // namespace util

#endif  // UTIL_BASE64_H_

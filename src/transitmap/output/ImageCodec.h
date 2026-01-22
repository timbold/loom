// Copyright 2026
// Image codec helpers

#ifndef TRANSITMAP_OUTPUT_IMAGECODEC_H_
#define TRANSITMAP_OUTPUT_IMAGECODEC_H_

#include <cstddef>
#include <string>
#include <vector>

namespace transitmapper {
namespace output {

struct RgbaImage {
  int width = 0;
  int height = 0;
  std::vector<unsigned char> rgba;
};

bool decodeImageToRgba(const std::string& format,
                       const unsigned char* data,
                       size_t len,
                       RgbaImage* out);

bool encodeRgbaToPng(const RgbaImage& image, std::vector<unsigned char>* out);

}  // namespace output
}  // namespace transitmapper

#endif  // TRANSITMAP_OUTPUT_IMAGECODEC_H_

// Copyright 2026
// Image codec helpers

#include "transitmap/output/ImageCodec.h"

#include <png.h>
#include <jpeglib.h>

#include <algorithm>
#include <csetjmp>
#include <cstring>
#include <cctype>
#include <string>
#include <vector>

#include "util/log/Log.h"
#include "util/Misc.h"

using util::ERROR;

namespace transitmapper {
namespace output {

namespace {

struct JpegErrorManager {
  jpeg_error_mgr pub;
  jmp_buf setjmp_buffer;
};

void jpegErrorExit(j_common_ptr cinfo) {
  JpegErrorManager* err =
      reinterpret_cast<JpegErrorManager*>(cinfo->err);
  longjmp(err->setjmp_buffer, 1);
}

struct PngReadState {
  const unsigned char* data = nullptr;
  size_t size = 0;
  size_t offset = 0;
};

void pngReadCallback(png_structp png_ptr, png_bytep outBytes,
                     png_size_t byteCountToRead) {
  PngReadState* state =
      reinterpret_cast<PngReadState*>(png_get_io_ptr(png_ptr));
  if (!state || state->offset + byteCountToRead > state->size) {
    png_error(png_ptr, "png read out of bounds");
    return;
  }
  std::memcpy(outBytes, state->data + state->offset, byteCountToRead);
  state->offset += byteCountToRead;
}

struct PngWriteState {
  std::vector<unsigned char>* out = nullptr;
};

void pngWriteCallback(png_structp png_ptr, png_bytep data,
                      png_size_t length) {
  PngWriteState* state =
      reinterpret_cast<PngWriteState*>(png_get_io_ptr(png_ptr));
  if (!state || !state->out) return;
  state->out->insert(state->out->end(), data, data + length);
}

void pngFlushCallback(png_structp png_ptr) { UNUSED(png_ptr); }

bool decodePngToRgba(const unsigned char* data, size_t len, RgbaImage* out) {
  if (!data || len == 0 || !out) return false;

  png_structp png_ptr =
      png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
  if (!png_ptr) return false;
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_read_struct(&png_ptr, nullptr, nullptr);
    return false;
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);
    return false;
  }

  PngReadState state;
  state.data = data;
  state.size = len;
  state.offset = 0;
  png_set_read_fn(png_ptr, &state, pngReadCallback);

  png_read_info(png_ptr, info_ptr);

  png_uint_32 width = 0;
  png_uint_32 height = 0;
  int bit_depth = 0;
  int color_type = 0;
  int interlace = 0;
  int compression = 0;
  int filter = 0;

  png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
               &interlace, &compression, &filter);

  if (bit_depth == 16) png_set_strip_16(png_ptr);
  if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png_ptr);
  if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
    png_set_expand_gray_1_2_4_to_8(png_ptr);
  if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
    png_set_tRNS_to_alpha(png_ptr);
  if (color_type == PNG_COLOR_TYPE_RGB ||
      color_type == PNG_COLOR_TYPE_GRAY ||
      color_type == PNG_COLOR_TYPE_PALETTE) {
    png_set_filler(png_ptr, 0xFF, PNG_FILLER_AFTER);
  }
  if (color_type == PNG_COLOR_TYPE_GRAY ||
      color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
    png_set_gray_to_rgb(png_ptr);
  }

  png_read_update_info(png_ptr, info_ptr);

  out->width = static_cast<int>(width);
  out->height = static_cast<int>(height);
  out->rgba.resize(out->width * out->height * 4);

  std::vector<png_bytep> rows(out->height);
  for (int y = 0; y < out->height; ++y) {
    rows[y] = out->rgba.data() + y * out->width * 4;
  }

  png_read_image(png_ptr, rows.data());
  png_read_end(png_ptr, nullptr);
  png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);

  return true;
}

bool decodeJpegToRgba(const unsigned char* data, size_t len, RgbaImage* out) {
  if (!data || len == 0 || !out) return false;

  jpeg_decompress_struct cinfo;
  JpegErrorManager jerr;
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = jpegErrorExit;
  jpeg_create_decompress(&cinfo);

  if (setjmp(jerr.setjmp_buffer)) {
    jpeg_destroy_decompress(&cinfo);
    return false;
  }

  jpeg_mem_src(&cinfo, data, static_cast<unsigned long>(len));
  if (jpeg_read_header(&cinfo, TRUE) != JPEG_HEADER_OK) {
    jpeg_destroy_decompress(&cinfo);
    return false;
  }

  jpeg_start_decompress(&cinfo);
  out->width = static_cast<int>(cinfo.output_width);
  out->height = static_cast<int>(cinfo.output_height);

  const int channels = static_cast<int>(cinfo.output_components);
  std::vector<unsigned char> row(out->width * channels);
  out->rgba.resize(out->width * out->height * 4);

  while (cinfo.output_scanline < cinfo.output_height) {
    unsigned char* rowPtr = row.data();
    jpeg_read_scanlines(&cinfo, &rowPtr, 1);
    int y = static_cast<int>(cinfo.output_scanline - 1);
    for (int x = 0; x < out->width; ++x) {
      size_t src = x * channels;
      size_t dst = (y * out->width + x) * 4;
      if (channels == 3) {
        out->rgba[dst + 0] = row[src + 0];
        out->rgba[dst + 1] = row[src + 1];
        out->rgba[dst + 2] = row[src + 2];
      } else if (channels == 1) {
        out->rgba[dst + 0] = row[src + 0];
        out->rgba[dst + 1] = row[src + 0];
        out->rgba[dst + 2] = row[src + 0];
      } else {
        out->rgba[dst + 0] = row[src + 0];
        out->rgba[dst + 1] = row[src + 1];
        out->rgba[dst + 2] = row[src + 2];
      }
      out->rgba[dst + 3] = 255;
    }
  }

  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  return true;
}

}  // namespace

bool decodeImageToRgba(const std::string& format,
                       const unsigned char* data,
                       size_t len,
                       RgbaImage* out) {
  if (!out) return false;
  std::string lower = format;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  if (lower == "png") return decodePngToRgba(data, len, out);
  if (lower == "jpg" || lower == "jpeg") return decodeJpegToRgba(data, len, out);

  LOG(ERROR) << "Unsupported mbtiles format: " << format;
  return false;
}

bool encodeRgbaToPng(const RgbaImage& image, std::vector<unsigned char>* out) {
  if (!out || image.width <= 0 || image.height <= 0 || image.rgba.empty())
    return false;

  png_structp png_ptr =
      png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
  if (!png_ptr) return false;
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_write_struct(&png_ptr, nullptr);
    return false;
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    return false;
  }

  PngWriteState state;
  state.out = out;
  png_set_write_fn(png_ptr, &state, pngWriteCallback, pngFlushCallback);

  png_set_IHDR(png_ptr, info_ptr, image.width, image.height, 8,
               PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png_ptr, info_ptr);

  std::vector<png_bytep> rows(image.height);
  for (int y = 0; y < image.height; ++y) {
    rows[y] = const_cast<unsigned char*>(
        image.rgba.data() + y * image.width * 4);
  }

  png_write_image(png_ptr, rows.data());
  png_write_end(png_ptr, nullptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  return true;
}

}  // namespace output
}  // namespace transitmapper

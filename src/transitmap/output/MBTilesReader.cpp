// Copyright 2026
// MBTiles reader

#include "transitmap/output/MBTilesReader.h"

#include <cstdlib>
#include <cstring>
#include <cctype>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "util/log/Log.h"

using util::ERROR;
using util::INFO;

namespace transitmapper {
namespace output {

MBTilesReader::MBTilesReader(const std::string& path) {
  if (sqlite3_open_v2(path.c_str(), &_db, SQLITE_OPEN_READONLY, nullptr) !=
      SQLITE_OK) {
    const char* err = _db ? sqlite3_errmsg(_db) : "unknown error";
    LOG(ERROR) << "Failed to open mbtiles " << path << ": " << err;
    if (_db) sqlite3_close(_db);
    _db = nullptr;
    return;
  }

  if (!readMetadata()) {
    LOG(ERROR) << "Failed to read mbtiles metadata";
  }

  const char* stmt =
      "SELECT tile_data FROM tiles WHERE zoom_level=?1 AND tile_column=?2 "
      "AND tile_row=?3";
  if (sqlite3_prepare_v2(_db, stmt, -1, &_tileStmt, nullptr) != SQLITE_OK) {
    LOG(ERROR) << "Failed to prepare tile query: " << sqlite3_errmsg(_db);
  }
}

MBTilesReader::~MBTilesReader() {
  if (_tileStmt) sqlite3_finalize(_tileStmt);
  if (_db) sqlite3_close(_db);
}

bool MBTilesReader::isOpen() const { return _db != nullptr; }

const MBTilesMetadata& MBTilesReader::metadata() const { return _meta; }

std::string MBTilesReader::readMetadataValue(const std::string& key) const {
  if (!_db) return "";
  sqlite3_stmt* stmt = nullptr;
  const char* sql = "SELECT value FROM metadata WHERE name=?1";
  if (sqlite3_prepare_v2(_db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
    return "";
  }

  sqlite3_bind_text(stmt, 1, key.c_str(), -1, SQLITE_TRANSIENT);
  std::string value;
  if (sqlite3_step(stmt) == SQLITE_ROW) {
    const unsigned char* text = sqlite3_column_text(stmt, 0);
    if (text) value = reinterpret_cast<const char*>(text);
  }
  sqlite3_finalize(stmt);
  return value;
}

bool MBTilesReader::readMetadata() {
  if (!_db) return false;
  std::string minZoom = readMetadataValue("minzoom");
  std::string maxZoom = readMetadataValue("maxzoom");
  std::string format = readMetadataValue("format");
  std::string scheme = readMetadataValue("scheme");

  if (!minZoom.empty()) _meta.minZoom = atoi(minZoom.c_str());
  if (!maxZoom.empty()) _meta.maxZoom = atoi(maxZoom.c_str());
  if (!format.empty()) _meta.format = format;
  if (!scheme.empty()) {
    std::string lower = scheme;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    _meta.scheme = lower;
  }

  LOGTO(INFO, std::cerr) << "MBTiles metadata: minzoom=" << _meta.minZoom
                         << " maxzoom=" << _meta.maxZoom
                         << " format=" << _meta.format
                         << " scheme=" << _meta.scheme;
  return true;
}

bool MBTilesReader::getTile(int z, int x, int y,
                            std::vector<unsigned char>* out) const {
  if (!_db || !_tileStmt) return false;
  if (!out) return false;

  int tileY = y;
  if (_meta.scheme != "xyz") {
    int n = 1 << z;
    tileY = (n - 1) - y;
  }

  sqlite3_reset(_tileStmt);
  sqlite3_clear_bindings(_tileStmt);
  sqlite3_bind_int(_tileStmt, 1, z);
  sqlite3_bind_int(_tileStmt, 2, x);
  sqlite3_bind_int(_tileStmt, 3, tileY);

  int rc = sqlite3_step(_tileStmt);
  if (rc != SQLITE_ROW) return false;

  const void* blob = sqlite3_column_blob(_tileStmt, 0);
  int bytes = sqlite3_column_bytes(_tileStmt, 0);
  if (!blob || bytes <= 0) return false;

  out->resize(bytes);
  std::memcpy(out->data(), blob, bytes);
  return true;
}

}  // namespace output
}  // namespace transitmapper

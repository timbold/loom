// Copyright 2026
// MBTiles reader

#ifndef TRANSITMAP_OUTPUT_MBTILESREADER_H_
#define TRANSITMAP_OUTPUT_MBTILESREADER_H_

#include <sqlite3.h>

#include <string>
#include <vector>

namespace transitmapper {
namespace output {

struct MBTilesMetadata {
  int minZoom = 0;
  int maxZoom = 22;
  std::string format = "png";
  std::string scheme = "tms";
};

class MBTilesReader {
 public:
  explicit MBTilesReader(const std::string& path);
  ~MBTilesReader();

  bool isOpen() const;
  const MBTilesMetadata& metadata() const;

  bool getTile(int z, int x, int y, std::vector<unsigned char>* out) const;

 private:
  bool readMetadata();
  std::string readMetadataValue(const std::string& key) const;

  sqlite3* _db = nullptr;
  sqlite3_stmt* _tileStmt = nullptr;
  MBTilesMetadata _meta;
};

}  // namespace output
}  // namespace transitmapper

#endif  // TRANSITMAP_OUTPUT_MBTILESREADER_H_

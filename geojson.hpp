#ifndef _geojson_hpp_
#define _geojson_hpp_

#include "geom.hpp"
#include "compactengine.hpp"

namespace complib {

  GeoCollection ReadGeoJSON(std::string filename);
  GeoCollection ReadGeoJSONFile(std::string filename);

  std::string OutScoreJSON(const GeoCollection &gc, const std::string id);
}

#endif
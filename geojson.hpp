#ifndef _geojson_hpp_
#define _geojson_hpp_

#include "geom.hpp"

namespace complib {

  std::string PrepGeoJSON(std::string geojson);

  GeoCollection ReadGeoJSON(std::string geojson);
  GeoCollection ReadGeoJSONFile(std::string filename);

  std::string OutScoreJSON(const GeoCollection &gc, const std::string id);
}

#endif
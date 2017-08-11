#ifndef _geojson_hpp_
#define _geojson_hpp_

#include "geom.hpp"

namespace complib {

  GeoCollection ReadGeoJSON(std::string filename);
  GeoCollection ReadGeoJSONFile(std::string filename);

}

#endif
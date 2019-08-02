#ifndef _wkt_hpp_
#define _wkt_hpp_

#include "geom.hpp"
#include <string>

namespace complib {

  GeoCollection ReadWKT(std::string wktstr);
  GeoCollection ReadWKTFile(std::string filename);
  std::string   GetWKT(const MultiPolygon &mp);
}


#endif

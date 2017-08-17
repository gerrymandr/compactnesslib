#ifndef _shapefill_hpp_
#define _shapefill_hpp_

#include <string>
#include "geom.hpp"

namespace complib {
  GeoCollection ReadShapefile(std::string filename);
  void WriteShapefile(const GeoCollection &mps, std::string filename, std::string layername);
}

#endif
#ifndef _shapefill_hpp_
#define _shapefill_hpp_

#include <string>
#include "geom.hpp"

namespace complib {
  GeoCollection ReadShapefile(std::string filename);
  void WriteShapefile(const GeoCollection &gc, const std::string filename);
  void WriteShapeScores(const GeoCollection &gc, const std::string filename);
}

#endif
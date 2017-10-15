#ifndef _neighbours_hpp_
#define _neighbours_hpp_

#include "geom.hpp"

namespace complib {
  void FindNeighbouringDistricts(GeoCollection &gc);
  void CalcParentOverlap(GeoCollection &subunits, GeoCollection &superunits);
}

#endif

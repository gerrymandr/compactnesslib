#ifndef _neighbours_hpp_
#define _neighbours_hpp_

#include "geom.hpp"

namespace complib {
  void FindNeighbouringDistricts(GeoCollection &gc);
  void FindExteriorDistricts(GeoCollection &subunits, const GeoCollection &superunits);
  void CalcParentOverlap(GeoCollection &subunits, GeoCollection &superunits);
}

#endif

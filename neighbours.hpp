#ifndef _neighbours_hpp_
#define _neighbours_hpp_

#include "geom.hpp"

namespace complib {
  void FindNeighbouringDistricts(GeoCollection &gc);

  void FindExteriorDistricts(
    GeoCollection &subunits,
    const GeoCollection &superunits,
    const int shrink,                //Amount by which to shrink superunits when trying to determine their children
    double border_dist_cutoff        //Subunits with at least one point within this distance of a superunit boundary are said to be part of the exterior set of that superunit
  );

  void CalcParentOverlap(GeoCollection &subunits, GeoCollection &superunits);
}

#endif

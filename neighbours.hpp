#ifndef _neighbours_hpp_
#define _neighbours_hpp_

#include "geom.hpp"

namespace complib {
  void FindNeighbouringDistricts(GeoCollection &gc);

  void FindExteriorDistricts(
    GeoCollection &subunits,
    GeoCollection &superunits,
    const double max_boundary_pt_dist, //Maximum distance between boundary points
    const double border_dist_cutoff    //Points within at least this distance of a superunit boundary are potentially included in that superunit
  );

  void CalcParentOverlap(GeoCollection &subunits, GeoCollection &superunits);
}

#endif

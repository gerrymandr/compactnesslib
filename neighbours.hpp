#ifndef _neighbours_hpp_
#define _neighbours_hpp_

#include "geom.hpp"

namespace complib {
  void FindNeighbouringDistricts(
    GeoCollection &gc,  
    const double max_neighbour_dist,     ///< Distance within which a units are considered to be neighbours.
    const double expand_bb_by            ///< Distance by which units' bounding boxes are expanded. Only districts with overlapping boxes are checked for neighbourness. Value should be >0.
  );

  void CalcParentOverlap(
    GeoCollection &subunits,
    GeoCollection &superunits,
    const double complete_inclusion_thresh, ///< A subunit with more fractional area than this in the parent are 100% included, all other potential parents are ignored
    const double not_included_thresh,       ///< A subunit with less fractional area than this in a parent disregards that parent
    const double edge_adjacency_dist,       ///< Distance within which a subunit is considered to be on the border of a superunit.  
    const bool   print_parent_columns       ///< If True, a column is printed for each possible parent indicating subunit inclusion fractions. Generally you should want this to be False.
  );

  //Finds external children who still have 100% overlap with their parent
  void FindExternalChildren(
    GeoCollection &subunits,
    GeoCollection &superunits,
    const double edge_adjacency_dist        ///< Distance within which a subunit is considered to be on the border of a superunit.
  );

}

#endif

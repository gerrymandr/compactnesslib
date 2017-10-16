#ifndef _SpIndex_hpp_
#define _SpIndex_hpp_

#include "geom.hpp"
#include "lib/spatialindex/SpatialIndex.h"
#include <memory>

namespace complib {

class SpIndexData;

typedef std::vector< std::pair<unsigned int, BoundingBox> > idbb;
  
class SpIndex {
 private:
  std::shared_ptr<SpIndexData> d;
  idbb boxes_to_insert;

  static constexpr double prop_fill_factor    = 0.7;
  static constexpr int    prop_index_capacity = 100;
  static constexpr int    prop_leaf_capacity  = 100;
  static constexpr int    prop_dimension      = 2;
  static constexpr auto   prop_tree_variant   = SpatialIndex::RTree::RV_RSTAR;

 public:
    /* creation of spatial index */

  SpIndex();
    explicit SpIndex ( const idbb &fi );

    void buildIndex();

    SpIndex& operator=( const SpIndex &other );

    /* operations */

    void insert( const unsigned int id, const BoundingBox &bb );
    void insertDeferred( const unsigned int id, const BoundingBox &bb );

    /* queries */

    std::vector<unsigned int> query( const BoundingBox &bb ) const;
    std::vector<unsigned int> query( const MultiPolygon &mp ) const;


//    friend class SpIndexData; // for access to featureInfo()
};

void AddToSpIndex(const MultiPolygon &mp, SpIndex &sp, const unsigned int id, const double expandby);

}

#endif


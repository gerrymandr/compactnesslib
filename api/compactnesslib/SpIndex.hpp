#ifndef _SpIndex_hpp_
#define _SpIndex_hpp_

#include "geom.hpp"
#include <memory>

#include <spatialindex/SpatialIndex.h>

namespace complib {

class SpIndexData;

typedef std::vector< std::pair<SpatialIndex::id_type, BoundingBox> > idbb;
  
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

    void insert( const SpatialIndex::id_type id, const BoundingBox &bb );
    void insertDeferred( const SpatialIndex::id_type id, const BoundingBox &bb );

    /* queries */

    std::vector<SpatialIndex::id_type> query( const BoundingBox &bb  ) const;
    std::vector<SpatialIndex::id_type> query( const MultiPolygon &mp ) const;
    std::vector<SpatialIndex::id_type> query( const Point2D &pt      ) const;


//    friend class SpIndexData; // for access to featureInfo()
};

void AddToSpIndex(const MultiPolygon &mp, SpIndex &sp, const SpatialIndex::id_type id, const double expandby);

}

#endif


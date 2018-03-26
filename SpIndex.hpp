#ifndef _SpIndex_hpp_
#define _SpIndex_hpp_

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "geom.hpp"

namespace complib {

template<class CoordinateType, class ValueType>
class SpIndex {
 private:
  typedef boost::geometry::model::point<CoordinateType, 2, boost::geometry::cs::cartesian> point;
  typedef boost::geometry::model::box<point> box;
  typedef std::pair<box, ValueType> value;
  typedef boost::geometry::index::rtree< value, boost::geometry::index::rstar<16> > rtree_t;
  std::vector<value> boxes_to_insert;
  rtree_t rtree;

 public:
  void addBox(const CoordinateType xmin, const CoordinateType ymin, const CoordinateType xmax, const CoordinateType ymax, const ValueType id);
  void addBoxDeferred(const CoordinateType xmin, const CoordinateType ymin, const CoordinateType xmax, const CoordinateType ymax, const ValueType id);
  void insert(const ValueType val, const BoundingBox& bb);
  std::vector<ValueType> query(const Point2D &xy) const;
  std::vector<ValueType> query(const MultiPolygon &bb) const;
  void buildIndex();
};

template<class CoordinateType, class ValueType>
void SpIndex<CoordinateType, ValueType>::addBox(
  const CoordinateType xmin, 
  const CoordinateType ymin, 
  const CoordinateType xmax, 
  const CoordinateType ymax, 
  const ValueType id
){
  box b(point(xmin,ymin),point(xmax,ymax));
  rtree.insert(std::make_pair(b, id));
}

template<class CoordinateType, class ValueType>
void SpIndex<CoordinateType, ValueType>::addBoxDeferred(
  const CoordinateType xmin, 
  const CoordinateType ymin, 
  const CoordinateType xmax, 
  const CoordinateType ymax, 
  const ValueType id
){
  box b(point(xmin,ymin),point(xmax,ymax));
  boxes_to_insert.push_back(std::make_pair(b, id));
}

template<class CoordinateType, class ValueType>
void SpIndex<CoordinateType, ValueType>::insert(const ValueType val, const BoundingBox& bb){
  addBox(bb.xmin(),bb.ymin(),bb.xmax(),bb.ymax(),val);
}

template<class CoordinateType, class ValueType>
void SpIndex<CoordinateType, ValueType>::buildIndex(){
  rtree = rtree_t(boxes_to_insert);
  boxes_to_insert.clear();
  boxes_to_insert.shrink_to_fit();
}

template<class CoordinateType, class ValueType>
std::vector<ValueType> SpIndex<CoordinateType, ValueType>::query(const Point2D &xy) const {
  box query_box(point(xy.x, xy.y), point(xy.x, xy.y));
  std::vector<value> result_s;
  rtree.query(
    boost::geometry::index::intersects(query_box),
    std::back_inserter(result_s)
  );

  std::vector<ValueType> ret;
  for(const auto &x: result_s)
    ret.emplace_back(x.second);

  return ret;
}

template<class CoordinateType, class ValueType>
std::vector<ValueType> SpIndex<CoordinateType, ValueType>::query(const MultiPolygon &mp) const {
  const auto bb = mp.bbox();
  box query_box(point(bb.xmin(), bb.ymin()), point(bb.xmax(), bb.ymax()));
  std::vector<value> result_s;
  rtree.query(
    boost::geometry::index::intersects(query_box),
    std::back_inserter(result_s)
  );

  std::vector<ValueType> ret;
  for(const auto &x: result_s)
    ret.emplace_back(x.second);

  return ret;
}


template<class CoordinateType, class ValueType>
void AddToSpIndex(const MultiPolygon &mp, SpIndex<CoordinateType, ValueType> &sp, const int id, const CoordinateType expand){
  const auto bb = mp.bbox();
  sp.addBoxDeferred(bb.xmin()-expand, bb.ymin()-expand, bb.xmax()+expand, bb.ymax()+expand, id);
}

}



#endif


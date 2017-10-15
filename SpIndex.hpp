#ifndef _SpIndex_hpp_
#define _SpIndex_hpp_

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "Point.hpp"
#include "Polygon.hpp"
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
  void addBox(const int xmin, const int ymin, const int xmax, const int ymax, const int id);
  void addBoxDeferred(const int xmin, const int ymin, const int xmax, const int ymax, const int id);
  int  queryPoint(const Point2D &xy) const;
  int  queryBox(const BoundingBox &bb) const;
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
void SpIndex<CoordinateType, ValueType>::buildIndex(){
  rtree = rtree_t(boxes_to_insert);
  boxes_to_insert.clear();
  boxes_to_insert.shrink_to_fit();
}

template<class CoordinateType, class ValueType>
std::vector<ValueType> SpIndex<CoordinateType, ValueType>::query(const Point2D &xy) const {
  box query_box(point(xy.x, xy.y), point(xy.x, xy.y));
  std::vector<ValueType> result_s;
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
  box query_box(point(bb.minx(), bb.miny()), point(bb.maxx(), bb.maxy()));
  std::vector<ValueType> result_s;
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
  sp.addBoxDeferred(bb.minx()-expand, bb.miny()-expand, bb.maxx()+expand, bb.maxy()+expand, id);
}

}



#endif


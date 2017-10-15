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
  int  query(const Point2D &xy) const;
  int  query(const BoundingBox &bb) const;
  void buildIndex();
};

void SpIndex::addBox(
  const int xmin, 
  const int ymin, 
  const int xmax, 
  const int ymax, 
  const int id
){
  box b(point(xmin,ymin),point(xmax,ymax));
  rtree.insert(std::make_pair(b, id));
}

void SpIndex::addBoxDeferred(
  const int xmin, 
  const int ymin, 
  const int xmax, 
  const int ymax, 
  const int id
){
  box b(point(xmin,ymin),point(xmax,ymax));
  boxes_to_insert.push_back(std::make_pair(b, id));
}

void SpIndex::buildIndex(){
  rtree = rtree_t(boxes_to_insert);
  boxes_to_insert.clear();
  boxes_to_insert.shrink_to_fit();
}

std::vector<ValueType> SpIndex::query(const Point2D &xy) const {
  box query_box(point(xy.x, xy.y), point(xy.x, xy.y));
  std::vector<ValueType> result_s;
  rtree.query(
    boost::geometry::index::intersects(query_box),
    std::back_inserter(result_s)
  );
  return result_s;
}

std::vector<ValueType> SpIndex::query(const MultiPolygon &mp) const {
  const auto bb = mp.bbox();
  box query_box(point(bb.xmin(), bb.ymin()), point(bb.xmax(), bb.ymax()));
  std::vector<ValueType> result_s;
  rtree.query(
    boost::geometry::index::intersects(query_box),
    std::back_inserter(result_s)
  );
  return result_s;
}


template<class CoordinateType, class ValueType>
void AddToSpIndex(const MultiPolygon &mp, SpIndex &sp, const int id){
  const auto bb = mp.bbox();
  sp.addBoxDeferred(bb.xmin(),bb.ymin(),bb.xmax(),bb.ymax(),id);
}

}



#endif


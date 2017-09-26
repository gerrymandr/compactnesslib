#ifndef _geom_hpp_
#define _geom_hpp_

#include <vector>
#include <cmath>
#include <set>
#include <string>
#include "Props.hpp"
#include "lib/clipper.hpp"
#include "lib/iterator_tpl.h"

namespace complib {

namespace cl = ClipperLib;

class Point2D;
class Polygon;
class Ring;
class MultiPolygon;
class Geometry;

typedef std::vector<Point2D>      Points;
typedef std::vector<Polygon>      Polygons;
typedef std::vector<Ring>         Rings;
typedef std::vector<MultiPolygon> MultiPolygons;

void PrintProps(const Props &ps);

class Point2D {
 public:
  double x;
  double y;
  Point2D() = default;
  Point2D(double x0, double y0);
};

class Ring {
 public:
  Points v;
  Ring() = default;
  Ring(const std::vector<Point2D> &ptvec);
  mutable std::vector<Point2D> hull;
  Ring getHull() const;
  mutable ClipperLib::Path clipper_paths;
  EXPOSE_STL_VECTOR(v);
};

class Polygon {
 public:
  Rings v;
  EXPOSE_STL_VECTOR(v);
};

class MultiPolygon {
 public:
  Polygons v;
  Props props;
  Scores scores;
  mutable Ring hull;
  const Ring& getHull() const;
  void toRadians();
  void toDegrees();
  MultiPolygon intersect(const MultiPolygon &b) const;
  mutable ClipperLib::Paths clipper_paths;
  void reverse();
  EXPOSE_STL_VECTOR(v);

  std::set<unsigned int> neighbours;
};

class GeoCollection {
 public:
  MultiPolygons v;
  std::string prj_str;
  void reverse();
  void correctWindingDirection();
  EXPOSE_STL_VECTOR(v);
};



double EuclideanDistance(const Point2D &a, const Point2D &b);



double area(const Ring &r);
double areaIncludingHoles(const Polygon &p);
double areaIncludingHoles(const MultiPolygon &mp);
double areaExcludingHoles(const MultiPolygon &mp);

double areaHoles(const Polygon &p);
double areaHoles(const MultiPolygon &mp);

double perim(const Ring &r);
double perimExcludingHoles(const Polygon &p);
double perimExcludingHoles(const MultiPolygon &mp);
double perimIncludingHoles(const Polygon &p);
double perimIncludingHoles(const MultiPolygon &mp);

double perimHoles(const Polygon &p);
double perimHoles(const MultiPolygon &mp);

double hullArea(const Ring &r);
double hullAreaOfOuter(const Polygon &p);
double hullAreaPolygonOuterRings(const MultiPolygon &mp);

double hullAreaOfHoles(const Polygon &p);
double hullAreaOfHoles(const MultiPolygon &mp);

double diameter(const Ring &r);
double diameterOuter(const Polygon &p);
double diameterOfEntireMultiPolygon(const MultiPolygon &mp);

unsigned holeCount(const Polygon &p);
unsigned polyCount(const MultiPolygon &mp);
unsigned holeCount(const MultiPolygon &mp);


const cl::Path& ConvertToClipper(const Ring &ring, const bool reversed);
const cl::Paths& ConvertToClipper(const MultiPolygon &mp, const bool reversed);

template<class T, class U>
double IntersectionArea(const T &a, const U &b) {
  const auto paths_a = ConvertToClipper(a,false);
  const auto paths_b = ConvertToClipper(b,false);

  cl::Clipper clpr;
  clpr.AddPaths(paths_a, cl::ptSubject, true);
  clpr.AddPaths(paths_b, cl::ptClip, true);
  cl::Paths solution;
  clpr.Execute(cl::ctIntersection, solution, cl::pftEvenOdd, cl::pftEvenOdd);

  double area = 0;
  for(const auto &path: solution)
    area += cl::Area(path);
  return area;
}



template<class T>
std::pair<Point2D, Point2D> MostDistantPoints(const T &geom){
  //We'll use the Convex Hull to find the two most distant points
  const auto &hull = geom.getHull();

  std::pair<unsigned int, unsigned int> idx_maxpts;
  double maxdist = 0;

  //TODO: There's a faster way to do this  
  for(unsigned int i=0;i<hull.size();i++)
  for(unsigned int j=i+1;j<hull.size();j++){
    const double dist = EuclideanDistance(hull.at(i),hull.at(j));
    if(dist>maxdist){
      idx_maxpts = std::make_pair(i,j);
      maxdist    = dist;
    }
  }

  return std::make_pair(hull.at(idx_maxpts.first), hull.at(idx_maxpts.second));
}

MultiPolygon GetBoundingCircle(const MultiPolygon &mp);
MultiPolygon GetBoundingCircleMostDistant(const MultiPolygon &mp);

}

#endif

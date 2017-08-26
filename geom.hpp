#ifndef _geom_hpp_
#define _geom_hpp_

#include <vector>
#include <cmath>
#include <string>
#include <memory>
#include "Props.hpp"
#include "clipper.h"

namespace complib {

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

class Ring : public Points {
 public:
  Ring() = default;
  Ring(std::vector<Point2D>::iterator first, std::vector<Point2D>::iterator last);
  mutable std::unique_ptr<Ring> hull;
  void getHull() const;
};

class Polygon : public Rings {

};

class MultiPolygon : public Polygons {
 public:
  Props props;
  Scores scores;
  void toRadians();
  void toDegrees();
  MultiPolygon intersect(const MultiPolygon &b) const;
  mutable ClipperLib::Paths clipper_paths;
};

class GeoCollection : public MultiPolygons {
 public:
  std::string prj_str;
};



double area(const Ring &r);
double areaOuter(const Polygon &p);
double areaOfPolygonsIncludingHoles(const MultiPolygon &mp);

double areaHoles(const Polygon &p);
double areaHoles(const MultiPolygon &mp);

double perim(const Ring &r);
double perimOuter(const Polygon &p);
double perimPolygonOuterRings(const MultiPolygon &mp);

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

double IntersectionArea(const MultiPolygon &a, const MultiPolygon &b);

}

#endif

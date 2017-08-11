#ifndef _geom_hpp_
#define _geom_hpp_

#include <vector>
#include <cmath>
#include <memory>
#include <map>
#include <string>
#include <any>

class Point2D;
class Polygon;
class Ring;
class MultiPolygon;
class Geometry;

typedef std::map<std::string,std::any> Props;

typedef std::pair<bool, double>   CachedValue;

typedef std::vector<std::shared_ptr<Geometry> > Geometries;
typedef std::vector<Point2D>      Points;
typedef std::vector<Polygon>      Polygons;
typedef std::vector<Ring>         Rings;
typedef std::vector<MultiPolygon> MultiPolygons;

static const double DEG_TO_RAD = M_PI/180.0;
static const double RAD_TO_DEG = 180.0/M_PI;

void PrintProps(const Props &ps);

class Geometry {
 protected:
  mutable bool valid = false;
 public:
  virtual double minX() const = 0;
  virtual double maxX() const = 0;
  virtual double minY() const = 0;
  virtual double maxY() const = 0;
  virtual double avgX() const = 0;
  virtual double avgY() const = 0;
  virtual double sumX() const = 0;
  virtual double sumY() const = 0;
  virtual unsigned points() const = 0;
  virtual double area () const = 0;
  virtual double perim() const = 0;
  virtual void toRadians() = 0;
  virtual void toDegrees() = 0;
  virtual bool containsPoint(const Point2D &xy) const = 0;
  virtual double hullArea() const = 0;
  virtual double diameter() const = 0;
  virtual void print() const = 0;
  void invalidate();
};

class Point2D : public Geometry {
 public:
  double x;
  double y;
  Point2D() = default;
  Point2D(double x0, double y0);
  double minX() const override;
  double maxX() const override;
  double minY() const override;
  double maxY() const override;
  double avgX() const override;
  double avgY() const override;
  double sumX() const override;
  double sumY() const override;
  unsigned points() const override;
  void toRadians() override;
  void toDegrees() override;
  double area() const override;
  double perim() const override;
  bool containsPoint(const Point2D &xy) const override;
  double hullArea() const override;
  double diameter() const override;
  void print() const override;
};

class Ring : public Geometry, public Points {
 public:
  double minX() const override;
  double maxX() const override;
  double minY() const override;
  double maxY() const override;
  double avgX() const override;
  double avgY() const override;
  double sumX() const override;
  double sumY() const override;
  unsigned points() const override;
  double area() const override;
  double perim() const override;
  void toRadians() override;
  void toDegrees() override;
  bool containsPoint(const Point2D &xy) const override;
  Ring getHull() const;
  double hullArea() const override;
  double diameter() const override;
  void print() const override;
};

class Polygon : public Geometry {
 public:
  Ring outer;
  Rings holes;
  double minX() const override;
  double maxX() const override;
  double minY() const override;
  double maxY() const override;
  double avgX() const override;
  double avgY() const override;
  double sumX() const override;
  double sumY() const override;
  unsigned points() const override;
  double area() const override;
  double perim() const override;
  void toRadians() override;
  void toDegrees() override;
  bool containsPoint(const Point2D &xy) const override;
  Ring getHull() const;
  double hullArea() const override;
  double diameter() const override;
  void print() const override;
};

class MultiPolygon : public Geometry, public Polygons {
 public:
  Props props;
  double minX() const override;
  double maxX() const override;
  double minY() const override;
  double maxY() const override;
  double avgX() const override;
  double avgY() const override;
  double sumX() const override;
  double sumY() const override;
  unsigned points() const override;
  double area() const override;
  double perim() const override;
  void toRadians() override;
  void toDegrees() override;
  bool containsPoint(const Point2D &xy) const override;
  double hullArea() const override;
  double diameter() const override;
  void print() const override;
};

class GeoCollection : public MultiPolygons {
 public:
  std::string srs;
};

#endif

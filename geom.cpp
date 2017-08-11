#include "geom.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include "doctest.h"


inline double EuclideanDistance(const Point2D &a, const Point2D &b){
  const double xd = (a.x-b.x);
  const double yd = (a.y-b.y);
  return std::sqrt(xd*xd+yd*yd);
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int FindOrientation(
  const Point2D &p,
  const Point2D &q,
  const Point2D &r
){
  const int val = (q.y-p.y) * (r.x-q.x) - (q.x-p.x) * (r.y-q.y);

  if (val == 0) 
    return 0;   //Colinear
  else if (val>0)
    return 1;   //Clockwise
  else
    return 2;   //Counter-clockwise
}








void Geometry::invalidate(){
  valid = false;
}


Point2D::Point2D(double x0, double y0) {
  x = x0;
  y = y0;
}

void Point2D::toRadians() {
  invalidate();

  assert(-180<=x && x<=180);
  assert(-90<=y && y<=90);

  x *= DEG_TO_RAD;
  y *= DEG_TO_RAD;

  assert(-M_PI<=x && x<=M_PI);
  assert(-M_PI/2<=y && y<=M_PI/2);
}

void Point2D::toDegrees() {
  invalidate();

  assert(-M_PI<=x && x<=M_PI);
  assert(-M_PI/2<=y && y<=M_PI/2);

  x *= RAD_TO_DEG;
  y *= RAD_TO_DEG;

  assert(-180<=x && x<=180);
  assert(-90<=y && y<=90);
}

double Point2D::minX()  const { return x; }
double Point2D::maxX()  const { return x; }
double Point2D::minY()  const { return y; }
double Point2D::maxY()  const { return y; }
double Point2D::area()  const { return 0; }
double Point2D::perim() const { return 0; }
bool   Point2D::containsPoint(const Point2D &xy) const { return x==xy.x && y==xy.y; }
double Point2D::hullArea() const { return 0; }
double Point2D::diameter() const { return 0; }

void Point2D::print() const {
  std::cout<<x<<" "<<y<<"\n";
}

double Ring::minX() const {
  static double minx = std::numeric_limits<double>::infinity();
  if(valid)
    return minx;

  for(const auto &p: *this)
    minx = std::min(p.x,minx);

  return minx;
}

double Ring::maxX() const {
  static double maxx = std::numeric_limits<double>::infinity();
  if(valid)
    return maxx;

  for(const auto &p: *this)
    maxx = std::min(p.x,maxx);
  
  return maxx;
}

double Ring::minY() const {
  static double miny = std::numeric_limits<double>::infinity();
  if(valid)
    return miny;

  for(const auto &p: *this)
    miny = std::min(p.y,miny);
  
  return miny;
}

double Ring::maxY() const {
  static double maxy = std::numeric_limits<double>::infinity();
  if(valid)
    return maxy;

  for(const auto &p: *this)
    maxy = std::min(p.y,maxy);
  
  return maxy;
}

double Ring::area() const {
  static double area = 0;
  if(valid)
    return area;

  const auto &self = *this;

  //The "shoelace" algorithm
  unsigned int j = size()-1;
  for(unsigned int i=0;i<size();i++){
    area += (self[j].x + self[i].x) * (self[j].y - self[i].y);
    j = i;
  }
  
  area = std::abs(area/2);

  return area;
}

double Ring::perim() const {
  static double perim = 0;
  if(valid)
    return perim;

  if(size()==1)
    return 0;

  const auto &self = *this;

  for(unsigned int i=0;i<size()-1;i++)
    perim += EuclideanDistance(self[i],self[i+1]);
  perim += EuclideanDistance(self.front(),self.back());

  return perim;
}

void   Ring::toRadians() {
  invalidate();
  for(auto &p: *this)
    p.toRadians();
}

void   Ring::toDegrees() {
  invalidate();
  for(auto &p: *this)
    p.toDegrees();
}

bool   Ring::containsPoint(const Point2D &xy) const {
  unsigned int i, j;
  int c = 0;
  auto &self = *this;
  for (i = 0, j = self.size()-1; i < self.size(); j = i++) {
    if ( ((self[i].y>xy.y) != (self[j].y>xy.y)) &&
     (xy.x < (self[j].x-self[i].x) * (xy.y-self[i].y) / (self[j].y-self[i].y) + self[i].x) )
       c = !c;
  }
  return c;
}


//Andrew's monotone chain convex hull algorithm
//https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
Ring Ring::getHull() const {
  if (size() < 3)
    throw std::runtime_error("There must be at least 3 points for a convex hull!");

  const auto cross = [&](const Point2D &a, const Point2D &b, const Point2D &o) {
    return (a.x-o.x)*(b.y-o.y) - (a.y-o.y)*(b.x-o.x);
  };

  const auto &self = *this;

  std::vector<unsigned int> idx(size(),0);
  for(unsigned int i=0;i<size();i++)
    idx[i] = i;

  std::sort(idx.begin(),idx.end(),[&](const unsigned int a, const unsigned int b){
    if(self[a].x==self[b].x)
      return self[a].y<self[b].y;
    return self[a].x<self[b].x;
  });

  //Lower half of hull
  Ring L;
  for(unsigned int i=0;i<size();i++){
    const auto &thisp = self[idx[i]];
    while(L.size()>=2 && cross(L[L.size()-2], L[L.size()-1], thisp)<=0)
      L.pop_back();
    L.push_back(thisp);
  }

  Ring U;
  for(int i=((signed int)size())-1;i>=0;i--){
    const auto &thisp = self[idx[i]];
    while(U.size()>=2 && cross(U[U.size()-2], U[U.size()-1], thisp)<=0)
      U.pop_back();
    U.push_back(thisp);
  }

  L.pop_back(); //Last point of L is first point of U
  U.pop_back(); //Last point of U is first point of L

  L.insert(L.end(),U.begin(),U.end());

  return L;
}

double Ring::hullArea() const {
  const auto hull = getHull();
  return hull.area();
}

double Ring::diameter() const {
  static double maxdist = 0;
  if(valid) return maxdist;

  const auto hull = getHull();

  for(unsigned int i=0;i<hull.size();i++)
  for(unsigned int j=i+1;j<hull.size();j++)
    maxdist = std::max(maxdist,EuclideanDistance(hull[i],hull[j]));
  return maxdist;
}

void Ring::print() const {
  for(const auto &p: *this)
    p.print();
}






double Polygon::minX() const {
  static double minx=std::numeric_limits<double>::infinity();
  if(valid) return minx;

  for(const auto &p: outer)
    minx = std::min(p.minX(),minx);

  return minx;
}

double Polygon::maxX() const {
  static double maxx=-std::numeric_limits<double>::infinity();
  if(valid) return maxx;

  for(const auto &p: outer)
    maxx = std::max(p.maxX(),maxx);

  return maxx;
}

double Polygon::minY() const {
  static double miny=std::numeric_limits<double>::infinity();
  if(valid) return miny;

  for(const auto &p: outer)
    miny=std::min(p.minY(),miny);
  return miny;
}

double Polygon::maxY() const {
  static double maxy=-std::numeric_limits<double>::infinity();
  if(valid) return maxy;

  for(const auto &p: outer)
    maxy=std::max(p.maxY(),maxy);
  return maxy;
}

void Polygon::toRadians() {
  invalidate();
  outer.toRadians();
  for(auto &h: holes)
    h.toRadians();
}

void Polygon::toDegrees() {
  invalidate();
  outer.toDegrees();
  for(auto &h: holes)
    h.toDegrees();
}

double Polygon::area() const {
  static double area = 0;
  if(valid) return area;

  area += outer.area();
  for(auto &h: holes)
    area -= h.area();

  return area;
}

double Polygon::perim() const {
  static double perim = 0;
  if(valid) return perim;

  perim += outer.perim();
  for(auto &h: holes)
    perim += h.perim();

  return perim;
}

//Info: https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
//Info: https://stackoverflow.com/questions/8721406/how-to-determine-if-a-point-is-inside-a-2d-convex-polygon/23223947#23223947
//Info: http://stackoverflow.com/a/2922778/752843
bool Polygon::containsPoint(const Point2D &xy) const {
  unsigned int i, j;
  int c = 0;
  for (i = 0, j = outer.size()-1; i < outer.size(); j = i++) {
    const auto &oi = outer[i];
    const auto &oj = outer[j];
    if ( ((oi.y>xy.y) != (oj.y>xy.y)) &&
     (xy.x < (oj.x-oi.x) * (xy.y-oi.y) / (oj.y-oi.y) + oi.x) )
       c = !c;
  }
  return c;
}

Ring Polygon::getHull() const {
  return outer.getHull();
}

double Polygon::hullArea() const {
  if(holes.size()>0)
    std::cerr<<"Warning: Taking the hull area of a polygon with holes!"<<std::endl;
  return outer.hullArea();
}

double Polygon::diameter() const {
  if(holes.size()>0)
    std::cerr<<"Warning: Taking the diameter of a polygon with holes!"<<std::endl;
  return outer.diameter();
}

void Polygon::print() const {
  outer.print();
}




double MultiPolygon::minX() const {
  static double minx = std::numeric_limits<double>::infinity();
  if(valid) return minx;

  for(const auto &p: *this)
    minx = std::min(p.minX(),minx);

  return minx;
}

double MultiPolygon::maxX() const {
  static double maxx = std::numeric_limits<double>::infinity();
  if(valid) return maxx;

  for(const auto &p: *this)
    maxx = std::min(p.maxX(),maxx);

  return maxx;
}

double MultiPolygon::minY() const {
  static double miny=std::numeric_limits<double>::infinity();
  if(valid) return miny;

  for(const auto &p: *this)
    miny = std::min(p.minY(),miny);

  return miny;
}

double MultiPolygon::maxY() const {
  static double maxy=std::numeric_limits<double>::infinity();
  if(valid) return maxy;

  for(const auto &p: *this)
    maxy = std::min(p.maxY(),maxy);

  return maxy;
}

double MultiPolygon::area() const {
  static double area=0;
  if(valid) return area;

  for(const auto &p: *this)
    area += p.area();

  return area;
}

double MultiPolygon::perim() const {
  static double perim=0;
  if(valid) return perim;

  for(const auto &p: *this)
    perim += p.perim();

  return perim;
}

void MultiPolygon::toRadians() {
  invalidate();
  for(auto &p: *this)
    p.toDegrees();
}

void MultiPolygon::toDegrees() {
  invalidate();
  for(auto &p: *this)
    p.toDegrees();
}

bool MultiPolygon::containsPoint(const Point2D &xy) const {
  for(const auto &p: *this)
    if(p.containsPoint(xy))
      return true;
  return false;
}

double MultiPolygon::hullArea() const {
  static double harea=0;
  if(valid) return harea;

  for(const auto &p: *this)
    harea += p.hullArea();

  return harea;
}

double MultiPolygon::diameter() const {
  static double diam = 0;
  if(valid) return diam;

  if(size()>1)
    std::cerr<<"Warning: Taking diameter of a multipolygon with more than one child!"<<std::endl;

  for(const auto &p: *this)
    diam += p.diameter();

  return diam;
}

void MultiPolygon::print() const {
  for(const auto &mp: *this){
    mp.print();
    std::cout<<"\n"<<std::endl;
  }
}







TEST_CASE("Point2D: Conversion to Radians"){
  Point2D p(-93, 45);
  p.toRadians();
  CHECK(p.x==doctest::Approx(-93*DEG_TO_RAD));
  CHECK(p.y==doctest::Approx(45*DEG_TO_RAD));
}



TEST_CASE("Point2D: Conversion to Degrees"){
  Point2D p(-93*DEG_TO_RAD, 45*DEG_TO_RAD);
  p.toDegrees();
  CHECK(p.x==doctest::Approx(-93));
  CHECK(p.y==doctest::Approx(45));
}











TEST_CASE("Polygon"){
  Polygon p;
  Point2D a(-93,45);
  Point2D b(-93,50);
  Point2D c(-90,50);
  Point2D d(-90,45);
  p.outer.push_back(a);
  p.outer.push_back(b);
  p.outer.push_back(c);
  p.outer.push_back(d);
  p.toRadians();
  b.toRadians();
  CHECK(p.exterior[1].x==doctest::Approx(b.x));
  CHECK(p.exterior[1].y==doctest::Approx(b.y));
  p.toDegrees();
  CHECK(p.exterior[2].x==doctest::Approx(c.x));
  CHECK(p.exterior[2].y==doctest::Approx(c.y));
  p.containsPoint(Point2D(-92,47));
  p.containsPoint(Point2D(-94,47));
}
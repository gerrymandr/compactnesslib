#include "geom.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <memory>
#include "clipper.hpp"
#include "doctest.h"

static const double DEG_TO_RAD = M_PI/180.0;
static const double RAD_TO_DEG = 180.0/M_PI;

namespace complib {


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



Point2D::Point2D(double x0, double y0) {
  x = x0;
  y = y0;
}



// bool   Ring::containsPoint(const Point2D &xy) const {
//   unsigned int i, j;
//   int c = 0;
//   auto &self = *this;
//   for (i = 0, j = self.size()-1; i < self.size(); j = i++) {
//     if ( ((self[i].y>xy.y) != (self[j].y>xy.y)) &&
//      (xy.x < (self[j].x-self[i].x) * (xy.y-self[i].y) / (self[j].y-self[i].y) + self[i].x) )
//        c = !c;
//   }
//   return c;
// }

Ring::Ring(std::vector<Point2D>::iterator first, std::vector<Point2D>::iterator last) :
  std::vector<Point2D>(first,last) {}


//Andrew's monotone chain convex hull algorithm
//https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
const Ring& Ring::getHull() const {
  if(hull!=nullptr)
    *hull;

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

  hull.reset(new Ring(L.begin(),L.end()));

  //NOTE: It may be necessary to reverse the ring to give a positive area

  return *hull;
}






void MultiPolygon::toRadians(){
  for(auto &poly: *this)
  for(auto &ring: poly)
  for(auto &pt: ring){
    pt.x *= DEG_TO_RAD;
    pt.y *= DEG_TO_RAD;
  }
}

void MultiPolygon::toDegrees(){
  for(auto &poly: *this)
  for(auto &ring: poly)
  for(auto &pt: ring){
    pt.x *= RAD_TO_DEG;
    pt.y *= RAD_TO_DEG;
  }
}

const Ring& MultiPolygon::getHull() const {
  if(!hull.empty())
    return hull;

  //Put all of the points into a ring 
  Ring temp;
  for(const auto &p: *this)
  for(const auto &r: p)
    temp.insert(temp.end(),r.begin(),r.end());
  temp.getHull();

  //Move temporary's hull into mp's hull
  std::swap(hull,*temp.hull);

  return hull;
}

void MultiPolygon::reverse() {
  for(auto &poly: *this){
    std::reverse(poly.at(0).begin(),poly.at(0).end());
    for(unsigned int i=1;i<poly.size();i++)
      std::reverse(poly.at(i).begin(),poly.at(i).end());
  }
}

void GeoCollection::reverse() {
  for(auto &mp: *this)
    mp.reverse();
}

void GeoCollection::correctWindingDirection(){
  if(areaExcludingHoles(this->at(0))<0){
    std::cerr<<"Reversed winding of polygons!"<<std::endl;
    reverse();
  }
}



double area(const Ring &r){
  double area = 0;

  //The "shoelace" algorithm
  unsigned int j = r.size()-1;
  for(unsigned int i=0;i<r.size();i++){
    area += (r[j].x + r[i].x) * (r[j].y - r[i].y);
    j = i;
  }
  
  area = area/2.;

  return std::abs(area);
}

//Produces the same answers as the foregoing, but with the opposite signedness of area
// double area2(const Ring &r){
//   double area = 0;

//   //The "shoelace" algorithm
//   unsigned int j = r.size()-1;
//   for(unsigned int i=0;i<r.size();i++){
//     area += (r[j].x * r[i].y) - (r[i].x * r[j].y);
//     j = i;
//   }
  
//   area = area/2.;

//   return area;
// }


double perim(const Ring &r){
  double perim = 0;

  if(r.size()==1)
    return 0;

  for(unsigned int i=0;i<r.size()-1;i++)
    perim += EuclideanDistance(r[i],r[i+1]);
  perim += EuclideanDistance(r.front(),r.back());

  return perim;
}

double hullArea(const Ring &r){
  r.getHull();
  return area(*r.hull);
}

double areaOuter(const Polygon &p){
  return area(p.at(0));
}

double areaHoles(const Polygon &p){
  return std::accumulate(p.begin()+1,p.end(),0.0,[](const double b, const Ring &r){ return b+area(r);});
}

double areaOfPolygonsIncludingHoles(const MultiPolygon &mp){
  return std::accumulate(mp.begin(),mp.end(),0.0,[](const double b, const Polygon &p){ return b+areaOuter(p);}); 
}

double areaHoles(const MultiPolygon &mp){
  return std::accumulate(mp.begin(),mp.end(),0.0,[](const double b, const Polygon &p){ return b+areaHoles(p);}); 
}

double perimOuter(const Polygon &p){
  return perim(p.at(0));
}

double perimHoles(const Polygon &p){
  return std::accumulate(p.begin()+1,p.end(),0.0,[](const double b, const Ring &r){ return b+perim(r);});
}

double perimPolygonOuterRings(const MultiPolygon &mp){
  return std::accumulate(mp.begin(),mp.end(),0.0,[](const double b, const Polygon &p){ return b+perimOuter(p);}); 
}

double perimHoles(const MultiPolygon &mp){
  return std::accumulate(mp.begin(),mp.end(),0.0,[](const double b, const Polygon &p){ return b+perimHoles(p);}); 
}


double hullAreaOuter(const Polygon &p){
  return hullArea(p.at(0));
}

double hullAreaPolygonOuterRings(const MultiPolygon &mp){
  return std::accumulate(mp.begin(),mp.end(),0.0,[](const double b, const Polygon &p){ return b+hullAreaOuter(p);}); 
}

double hullAreaOfHoles(const Polygon &p){
  return std::accumulate(p.begin()+1,p.end(),0.0,[](const double b, const Ring &r){ return b+hullArea(r);});
}

double hullAreaOfHoles(const MultiPolygon &mp){
  return std::accumulate(mp.begin(),mp.end(),0.0,[](const double b, const Polygon &p){ return b+hullAreaOfHoles(p);}); 
}

double diameter(const Ring &r){
  const auto &hull = r.getHull();

  //TODO: There's a faster way to do this
  double maxdist = 0;
  for(unsigned int i=0;i<hull.size();i++)
  for(unsigned int j=i+1;j<hull.size();j++)
    maxdist = std::max(maxdist,EuclideanDistance(hull.at(i),hull.at(j)));
  return maxdist;
}

double diameterOuter(const Polygon &p){
  return diameter(p.at(0));
}

double diameterOfEntireMultiPolygon(const MultiPolygon &mp){
  return diameter(mp.getHull());
}



const cl::Path& ConvertToClipper(const Ring &ring, const bool reversed){
  if(!ring.clipper_paths.empty())
    return ring.clipper_paths;

  cl::Path path;
  if(!reversed){
    for(const auto &pt: ring)
      path.emplace_back((long long)pt.x,(long long)pt.y);
  } else {
    for(auto pt=ring.rbegin();pt!=ring.rend();pt++)
      path.emplace_back((long long)pt->x,(long long)pt->y);    
  }

  std::swap(ring.clipper_paths,path);

  return ring.clipper_paths;
}


const cl::Paths& ConvertToClipper(const MultiPolygon &mp, const bool reversed) {
  if(!mp.clipper_paths.empty())
    return mp.clipper_paths;

  cl::Paths paths;

  for(const auto &poly: mp){
    //Send in outer perimter
    paths.push_back(ConvertToClipper(poly.at(0), reversed));

    //Send in the holes
    for(unsigned int i=1;i<poly.size();i++)
      paths.push_back(ConvertToClipper(poly.at(i), !reversed));
  }

  std::swap(mp.clipper_paths, paths);

  return mp.clipper_paths;
}







}

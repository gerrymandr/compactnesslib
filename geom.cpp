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
#include "lib/clipper.hpp"
#include "lib/doctest.h"
#include "lib/miniball.hpp"

static const double DEG_TO_RAD = M_PI/180.0;
static const double RAD_TO_DEG = 180.0/M_PI;

namespace complib {


double EuclideanDistance(const Point2D &a, const Point2D &b){
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

Ring::Ring(const std::vector<Point2D> &ptvec) {
  v = ptvec;
}


//Andrew's monotone chain convex hull algorithm
//https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
Ring Ring::getHull() const {
  if(!hull.empty())
    return hull;

  if (v.size() < 3)
    throw std::runtime_error("There must be at least 3 points for a convex hull!");

  const auto cross = [&](const Point2D &a, const Point2D &b, const Point2D &o) {
    return (a.x-o.x)*(b.y-o.y) - (a.y-o.y)*(b.x-o.x);
  };

  std::vector<unsigned int> idx(v.size(),0);
  for(unsigned int i=0;i<v.size();i++)
    idx[i] = i;

  std::sort(idx.begin(),idx.end(),[&](const unsigned int a, const unsigned int b){
    if(v[a].x==v[b].x)
      return v[a].y<v[b].y;
    return v[a].x<v[b].x;
  });

  //Lower half of hull
  std::vector<Point2D> L;
  for(unsigned int i=0;i<v.size();i++){
    const auto &thisp = v[idx[i]];
    while(L.size()>=2 && cross(L[L.size()-2], L[L.size()-1], thisp)<=0)
      L.pop_back();
    L.push_back(thisp);
  }

  std::vector<Point2D> U;
  for(int i=((signed int)v.size())-1;i>=0;i--){
    const auto &thisp = v[idx[i]];
    while(U.size()>=2 && cross(U[U.size()-2], U[U.size()-1], thisp)<=0)
      U.pop_back();
    U.push_back(thisp);
  }

  L.pop_back(); //Last point of L is first point of U
  U.pop_back(); //Last point of U is first point of L

  L.insert(L.end(),U.begin(),U.end());

  //Close the ring
  L.push_back(L.front());

  hull = L;

  //NOTE: It may be necessary to reverse the ring to give a positive area

  return hull;
}






void MultiPolygon::toRadians(){
  for(auto &poly: v)
  for(auto &ring: poly)
  for(auto &pt: ring){
    pt.x *= DEG_TO_RAD;
    pt.y *= DEG_TO_RAD;
  }
}

void MultiPolygon::toDegrees(){
  for(auto &poly: v)
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
  for(const auto &poly: v)
  for(const auto &ring: poly)
    temp.v.insert(temp.end(),ring.begin(),ring.end());
  temp.getHull();

  //Move temporary's hull into mp's hull
  std::swap(hull.v,temp.hull);

  return hull;
}

void MultiPolygon::reverse() {
  for(auto &poly: v){
    std::reverse(poly.at(0).begin(),poly.at(0).end());
    for(unsigned int i=1;i<poly.size();i++)
      std::reverse(poly.at(i).begin(),poly.at(i).end());
  }
}

void GeoCollection::reverse() {
  for(auto &mp: v)
    mp.reverse();
}

void GeoCollection::correctWindingDirection(){
  if(areaExcludingHoles(v.at(0))<0){
    std::cerr<<"Reversed winding of polygons!"<<std::endl;
    reverse();
  }
}



double area(const Ring &r){
  double area = 0;

  if(r.size()<3)
    return 0;

  //The "shoelace" algorithm
  unsigned int j = r.size()-1;
  for(unsigned int i=0;i<r.size();i++){
    //std::cerr<<"j="<<j<<", i="<<i<<", size="<<r.size();
    //std::cerr<<", xj="<<r[j].x<<", yj="<<r[j].y<<std::endl;
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
  return area(r.hull);
}

double areaIncludingHoles(const Polygon &p){
  return area(p.at(0));
}

double areaHoles(const Polygon &p){
  return std::accumulate(p.begin()+1,p.end(),0.0,[](const double b, const Ring &r){ return b+area(r);});
}

double areaIncludingHoles(const MultiPolygon &mp){
  return std::accumulate(mp.begin(),mp.end(),0.0,[](const double b, const Polygon &p){ return b+areaIncludingHoles(p);}); 
}

double areaExcludingHoles(const MultiPolygon &mp){
  return areaIncludingHoles(mp)-areaHoles(mp);
}

double areaHoles(const MultiPolygon &mp){
  return std::accumulate(mp.begin(),mp.end(),0.0,[](const double b, const Polygon &p){ return b+areaHoles(p);}); 
}

double perimExcludingHoles(const Polygon &p){
  return perim(p.at(0));
}

double perimHoles(const Polygon &p){
  return std::accumulate(p.begin()+1,p.end(),0.0,[](const double b, const Ring &r){ return b+perim(r);});
}

double perimExcludingHoles(const MultiPolygon &mp){
  return std::accumulate(mp.begin(),mp.end(),0.0,[](const double b, const Polygon &p){ return b+perimExcludingHoles(p);}); 
}

double perimIncludingHoles(const Polygon &p){
  return std::accumulate(p.begin(),p.end(),0.0,[](const double b, const Ring &r){ return b+perim(r);});
}

double perimIncludingHoles(const MultiPolygon &mp){
  return std::accumulate(mp.begin(),mp.end(),0.0,[](const double b, const Polygon &p){ return b+perimIncludingHoles(p);}); 
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



MultiPolygon GetBoundingCircle(const MultiPolygon &mp){
  //Number of unique points from which to construct the circle. The circle will
  //have one more point of than this in order to form a closed ring).
  const int CIRCLE_PT_COUNT = 1000;

  std::vector< std::vector<double> > pts;
  for(const auto &poly: mp)
  for(const auto &ring: poly)
  for(const auto &pt: ring)
    pts.push_back(std::vector<double>({{pt.x,pt.y}}));

  // define the types of iterators through the points and their coordinates
  // ----------------------------------------------------------------------
  typedef std::vector<std::vector<double> >::const_iterator PointIterator; 
  typedef std::vector<double>::const_iterator CoordIterator;

  // create an instance of Miniball
  // ------------------------------
  typedef Miniball::
    Miniball <Miniball::CoordAccessor<PointIterator, CoordIterator> > 
    MB;

  MB mb (2, pts.begin(), pts.end());
  
  const Point2D midpt(mb.center()[0], mb.center()[1]);
  const double radius = std::sqrt(mb.squared_radius());
  //"Computation time was "<< mb.get_time() << " seconds\n";

  MultiPolygon circle;
  circle.emplace_back();             //Make a polygon
  circle.back().emplace_back();      //Make a ring
  auto &ring = circle.back().back(); //Get the ring

  //Make a "circle"
  for(int i=0;i<CIRCLE_PT_COUNT;i++)
    ring.emplace_back(
      midpt.x+radius*std::cos(-2*M_PI*i/(double)CIRCLE_PT_COUNT),
      midpt.y+radius*std::sin(-2*M_PI*i/(double)CIRCLE_PT_COUNT)
    );
  //Close the "circle"
  ring.push_back(ring.front());

  return circle;
}



MultiPolygon GetBoundingCircleMostDistant(const MultiPolygon &mp){
  //Number of unique points from which to construct the circle. The circle will
  //have one more point of than this in order to form a closed ring).
  const int CIRCLE_PT_COUNT = 1000;

  const auto dist_pts = MostDistantPoints(mp);

  const Point2D &mpa = dist_pts.first;
  const Point2D &mpb = dist_pts.second;

  const Point2D midpt( (mpa.x+mpb.x)/2. , (mpa.y+mpb.y)/2. );
  const auto radius = EuclideanDistance(mpa,mpb)/2.;

  MultiPolygon circle;
  circle.emplace_back();             //Make a polygon
  circle.back().emplace_back();      //Make a ring
  auto &ring = circle.back().back(); //Get the ring

  //Make a "circle"
  for(int i=0;i<CIRCLE_PT_COUNT;i++)
    ring.emplace_back(
      midpt.x+radius*std::cos(-2*M_PI*i/(double)CIRCLE_PT_COUNT),
      midpt.y+radius*std::sin(-2*M_PI*i/(double)CIRCLE_PT_COUNT)
    );
  //Close the "circle"
  ring.push_back(ring.front());

  return circle;
}





}

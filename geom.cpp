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

double EuclideanDistanceSquared(const Point2D &a, const Point2D &b){
  const double xd = (a.x-b.x);
  const double yd = (a.y-b.y);
  return xd*xd+yd*yd;
}

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


BoundingBox::BoundingBox(double minx, double miny, double maxx, double maxy){
  min[0] = minx;
  min[1] = miny;
  max[0] = maxx;
  max[1] = maxy;
}

double& BoundingBox::xmin()       { return min[0]; }
double& BoundingBox::ymin()       { return min[1]; }
double& BoundingBox::xmax()       { return max[0]; }
double& BoundingBox::ymax()       { return max[1]; }
double  BoundingBox::xmin() const { return min[0]; }
double  BoundingBox::ymin() const { return min[1]; }
double  BoundingBox::xmax() const { return max[0]; }
double  BoundingBox::ymax() const { return max[1]; }
BoundingBox& BoundingBox::expand(const double expandby) {
  xmin() -= expandby;
  ymin() -= expandby;
  xmax() += expandby;
  ymax() += expandby;
  return *this;
}

BoundingBox BoundingBox::operator+(const BoundingBox &b) const {
  BoundingBox temp;
  temp.xmin() = std::min( xmin(), b.xmin() );
  temp.ymin() = std::min( ymin(), b.ymin() );
  temp.xmax() = std::max( xmax(), b.xmax() );
  temp.ymax() = std::max( ymax(), b.ymax() );
  return temp;
}

BoundingBox& BoundingBox::operator+=(const BoundingBox& b){
  xmin() = std::min( xmin(), b.xmin() );
  ymin() = std::min( ymin(), b.ymin() );
  xmax() = std::max( xmax(), b.xmax() );
  ymax() = std::max( ymax(), b.ymax() );
  return *this;
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

BoundingBox MultiPolygon::bbox() const {
  std::cerr<<"Deprecated in favour of compactnesslib::bbox()"<<std::endl;
  BoundingBox bb;
  for(const auto &p: *this)
  for(const auto &r: p)
  for(const auto &pt: r){
    bb.xmin() = std::min(bb.xmin(),pt.x);
    bb.xmax() = std::max(bb.xmax(),pt.x);
    bb.ymin() = std::min(bb.ymin(),pt.y);
    bb.ymax() = std::max(bb.ymax(),pt.y);
  }

  return bb;
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


/**
  Creates a clipper path for each element of the geocollection. This speeds
  processing later since these paths will be used to find intersection areas
  and other vector operations.
*/
void GeoCollection::clipperify() {
  for(unsigned int i=0;i<v.size();i++)
    v[i].clipper_paths = ConvertToClipper(v[i], false);
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

unsigned holeCount(const Polygon &p){
  return p.size()-1;
}

unsigned polyCount(const MultiPolygon &mp){
  return mp.size();
}

unsigned holeCount(const MultiPolygon &mp){
  return std::accumulate(mp.begin(), mp.end(), 0, [](const double b, const Polygon &p){ return b+holeCount(p);});
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



cl::Paths ConvertToClipper(const Ring &ring, const bool reversed){
  cl::Paths clipper_paths(1);

  auto &path = clipper_paths.at(0);

  if(!reversed){
    for(const auto &pt: ring)
      path.emplace_back((long long)pt.x,(long long)pt.y);
  } else {
    for(auto pt=ring.rbegin();pt!=ring.rend();pt++)
      path.emplace_back((long long)pt->x,(long long)pt->y);    
  }

  return clipper_paths;
}


cl::Paths ConvertToClipper(const MultiPolygon &mp, const bool reversed) {
  cl::Paths clipper_paths;

  for(const auto &poly: mp){
    //Send in outer perimter
    clipper_paths.push_back(ConvertToClipper(poly.at(0), reversed).front());

    //Send in the holes
    for(unsigned int i=1;i<poly.size();i++)
      clipper_paths.push_back(ConvertToClipper(poly.at(i), !reversed).front());
  }

  return clipper_paths;
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



cl::Paths BufferPath(const cl::Paths &paths, const int pad_amount){
  cl::Paths result;
  cl::ClipperOffset co;
  co.AddPaths(paths, cl::jtRound,  cl::etClosedPolygon);
  //co.AddPaths(paths, cl::jtSquare, cl::etClosedPolygon); //jtSquare runs about twice as fast as jtRound
  co.Execute(result, pad_amount);
  return result;
}



template<>
unsigned PointCount<Ring>(const Ring &r){
  return r.size();
}


Point2D CentroidPTSH(const MultiPolygon &mp){
  Point2D centroid(0,0);
  unsigned int ptcount = 0;
  for(const auto &poly: mp)
  for(const auto &pt: poly.at(0)){
    centroid.x += pt.x;
    centroid.y += pt.y;
    ptcount++;
  }

  centroid.x /= ptcount;
  centroid.y /= ptcount;

  return centroid;
}










////////////////////////////////////////////////////////////////////////////////

//https://stackoverflow.com/a/11427699/752843

// minimum distance (squared) between vertices, i.e. minimum segment length (squared)
const double EPSILON_MIN_VERTEX_DISTANCE_SQUARED = 0.00000001;

// An arbitrary tiny epsilon.  If you use float instead of double, you'll probably want to change this to something like 1E-7f
const double EPSILON_TINY = 1.0E-14;

// Arbitrary general epsilon.  Useful for places where you need more "slop" than EPSILON_TINY (which is most places).
// If you use float instead of double, you'll likely want to change this to something like 1.192092896E-04
const double EPSILON_GENERAL = 1.192092896E-07;

bool AreValuesEqual(
  const double val1, 
  const double val2, 
  const double tolerance
){
  return (val1 >= (val2 - tolerance) && val1 <= (val2 + tolerance));
}


double PointToPointDistanceSquared(
  const double p1x, const double p1y, 
  const double p2x, const double p2y
){
  double dx = p2x - p1x;
  double dy = p2y - p1y;
  return (dx * dx) + (dy * dy);
}


double PointSegmentDistanceSquared( 
  const double px,  const double py,
  const double p1x, const double p1y,
  const double p2x, const double p2y,
  double& t,
  double& qx, double& qy)
{
  const double dx   = p2x - p1x;
  const double dy   = p2y - p1y;
  const double dp1x = px - p1x;
  const double dp1y = py - p1y;

  const double segLenSquared = (dx * dx) + (dy * dy);
  
  if (AreValuesEqual(segLenSquared, 0.0, EPSILON_MIN_VERTEX_DISTANCE_SQUARED)){
    // segment is a point.
    qx = p1x;
    qy = p1y;
    t = 0.0;
    return ((dp1x * dp1x) + (dp1y * dp1y));
  } else {
    t = ((dp1x * dx) + (dp1y * dy)) / segLenSquared;
    if (t <= EPSILON_TINY){
      // intersects at or to the "left" of first segment vertex (p1x, p1y).  If t is approximately 0.0, then
      // intersection is at p1.  If t is less than that, then there is no intersection (i.e. p is not within
      // the 'bounds' of the segment)
      if (t >= -EPSILON_TINY){
        // intersects at 1st segment vertex
        t = 0.0;
      }
      // set our 'intersection' point to p1.
      qx = p1x;
      qy = p1y;
      // Note: If you wanted the ACTUAL intersection point of where the projected lines would intersect if
      // we were doing PointLineDistanceSquared, then qx would be (p1x + (t * dx)) and qy would be (p1y + (t * dy)).
    } else if (t >= (1.0 - EPSILON_TINY)) {
      // intersects at or to the "right" of second segment vertex (p2x, p2y).  If t is approximately 1.0, then
      // intersection is at p2.  If t is greater than that, then there is no intersection (i.e. p is not within
      // the 'bounds' of the segment)
      if (t <= (1.0 + EPSILON_TINY)){
        // intersects at 2nd segment vertex
        t = 1.0;
      }
      qx = p2x;
      qy = p2y;
      // Note: If you wanted the ACTUAL intersection point of where the projected lines would intersect if
      // we were doing PointLineDistanceSquared, then qx would be (p1x + (t * dx)) and qy would be (p1y + (t * dy)).
    } else {
      // The projection of the point to the point on the segment that is perpendicular succeeded and the point
      // is 'within' the bounds of the segment.  Set the intersection point as that projected point.
      qx = ((1.0 - t) * p1x) + (t * p2x);
      qy = ((1.0 - t) * p1y) + (t * p2y);
      // for debugging
      //ASSERT(AreValuesEqual(qx, p1x + (t * dx), EPSILON_TINY));
      //ASSERT(AreValuesEqual(qy, p1y + (t * dy), EPSILON_TINY));
    }
    // return the squared distance from p to the intersection point.
    double dpqx = px - qx;
    double dpqy = py - qy;
    return ((dpqx * dpqx) + (dpqy * dpqy));
  }
}


double SegmentSegmentDistanceSquared(   
  const double p1x, const double p1y,
  const double p2x, const double p2y,
  const double p3x, const double p3y,
  const double p4x, const double p4y,
  double& qx, double& qy
){
  // check to make sure both segments are long enough (i.e. verts are farther apart than minimum allowed vert distance).
  // If 1 or both segments are shorter than this min length, treat them as a single point.
  const double segLen12Squared = PointToPointDistanceSquared(p1x, p1y, p2x, p2y);
  const double segLen34Squared = PointToPointDistanceSquared(p3x, p3y, p4x, p4y);
  double t = 0.0;
  double minDist2 = 1E+38;
  if (segLen12Squared <= EPSILON_MIN_VERTEX_DISTANCE_SQUARED){
    qx = p1x;
    qy = p1y;
    if (segLen34Squared <= EPSILON_MIN_VERTEX_DISTANCE_SQUARED){
      // point to point
      minDist2 = PointToPointDistanceSquared(p1x, p1y, p3x, p3y);
    } else {
      // point - seg
      minDist2 = PointSegmentDistanceSquared(p1x, p1y, p3x, p3y, p4x, p4y, t, qx, qy);
    }
    return minDist2;
  } else if (segLen34Squared <= EPSILON_MIN_VERTEX_DISTANCE_SQUARED) {
    // seg - point
    minDist2 = PointSegmentDistanceSquared(p3x, p3y, p1x, p1y, p2x, p2y, t, qx, qy);
    return minDist2;
  }

  // if you have a point class and/or methods to do cross products, you can use those here.
  // This is what we're actually doing:
  // Point2D delta43(p4x - p3x, p4y - p3y);    // dir of p3 -> p4
  // Point2D delta12(p1x - p2x, p1y - p2y);    // dir of p2 -> p1
  // double d = delta12.Cross2D(delta43);
  double d = ((p4y - p3y) * (p1x - p2x)) - ((p1y - p2y) * (p4x - p3x));
  bool bParallel = AreValuesEqual(d, 0.0, EPSILON_GENERAL);

  if (!bParallel){
    // segments are not parallel.  Check for intersection.
    // Point2D delta42(p4x - p2x, p4y - p2y);    // dir of p2 -> p4
    // t = 1.0 - (delta42.Cross2D(delta43) / d);
    t = 1.0 - ((((p4y - p3y) * (p4x - p2x)) - ((p4y - p2y) * (p4x - p3x))) / d);
    double seg12TEps = sqrt(EPSILON_MIN_VERTEX_DISTANCE_SQUARED / segLen12Squared);
    if (t >= -seg12TEps && t <= (1.0 + seg12TEps)){
      // inside [p1,p2].   Segments may intersect.
      // double s = 1.0 - (delta12.Cross2D(delta42) / d);
      double s = 1.0 - ((((p4y - p2y) * (p1x - p2x)) - ((p1y - p2y) * (p4x - p2x))) / d);
      double seg34TEps = sqrt(EPSILON_MIN_VERTEX_DISTANCE_SQUARED / segLen34Squared);
      if (s >= -seg34TEps && s <= (1.0 + seg34TEps)){
        // segments intersect!
        minDist2 = 0.0;
        qx = ((1.0 - t) * p1x) + (t * p2x);
        qy = ((1.0 - t) * p1y) + (t * p2y);
        // for debugging
        //double qsx = ((1.0 - s) * p3x) + (s * p4x);
        //double qsy = ((1.0 - s) * p3y) + (s * p4y);
        //ASSERT(AreValuesEqual(qx, qsx, EPSILON_MIN_VERTEX_DISTANCE_SQUARED));
        //ASSERT(AreValuesEqual(qy, qsy, EPSILON_MIN_VERTEX_DISTANCE_SQUARED));
        return minDist2;
      }
    }
  }

  // Segments do not intersect.   Find closest point and return dist.   No other way at this
  // point except to just brute-force check each segment end-point vs opposite segment.  The
  // minimum distance of those 4 tests is the closest point.
  double tmpQx, tmpQy, tmpD2;
  minDist2 = PointSegmentDistanceSquared(p3x, p3y, p1x, p1y, p2x, p2y, t, qx, qy);
  tmpD2 = PointSegmentDistanceSquared(p4x, p4y, p1x, p1y, p2x, p2y, t, tmpQx, tmpQy);
  if (tmpD2 < minDist2){
    qx = tmpQx;
    qy = tmpQy;
    minDist2 = tmpD2;
  }
  tmpD2 = PointSegmentDistanceSquared(p1x, p1y, p3x, p3y, p4x, p4y, t, tmpQx, tmpQy);
  if (tmpD2 < minDist2){
    qx = p1x;
    qy = p1y;
    minDist2 = tmpD2;
  }
  tmpD2 = PointSegmentDistanceSquared(p2x, p2y, p3x, p3y, p4x, p4y, t, tmpQx, tmpQy);
  if (tmpD2 < minDist2){
    qx = p2x;
    qy = p2y;
    minDist2 = tmpD2;
  }

  return minDist2;
}

////////////////////////////////////////////////////////////////////////////////

double SegmentSegmentDistanceSquared(
  const Point2D &a1, const Point2D &a2,
  const Point2D &b1, const Point2D &b2
){
  double closest_x_on_1;
  double closest_y_on_1;
  return SegmentSegmentDistanceSquared(
    a1.x, a1.y,
    a2.x, a2.y,
    b1.x, b2.y,
    b2.x, b2.y,
    closest_x_on_1, closest_y_on_1
  );
}



//Info: https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
//Info: https://stackoverflow.com/questions/8721406/how-to-determine-if-a-point-is-inside-a-2d-convex-polygon/23223947#23223947
//Info: http://stackoverflow.com/a/2922778/752843
bool ContainsPoint(const Ring &ring, const Point2D &pt){
  const auto &ringpts = ring.v;
  unsigned int i, j;
  int c = 0;
  for (i = 0, j = ringpts.size()-1; i < ringpts.size(); j = i++) {
    if ( ((ringpts[i].y>pt.y) != (ringpts[j].y>pt.y)) &&
     (pt.x < (ringpts[j].x-ringpts[i].x) * (pt.y-ringpts[i].y) / (ringpts[j].y-ringpts[i].y) + ringpts[i].x) )
       c = !c;
  }
  return c;
}

bool ContainsPoint(const Polygon      &poly,  const Point2D &pt){
  for(const auto &ring: poly.v){
    if(ContainsPoint(ring,pt))
      return true;
  }
  return false;
}

bool ContainsPoint(const MultiPolygon &mp, const Point2D &pt){
  for(const auto &poly: mp.v){
    if(ContainsPoint(poly,pt))
      return true;
  }
  return false;
}







BoundingBox bbox(const Ring         &r ){
  BoundingBox bb;
  for(const auto &pt: r){
    bb.xmin() = std::min(bb.xmin(),pt.x);
    bb.xmax() = std::max(bb.xmax(),pt.x);
    bb.ymin() = std::min(bb.ymin(),pt.y);
    bb.ymax() = std::max(bb.ymax(),pt.y);
  }

  return bb;
}

BoundingBox bbox(const Polygon      &p ){
  BoundingBox bb;
  for(const auto &r: p)
    bb += bbox(r);

  return bb;
}

BoundingBox bbox(const MultiPolygon &mp){
  BoundingBox bb;
  for(const auto &p: mp)
    bb += bbox(p);

  return bb;
}







}

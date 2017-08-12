#ifndef DOCTEST_CONFIG_DISABLE
  #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#endif
#include "doctest.h"
#include "compactengine.hpp"
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <iostream>

namespace complib {


  // To find orientation of ordered triplet (p, q, r).
  // The function returns following values
  // 0 --> p, q and r are colinear
  // 1 --> Clockwise
  // 2 --> Counterclockwise
  int FindOrientation(
    const double px,
    const double py,
    const double qx,
    const double qy,
    const double rx,
    const double ry
  ){
    const int val = (qy-py) * (rx-qx) - (qx-px) * (ry-qy);

    if (val == 0) 
      return 0;   //Colinear
    else if (val>0)
      return 1;   //Clockwise
    else
      return 2;   //Counter-clockwise
  }
   
  void FindConvexHull(
    const double *const  x,
    const double *const  y,
    const PointCount     N,
    std::vector<double> &hx,
    std::vector<double> &hy
  ){
    if (N < 3)
      throw std::runtime_error("There must be at least 3 points for a convex hull!");

    // Find the leftmost point
    unsigned int l = 0;
    for (unsigned int i=1;i<N;i++)
      if (x[i] < x[l])
        l = i;

    // Start from leftmost point, keep moving counterclockwise
    // until reach the start point again.  This loop runs O(h)
    // times where h is number of points in result or output.
    unsigned int p = l;
    unsigned int q;
    do{
      // Add current point to result
      hx.emplace_back(x[p]);
      hy.emplace_back(y[p]);

      // Search for a point 'q' such that orientation(p, x,
      // q) is counterclockwise for all points 'x'. The idea
      // is to keep track of last visited most counterclock-
      // wise point in q. If any point 'i' is more counterclock-
      // wise than q, then update q.
      q = (p+1)%N;
      for (int i = 0; i < N; i++){
        // If i is more counterclockwise than current q, then
        // update q
        if (FindOrientation(x[p],y[p], x[i],y[i], x[q],y[q]) == 2)
          q = i;
      }

      // Now q is the most counterclockwise with respect to p
      // Set p as q for next iteration, so that q is added to
      // result 'hull'
      p = q;
    } while (p != l);  // While we don't come to first point
  }


  //if return >0 then point R is upper of lineseg PQ else down
  double rcCross(
    const double Px,
    const double Py,
    const double Qx,
    const double Qy,
    const double Rx,
    const double Ry
  ){
    return  (Qx-Px)*(Ry-Py) - (Rx-Px)*(Qy-Py);
  }

  void rcHull(
    std::vector<Point2D> &xy,
    std::vector<Point2D> &L,
    std::vector<Point2D> &U
  ){
    PointCount j=0;
    PointCount k=0;

    std::sort(xy.begin(),xy.end(), [](const Point2D &a, const Point2D &b){ return a.x<b.x; });

    const auto N = xy.size();

    L.resize(2*N);
    U.resize(2*N);

    for(int i=0;i<xy.size();i++){
      while(j>=2 && rcCross(L[j-2].x,L[j-2].y, L[j-1].x,L[j-1].y, xy[i].x,xy[i].y)<=0)//p[i] is making right turn we need left turn
        j--;
      while(k>=2 && rcCross(U[k-2].x,U[k-2].y, U[k-1].x,U[k-1].y, xy[i].x,xy[i].y)>=0)//p[i] is making left turn we need right
        k--;
      L[j++] = xy[i];
      U[k++] = xy[i];
    }

    L.resize(j);
    U.resize(k);
  }
 
  double pFarthestPoints(
    const double *const x,
    const double *const y,
    PointCount N
  ){
    double t;
    double n;
    double k1;
    double k2;

    std::vector<Point2D> xy;
    for(unsigned int i=0;i<N;i++)
      xy.emplace_back(x[i],y[i]);

    std::vector<Point2D> U;
    std::vector<Point2D> L;

    rcHull(xy,L,U);
            
    PointCount i = 0;
    PointCount j = L.size()-1;
    PointCount m = U.size()-1;
    double dist  = -1;
    while(i<m || j>0){
      dist = std::max(dist,EuclideanDistance(U[i],L[j]));
      if(i==m)
        j--;
      else if(j==0)
        i++;
      else {
        if ( (U[i+1].y-U[i].y) * (L[j].x-L[j-1].x) > (L[j].y-L[j-1].y) * (U[i+1].x-U[i].x) )
          i++;
        else
          j--;
      }
    }

    return dist;
  }

















  double pPerimeter(const double *const x, const double *const y, const PointCount N){
    double dist = 0;
    for(int i=0;i<N-1;i++){
      dist += EuclideanDistance(x[i],y[i],x[i+1],y[i+1]);
    }
    dist += EuclideanDistance(x[0],y[0],x[N-1],y[N-1]);
    return dist;
  }

  double pPolygonArea(const double *const x, const double *const y, const PointCount N){
    double area    = 0;
    unsigned int j = N-1;
    for(unsigned int i=0;i<N;i++){
      area += (x[j] + x[i]) * (y[j] - y[i]);
      j = i;
    }
    area /=2;
    return std::abs(area);
  }









  double pScorePolsbyPopper(const double *const x, const double *const y, const PointCount N){
    const double area  = pPolygonArea(x,y,N);
    const double perim = pPerimeter(x,y,N);
    return 4*M_PI*area/perim/perim;
  }

  double pScoreSchwartzberg(const double *const x, const double *const y, const PointCount N){
    const double perim  = pPerimeter(x,y,N);
    const double area   = pPolygonArea(x,y,N);
    const double radius = std::sqrt(area/M_PI);
    const double circum = 2*M_PI*radius;
    return circum/perim;
  }

  double pScoreConvexHull(const double *const x, const double *const y, const PointCount N){
    std::vector<double> hx;
    std::vector<double> hy;
    FindConvexHull(x,y,N,hx,hy);
    const double area      = pPolygonArea(x,y,N);
    const double hull_area = PolygonArea(hx,hy);
    return area/hull_area;
  }

  double pScoreReock(
    const double *const x,
    const double *const y,
    const PointCount N
  ){
    const double area      = pPolygonArea(x,y,N);
    const double radius    = pFarthestPoints(x,y,N)/2.0;
    const double circ_area = M_PI*radius*radius;

    return area/circ_area;
  }






  double Perimeter(const std::vector<double> &x, const std::vector<double> &y){
    return pPerimeter(x.data(),y.data(),x.size());
  }

  double PolygonArea(const std::vector<double> &x, const std::vector<double> &y){
    return pPolygonArea(x.data(),y.data(),x.size());
  }

  double ScorePolsbyPopper(const std::vector<double> &x, const std::vector<double> &y){
    return pScorePolsbyPopper(x.data(),y.data(),x.size());
  }

  double ScoreSchwartzberg(const std::vector<double> &x, const std::vector<double> &y){
    return pScoreSchwartzberg(x.data(),y.data(),x.size());
  }

  double ScoreConvexHull(const std::vector<double> &x, const std::vector<double> &y){
    return pScoreConvexHull(x.data(),y.data(),x.size());
  }

  double ScoreReock(const std::vector<double> &x, const std::vector<double> &y){
    return pScoreReock(x.data(),y.data(),x.size());
  }


  std::vector< std::pair<double,double> > Multifier(
    const std::vector<double> &x,
    const std::vector<double> &y,
    const std::vector<double> &id,
    std::function<double(const double *const x, const double *const y, const PointCount N)> func
  ){
    std::vector< std::pair<double,double> > ret;
    unsigned int start = 0;
    for(unsigned int end = 0;end<x.size();end++){
      if(id[start]!=id[end]){
        ret.emplace_back(id[start],func(x.data()+start,y.data()+start,end-start));
        start = end;
      }
    }
    ret.emplace_back(id[start],func(x.data()+start,y.data()+start,x.size()-start));
    return ret; 
  }


  std::vector< std::pair<double,double> > PerimeterMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id){
    return Multifier(x, y, id, pPerimeter);
  }

  std::vector< std::pair<double,double> > PolygonAreaMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id){
    return Multifier(x, y, id, pPolygonArea);
  }

  std::vector< std::pair<double,double> > ScorePolsbyPopperMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id){
    return Multifier(x, y, id, pScorePolsbyPopper);
  }

  std::vector< std::pair<double,double> > ScoreSchwartzbergMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id){
    return Multifier(x, y, id, pScoreSchwartzberg);
  }

  std::vector< std::pair<double,double> > ScoreConvexHullMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id){
    return Multifier(x, y, id, pScoreConvexHull);
  }

  std::vector< std::pair<double,double> > ScoreReockMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id){
    return Multifier(x, y, id, pScoreReock);
  }

}




TEST_CASE("CountTEST"){
  std::vector<double> x  = {{1,3,3,1, 3,5,5,3, 1,3,3,1, 3,5,5,3}};
  std::vector<double> y  = {{1,1,3,3, 1,1,3,3, 3,3,5,5, 3,3,5,5}};
  std::vector<double> id = {{1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4}};

  SUBCASE("Area"){
    CHECK(complib::pPolygonArea(x.data(),y.data(),4)==4);
  }

  SUBCASE("ID count"){
    auto ret = complib::PolygonAreaMulti(x,y,id);
    CHECK(ret.size()==4);
    for(const auto &r: ret){
      std::cout<<r.first<<" "<<r.second<<std::endl;
    }
  }
}
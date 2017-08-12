#include "GeoStuff.hpp"
#include "doctest.h"
#include <cmath>
#include <iostream>
#include <cassert>

const double RAD_TO_DEG = 180.0/M_PI;
const double Rearth = 6371; //km

double EuclideanDistance(
  const double x1,
  const double y1, 
  const double x2, 
  const double y2
){
  const double xd = x2-x1;
  const double yd = y2-y1;
  return std::sqrt(xd*xd+yd*yd);
}

double EuclideanDistance2(
  const double x1,
  const double y1, 
  const double x2, 
  const double y2
){
  const double xd = x2-x1;
  const double yd = y2-y1;
  return xd*xd+yd*yd;
}


double EuclideanDistance(const Point2D &a, const Point2D &b){
  return EuclideanDistance(a.x,a.y,b.x,b.y);
}

double EuclideanDistance2(const Point2D &a, const Point2D &b){
  return EuclideanDistance2(a.x,a.y,b.x,b.y);
}

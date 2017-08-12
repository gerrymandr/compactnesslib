#ifndef _geostuff_hpp_
#define _geostuff_hpp_

#include "Point.hpp"

double EuclideanDistance(const Point2D &a, const Point2D &b);
double EuclideanDistance2(const Point2D &a, const Point2D &b);

double EuclideanDistance(const double lon1, const double lat1, const double lon2, const double lat2);
double EuclideanDistance2(const double lon1, const double lat1, const double lon2, const double lat2);

#endif

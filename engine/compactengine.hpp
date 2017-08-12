#ifndef _compactengine_hpp_
#define _compactengine_hpp_

#include "GeoStuff.hpp"
#include "Point.hpp"
#include <vector>

typedef unsigned int PointCount;

namespace complib {
  double Perimeter(const std::vector<double> &x, const std::vector<double> &y);
  double PolygonArea(const std::vector<double> &x, const std::vector<double> &y);
  double ScorePolsbyPopper(const std::vector<double> &x, const std::vector<double> &y);
  double ScoreSchwartzberg(const std::vector<double> &x, const std::vector<double> &y);
  double ScoreConvexHull(const std::vector<double> &x, const std::vector<double> &y);
  double ScoreReock(const std::vector<double> &x, const std::vector<double> &y);

  std::vector< std::pair<double,double> > PerimeterMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id);
  std::vector< std::pair<double,double> > PolygonAreaMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id);
  std::vector< std::pair<double,double> > ScorePolsbyPopperMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id);
  std::vector< std::pair<double,double> > ScoreSchwartzbergMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id);
  std::vector< std::pair<double,double> > ScoreConvexHullMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id);
  std::vector< std::pair<double,double> > ScoreReockMulti(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &id);
}

#endif
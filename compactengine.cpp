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

double ScorePolsbyPopper(const MultiPolygon &mp){
  const double area  = mp.area();
  const double perim = mp.perim();
  return 4*M_PI*area/perim/perim;
}

double ScoreSchwartzberg(const MultiPolygon &mp){
  const double perim  = mp.perim();
  const double area   = mp.area();
  const double radius = std::sqrt(area/M_PI);
  const double circum = 2*M_PI*radius;
  return circum/perim;
}

double ScoreConvexHull(const MultiPolygon &mp){
  const double area      = mp.area();
  const double hull_area = mp.hullArea();
  return area/hull_area;
}

double ScoreReock(const MultiPolygon &mp){
  const double area      = mp.area();
  const double radius    = mp.diameter()/2;
  const double circ_area = M_PI*radius*radius;

  return area/circ_area;
}

void CalculateAllScores(GeoCollection &mps){
  #pragma omp parallel for
  for(unsigned int i=0;i<mps.size();i++){
    auto &mp = mps[i];
    mp.props["perim"]      = mp.perim();
    mp.props["area"]       = mp.area();
    mp.props["PolsbyPopp"] = complib::ScorePolsbyPopper(mp);
    mp.props["Schwartzbe"] = complib::ScoreSchwartzberg(mp);
    mp.props["ConvexHull"] = complib::ScoreConvexHull  (mp);
    mp.props["Reock"]      = complib::ScoreReock       (mp);
  }
}

}


/*
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
*/
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


//TODO: Use "https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html"
double ScoreReock(const MultiPolygon &mp){
  const double area      = mp.area();
  const double radius    = mp.diameter()/2;
  const double circ_area = M_PI*radius*radius;

  return area/circ_area;
}

void CalculateAllScores(GeoCollection &mps){
  CalculateListOfScores(mps, score_names);
}


void CalculateScoreFromString(MultiPolygon &mp, const std::string score){
  if(score=="perim")
    mp.scores[score] = mp.perim();
  else if(score=="area")
    mp.scores[score] = mp.area();
  else if(score=="PolsbyPopp")
    mp.scores[score] = complib::ScorePolsbyPopper(mp);
  else if(score=="Schwartzbe")
    mp.scores[score] = complib::ScoreSchwartzberg(mp);
  else if(score=="ConvexHull")
    mp.scores[score] = complib::ScoreConvexHull  (mp);
  else if(score=="Reock")
    mp.scores[score] = complib::ScoreReock       (mp);
  else 
    throw std::runtime_error("Unrecognized score name '" + score + "'!");
}

void CalculateListOfScores(GeoCollection &gc, const std::vector<std::string> &score_list){
  #pragma omp parallel for
  for(unsigned int i=0;i<gc.size();i++){
    for(const auto &s: score_list)
      CalculateScoreFromString(gc[i],s);
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
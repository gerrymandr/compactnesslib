#include "doctest.h"
#include "compactengine.hpp"
#include "geom.hpp"
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <iostream>

namespace complib {

double ScorePolsbyPopper(const MultiPolygon &mp){
  const double area  = areaIncludingHoles(mp);
  const double perim = perimExcludingHoles(mp);
  return 4*M_PI*area/perim/perim;
}

double ScoreSchwartzberg(const MultiPolygon &mp){
  const double area   = areaIncludingHoles(mp);
  const double perim  = perimExcludingHoles(mp);
  const double radius = std::sqrt(area/M_PI);
  const double circum = 2*M_PI*radius;
  return circum/perim;
}

double ScoreConvexHull(const MultiPolygon &mp){
  const double area      = areaIncludingHoles(mp);
  const double hull_area = hullAreaPolygonOuterRings(mp);
  return area/hull_area;
}

double ScoreConvexHullPTB(const MultiPolygon &mp, const MultiPolygon &border){
  const double area      = areaIncludingHoles(mp);
  const double hull_area = IntersectionArea(mp.getHull(),border);
  return area/hull_area;
}


//TODO: Use "https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html"
double ScoreReock(const MultiPolygon &mp){
  const double area      = areaIncludingHoles(mp);
  const double radius    = diameterOfEntireMultiPolygon(mp)/2;
  const double circ_area = M_PI*radius*radius;

  return area/circ_area;
}

void CalculateAllScores(GeoCollection &mps){
  CalculateListOfScores(mps, score_names);
}


void CalculateScoreFromString(MultiPolygon &mp, const std::string score){
  if(score=="perim")
    mp.scores[score] = perimExcludingHoles(mp);
  else if(score=="area")
    mp.scores[score] = areaIncludingHoles(mp);
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

void CalculateListOfScores(GeoCollection &gc, std::vector<std::string> score_list){
  if(score_list.empty())
    score_list = score_names;
  else if(score_list.size()==1 && score_list.at(0)=="all")
    score_list = score_names;

  for(unsigned int i=0;i<gc.size();i++){
    for(const auto &s: score_list)
      CalculateScoreFromString(gc[i],s);
  }
}

}

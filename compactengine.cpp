#include "doctest.h"
#include "compactengine.hpp"
#include "geom.hpp"
#include <cmath>
#include <vector>
#include <stdexcept>

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

void CalculateAllUnboundedScores(GeoCollection &mps){
  CalculateListOfUnboundedScores(mps, getListOfUnboundedScores());
}

void CalculateListOfUnboundedScores(GeoCollection &gc, std::vector<std::string> score_list){
  if(score_list.empty())
    score_list = getListOfUnboundedScores();
  else if(score_list.size()==1 && score_list.at(0)=="all")
    score_list = getListOfUnboundedScores();

  for(unsigned int i=0;i<gc.size();i++){
    for(const auto &sn: score_list){
      if(!unbounded_score_map.count(sn))
        throw std::runtime_error("Unrecognized score name '" + sn + "'!");
      gc[i].scores[sn] = unbounded_score_map.at(sn)(gc[i]);
    }
  }
}

const std::vector<std::string>& getListOfUnboundedScores(){
  static std::vector<std::string> score_names;
  for(const auto &kv: unbounded_score_map)
    score_names.push_back(kv.first);
  return score_names;
}

const score_map_t unbounded_score_map({
  // {"areaAH",            areaIncludingHoles},
  // {"perimSH",           perimExcludingHoles},
  // {"ScorePolsbyPopper", ScorePolsbyPopper},
  // {"ScoreSchwartzberg", ScoreSchwartzberg},
  {"ScoreConvexHull",   ScoreConvexHull},
  {"ScoreReock",        ScoreReock}
});

}

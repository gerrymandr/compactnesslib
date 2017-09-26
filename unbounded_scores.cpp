#include "unbounded_scores.hpp"
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

double ScoreConvexHullPS(const MultiPolygon &mp){
  const double area      = areaIncludingHoles(mp);
  const double hull_area = hullAreaPolygonOuterRings(mp);
  return area/hull_area;
}

double ScoreConvexHullPT(const MultiPolygon &mp){
  const double area_mp   = areaIncludingHoles(mp);
  const double hull_area = area(mp.getHull());
  return area_mp/hull_area;
}

//TODO: Use "https://people.inf.ethz.ch/gaertner/subdir/software/miniball.html"
double ScoreReockPT(const MultiPolygon &mp){
  const double area      = areaIncludingHoles(mp);
  const double radius    = diameterOfEntireMultiPolygon(mp)/2;
  const double circ_area = M_PI*radius*radius;

  return area/circ_area;
}

double ScoreReockPS(const MultiPolygon &mp){
  const double area = areaIncludingHoles(mp);
  double circ_area  = 0;
  for(const auto &poly: mp)
    circ_area += M_PI*std::pow(diameterOuter(poly)/2.0,2.0);

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
      if(unbounded_score_map.count(sn))
        gc[i].scores[sn] = unbounded_score_map.at(sn)(gc[i]);
    }
  }
}

const std::vector<std::string>& getListOfUnboundedScores(){
  static std::vector<std::string> score_names;
  if(!score_names.empty())
    return score_names;
  for(const auto &kv: unbounded_score_map)
    score_names.push_back(kv.first);
  return score_names;
}

const unbounded_score_map_t unbounded_score_map({
  {"areaAH",     [](const MultiPolygon &mp) { return areaIncludingHoles(mp);  }},
  {"perimSH",    [](const MultiPolygon &mp) { return perimExcludingHoles(mp); }},
  {"HoleCount",  [](const MultiPolygon &mp) { return holeCount(mp); }},
  {"PolyCount",  [](const MultiPolygon &mp) { return polyCount(mp); }},
  {"PolsbyPopp", ScorePolsbyPopper},
  {"Schwartzbe", ScoreSchwartzberg},
  {"CvxHullPS",  ScoreConvexHullPS},
  {"CvxHullPT",  ScoreConvexHullPT},
  {"ReockPT",    ScoreReockPT},
  {"ReockPS",    ScoreReockPS}
});

}

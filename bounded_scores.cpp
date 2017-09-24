#include "bounded_scores.hpp"
#include "geom.hpp"
#include <cmath>
#include <vector>
#include <stdexcept>
#include <unordered_map>
#include "lib/clipper.hpp"

#include <sstream>  //TODO
#include <iostream> //TODO

namespace complib {

namespace cl = ClipperLib;

double ScoreConvexHullPTB(const MultiPolygon &mp, const MultiPolygon &border){
  const double area      = areaIncludingHoles(mp);
  const double hull_area = IntersectionArea(mp.getHull(),border);
  double ratio = area/hull_area;
  if(ratio>1)
    ratio = 1;
  return ratio;
}



double ScoreReockPTB(const MultiPolygon &mp, const MultiPolygon &border){
  const auto   circle = GetBoundingCircle(mp);
  const auto   iarea  = IntersectionArea(circle, border);
  const double area   = areaIncludingHoles(mp);

  double ratio = area/iarea;
  if(ratio>1)
    ratio = 1;
  return ratio;
}



void CalculateAllBoundedScores(
  GeoCollection &subunits,
  const GeoCollection &superunits,
  const std::string join_on
){
  CalculateListOfBoundedScores(subunits, superunits, join_on, getListOfBoundedScores());
}



void CalculateListOfBoundedScores(
  GeoCollection &subunits,
  const GeoCollection &superunits,
  const std::string join_on,
  std::vector<std::string> score_list
){
  for(const auto &mp: subunits)
    if(!mp.props.count(join_on))
      throw std::runtime_error("At least one subunit was missing the joining attribute!");

  //A quick was to access superunits based on their key
  std::unordered_map<std::string, const MultiPolygon *> su_key;
  for(const auto &mp: superunits){
    if(!mp.props.count(join_on))
      throw std::runtime_error("At least one superunit was missing the joining attribute!");
    if(su_key.count(mp.props.at(join_on)))
      throw std::runtime_error("More than one superunit had the same key!");
    su_key[mp.props.at(join_on)] = &mp;
  }

  if(score_list.empty())
    score_list = getListOfBoundedScores();
  else if(score_list.size()==1 && score_list.at(0)=="all")
    score_list = getListOfBoundedScores();

  for(auto& su: subunits){
    for(const auto &sn: score_list){
      if(bounded_score_map.count(sn))
        su.scores[sn] = bounded_score_map.at(sn)(su,*su_key.at(su.props.at(join_on)));
    }
  }
}



const std::vector<std::string>& getListOfBoundedScores(){
  static std::vector<std::string> score_names;
  if(!score_names.empty())
    return score_names;
  for(const auto &kv: bounded_score_map)
    score_names.push_back(kv.first);
  return score_names;
}



const bounded_score_map_t bounded_score_map({
  {"CvxHullPTB", ScoreConvexHullPTB},
  {"ReockPTB",   ScoreReockPTB}
});

}

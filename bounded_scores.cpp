#include "doctest.h"
#include "bounded_scores.hpp"
#include "geom.hpp"
#include <cmath>
#include <vector>
#include <stdexcept>
#include <unordered_map>

namespace complib {

double ScoreConvexHullPTB(const MultiPolygon &mp, const MultiPolygon &border){
  const double area      = areaIncludingHoles(mp);
  const double hull_area = IntersectionArea(mp.getHull(),border);
  return area/hull_area;
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
  for(const auto &mp: subunits.v)
    if(!mp.props.count(join_on))
      throw std::runtime_error("At least one subunit was missing the joining attribute!");

  //A quick was to access superunits based on their key
  std::unordered_map<std::string, const MultiPolygon *> su_key;
  for(const auto &mp: superunits.v){
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

  for(auto& su: subunits.v){
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
  {"ConvexHullPTB", ScoreConvexHullPTB}
});

}

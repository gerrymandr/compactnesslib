#ifndef _bounded_scores_hpp_
#define _bounded_scores_hpp_

#include "geojson.hpp"
#include "geom.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>

namespace complib {
  const std::vector<std::string>& getListOfBoundedScores();

  double ScoreConvexHullPTB        (const MultiPolygon &mp, const MultiPolygon &border);

  void CalculateAllBoundedScores(
    GeoCollection &subunits,
    const GeoCollection &superunits,
    const std::string join_on
  );

  void CalculateListOfBoundedScores(
    GeoCollection &subunits,
    const GeoCollection &superunits,
    const std::string join_on,
    std::vector<std::string> score_list
  );

  typedef std::unordered_map<std::string, std::function<double(const MultiPolygon &subunit, const MultiPolygon &superunit)> > bounded_score_map_t;
  extern const bounded_score_map_t bounded_score_map;
}

#endif
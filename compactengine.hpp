#ifndef _compactengine_hpp_
#define _compactengine_hpp_

#include "geojson.hpp"
#include "geom.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>

namespace complib {
  const std::vector<std::string>& getListOfUnboundedScores();

  double ScorePolsbyPopper           (const MultiPolygon &mp);
  double ScoreSchwartzberg           (const MultiPolygon &mp);
  double ScoreConvexHull             (const MultiPolygon &mp);
  double ScoreConvexHullPTB          (const MultiPolygon &mp, const MultiPolygon &border);
  double ScoreReock                  (const MultiPolygon &mp);
  void CalculateAllUnboundedScores   (GeoCollection &mps);
  void CalculateListOfUnboundedScores(GeoCollection &gc, std::vector<std::string> score_list);

  typedef std::unordered_map<std::string, std::function<double(const MultiPolygon &mp)> > score_map_t;
  extern const score_map_t unbounded_score_map;
}

#endif
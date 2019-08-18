#ifndef _unbounded_scores_hpp_
#define _unbounded_scores_hpp_

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
  double ScoreConvexHullPT           (const MultiPolygon &mp);
  double ScoreReockPT                (const MultiPolygon &mp);
  double ScoreConvexHullPS           (const MultiPolygon &mp);
  double ScoreReockPS                (const MultiPolygon &mp);
  void CalculateAllUnboundedScores   (GeoCollection &mps);
  void CalculateListOfUnboundedScores(GeoCollection &gc, std::vector<std::string> score_list);

  typedef std::unordered_map<std::string, std::function<double(const MultiPolygon &mp)> > unbounded_score_map_t;
  extern const unbounded_score_map_t unbounded_score_map;
}

#endif
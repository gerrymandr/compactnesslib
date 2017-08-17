#ifndef _compactengine_hpp_
#define _compactengine_hpp_

#include "geojson.hpp"
#include "geom.hpp"
#include <string>
#include <vector>

namespace complib {
  const std::vector<std::string> score_names = {{"perim","area","PolsbyPopp","Schwartzbe","ConvexHull","Reock"}};

  double ScorePolsbyPopper(const MultiPolygon &mp);
  double ScoreSchwartzberg(const MultiPolygon &mp);
  double ScoreConvexHull  (const MultiPolygon &mp);
  double ScoreReock       (const MultiPolygon &mp);
  void CalculateAllScores(GeoCollection &mps);
}

#endif
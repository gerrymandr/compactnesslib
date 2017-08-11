#ifndef _compactengine_hpp_
#define _compactengine_hpp_

#include "geom.hpp"
#include <string>

namespace complib {
  double ScorePolsbyPopper(const MultiPolygon &mp);
  double ScoreSchwartzberg(const MultiPolygon &mp);
  double ScoreConvexHull  (const MultiPolygon &mp);
  double ScoreReock       (const MultiPolygon &mp);
  void CalculateAllScores(GeoCollection &mps);
  GeoCollection ReadGeoJSONFile(std::string filename);
  GeoCollection ReadGeoJSON(std::string filename);
}

#endif
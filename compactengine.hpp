#ifndef _compactengine_hpp_
#define _compactengine_hpp_

#include "geom.hpp"
#include "shapefile.hpp"

namespace complib {
  double ScorePolsbyPopper(const MultiPolygon &mp);
  double ScoreSchwartzberg(const MultiPolygon &mp);
  double ScoreConvexHull  (const MultiPolygon &mp);
  double ScoreReock       (const MultiPolygon &mp);
  GeoCollection ReadShapefile(std::string filename, std::string layername);
  void WriteShapefile(const GeoCollection &mps, std::string filename, std::string layername);
  void CalculateAllScores(GeoCollection &mps);
}

#endif
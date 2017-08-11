#ifndef _compactengine_hpp_
#define _compactengine_hpp_

#include "geom.hpp"

namespace complib {
  double ScorePolsbyPopper(const MultiPolygon &mp);
  double ScoreSchwartzberg(const MultiPolygon &mp);
  double ScoreConvexHull  (const MultiPolygon &mp);
  double ScoreReock       (const MultiPolygon &mp);
}

#endif
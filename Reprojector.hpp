#ifndef _reprojector_hpp_
#define _reprojector_hpp_

#include <string>
#include "geom.hpp"

namespace complib {

class Reprojector {
 public:
  std::string in_projection;
  std::string out_projection;
  void operator()(Point2D &p) const;
};

std::string GetBestEffortProjection(const MultiPolygon &mp);

void ReprojectPoint2D(Point2D &x, const Reprojector &t);
void ReprojectRing(Ring &x, const Reprojector &t);
void ReprojectPolygon(Polygon &x, const Reprojector &t);
void ReprojectMultiPolygon(MultiPolygon &x, const Reprojector &t);

}

#endif

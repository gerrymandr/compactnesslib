#include "Reprojector.hpp"

/*
void Reprojector::operator()(Point2D &p) const {

}


std::string GetBestEffortProjection(const MultiPolygon &mp){
  const double avgx = mp.avgX();
  const double avgy = mp.avgY();

  
}


void ReprojectPoint2D(Point2D &x, const Reprojector &rp){
  rp(x);
}

void ReprojectRing(Ring &x, const Reprojector &rp){
  for(auto &p: x)
    rp(p);
}

void ReprojectPolygon(Polygon &x, const Reprojector &rp){
  ReprojectRing(x.outer, rp);
  for(auto &h: x.holes)
    ReprojectRing(h,rp);
}

void ReprojectMultiPolygon(MultiPolygon &x, const Reprojector &rp){
  for(auto &p: x)
    ReprojectPolygon(p,rp);
}


*/
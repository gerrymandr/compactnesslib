#ifndef _point_hpp_
#define _point_hpp_

class Point3D;

class Point2D {
 public:
  double x;
  double y;
  Point2D() = default;
  Point2D(double x0, double y0);
  Point2D& toRadians();
  Point2D& toDegrees();
  Point2D& rotateTheta(const double rtheta);
};

#endif

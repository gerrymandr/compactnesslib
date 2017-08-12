#include "Point.hpp"
#include <cmath>
#include <cassert>
#include <cstdlib>
#include "doctest.h"
#include <vector>
#include <stdexcept>

const double DEG_TO_RAD = M_PI/180.0;
const double RAD_TO_DEG = 180.0/M_PI;


TEST_CASE("Point2D: Constructor"){
  Point2D p;
}

Point2D::Point2D(double x0, double y0) {
  x = x0;
  y = y0;
}

Point2D& Point2D::toRadians() {
  assert(-180<=x && x<=180);
  assert(-90<=y && y<=90);

  x *= DEG_TO_RAD;
  y *= DEG_TO_RAD;

  assert(-M_PI<=x && x<=M_PI);
  assert(-M_PI/2<=y && y<=M_PI/2);
  return *this;
}

TEST_CASE("Point2D: Conversion to Radians"){
  Point2D p(-93, 45);
  p.toRadians();
  CHECK(p.x==doctest::Approx(-93*DEG_TO_RAD));
  CHECK(p.y==doctest::Approx(45*DEG_TO_RAD));
}



Point2D& Point2D::toDegrees() {
  assert(-M_PI<=x && x<=M_PI);
  assert(-M_PI/2<=y && y<=M_PI/2);

  x *= RAD_TO_DEG;
  y *= RAD_TO_DEG;

  assert(-180<=x && x<=180);
  assert(-90<=y && y<=90);
  return *this;
}

TEST_CASE("Point2D: Conversion to Degrees"){
  Point2D p(-93*DEG_TO_RAD, 45*DEG_TO_RAD);
  p.toDegrees();
  CHECK(p.x==doctest::Approx(-93));
  CHECK(p.y==doctest::Approx(45));
}

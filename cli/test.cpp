#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../compactnesslib.hpp"
#include "../doctest.h"
#include <cmath>

using namespace complib;

TEST_CASE("Square Test"){
  const std::string rect2by2 = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[0,0],[2,0],[2,2],[0,2],[0,0]]]}}]}";

  auto gc = ReadGeoJSON(rect2by2);

  auto &mp = gc.at(0);

  CHECK(areaOfPolygonsIncludingHoles(mp)==4);
  CHECK(perimPolygonOuterRings(mp)==8);
  CHECK(hullAreaPolygonOuterRings(mp)==4);
  CHECK(hullAreaOfHoles(mp)==0);
  CHECK(areaHoles(mp)==0);
  CHECK(perimHoles(mp)==0);
  CHECK(diameterOfEntireMultiPolygon(mp)==doctest::Approx(std::sqrt(2*2+2*2)));
}

TEST_CASE("Circle"){
  MultiPolygon mp;
  mp.emplace_back();             //Make a polygon
  mp.back().emplace_back();      //Make a ring
  auto &ring = mp.back().back(); //Get the ring

  //Make a "circle"
  const double density=1000;
  for(double i=0;i<density;i++)
    ring.emplace_back(33*std::cos(2*M_PI*i/density),33*std::sin(2*M_PI*i/density));
  //Close the "circle"
  ring.emplace_back(33*std::cos(2*M_PI*0),33*std::sin(2*M_PI*0));

  CHECK(areaOfPolygonsIncludingHoles(mp)==doctest::Approx(M_PI*33*33));
  CHECK(perimPolygonOuterRings(mp)==doctest::Approx(2*M_PI*33));
  CHECK(areaHoles(mp)==0);
  CHECK(perimHoles(mp)==0);
  CHECK(diameterOfEntireMultiPolygon(mp)==doctest::Approx(2*33));
  CHECK(ScorePolsbyPopper(mp)==doctest::Approx(1.0));
  CHECK(ScoreSchwartzberg(mp)==doctest::Approx(1.0));
  CHECK(ScoreConvexHull(mp)==doctest::Approx(1.0));
  CHECK(ScoreReock(mp)==doctest::Approx(1.0));
}
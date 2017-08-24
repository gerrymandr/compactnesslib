#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../compactnesslib.hpp"
#include "../doctest.h"
#include <cmath>

using namespace complib;

TEST_CASE("Square Test"){
  const std::string rect2by2 = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[0,0],[2,0],[2,2],[0,2],[0,0]]]}}]}";

  auto gc = ReadGeoJSON(rect2by2);

  auto &mp = gc.at(0);

  SUBCASE("Area"){
    CHECK(areaOfPolygonsIncludingHoles(mp)==4);
    CHECK(perimPolygonOuterRings(mp)==8);
    CHECK(hullAreaPolygonOuterRings(mp)==4);
    CHECK(hullAreaOfHoles(mp)==0);
    CHECK(areaHoles(mp)==0);
    CHECK(perimHoles(mp)==0);
    CHECK(diameterOfEntireMultiPolygon(mp)==doctest::Approx(std::sqrt(2*2+2*2)));
  }
}

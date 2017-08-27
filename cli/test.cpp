#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../compactnesslib.hpp"
#include "../doctest.h"
#include <cmath>
#include <algorithm>

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
    ring.emplace_back(33*std::cos(-2*M_PI*i/density),33*std::sin(-2*M_PI*i/density));
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

TEST_CASE("Name lenth"){
  //Score names can't exceed 10 characters due to shapefile limitations
  for(auto &sn: score_names)
    CHECK(sn.size()<=10);
}

TEST_CASE("Intersection areas"){
  const std::string inita = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[0,0],[4,0],[4,4],[0,4],[0,0]]]}}]}";
  const std::string initb = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[1,1],[3,1],[3,3],[1,3],[1,1]]]}}]}";

  auto gca = ReadGeoJSON(inita);
  auto gcb = ReadGeoJSON(initb);

  SUBCASE("Area forward"){
    CHECK(IntersectionArea(gca[0],gcb[0])==4);
  }

  SUBCASE("Area backward"){
    std::reverse(gca.back().back().begin(),gca.back().back().end());
    std::reverse(gcb.back().back().begin(),gcb.back().back().end());
    CHECK(IntersectionArea(gca[0],gcb[0])==4);
  }
}

TEST_CASE("Intersection areas"){
  const std::string inita = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[0,2],[2,0],[2,2],[0,2],[0,0]]]}}]}";
  const std::string initb = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[1,1],[3,1],[3,3],[1,3],[1,1]]]}}]}";

  auto gca = ReadGeoJSON(inita);
  auto gcb = ReadGeoJSON(initb);

  SUBCASE("Area forward"){
    CHECK(IntersectionArea(gca[0],gcb[0])==1);
  }

  SUBCASE("Area backward"){
    std::reverse(gca.back().back().begin(),gca.back().back().end());
    std::reverse(gcb.back().back().begin(),gcb.back().back().end());
    CHECK(IntersectionArea(gca[0],gcb[0])==4);
  }
}







/*
TEST_CASE("Polygon"){
  Polygon p;
  Point2D a(-93,45);
  Point2D b(-93,50);
  Point2D c(-90,50);
  Point2D d(-90,45);
  p.outer.push_back(a);
  p.outer.push_back(b);
  p.outer.push_back(c);
  p.outer.push_back(d);
  p.toRadians();
  b.toRadians();
  CHECK(p.exterior[1].x==doctest::Approx(b.x));
  CHECK(p.exterior[1].y==doctest::Approx(b.y));
  p.toDegrees();
  CHECK(p.exterior[2].x==doctest::Approx(c.x));
  CHECK(p.exterior[2].y==doctest::Approx(c.y));
  p.containsPoint(Point2D(-92,47));
  p.containsPoint(Point2D(-94,47));
}*/
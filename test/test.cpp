#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../compactnesslib.hpp"
#include "../lib/doctest.h"
#include <cmath>
#include <map>
#include <iostream>
#include <algorithm>

using namespace complib;

TEST_CASE("Data test"){
  auto gc = complib::ReadShapefile("test_data/cb_2015_us_cd114_20m.shp");
  for(const auto &mp: gc)
    CHECK(areaExcludingHoles(mp)>0);
  for(const auto &mp: gc)
    CHECK(IntersectionArea(mp,mp.getHull())>0);
  CHECK(gc.size()==216);
  for(const auto &mp: gc){
    if(mp.props.at("GEOID")=="1307"){
      CHECK(perimIncludingHoles(mp)==doctest::Approx(171105));
      CHECK(areaIncludingHoles(mp)==doctest::Approx(1048787000));
    } else if(mp.props.at("GEOID")=="1226"){
      CHECK(perimIncludingHoles(mp)==doctest::Approx(776915));
      CHECK(areaIncludingHoles(mp)==doctest::Approx(7659883000));
    }
  }
}

TEST_CASE("Square Test"){
  const std::string rect2by2 = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[0,0],[2,0],[2,2],[0,2],[0,0]]]}}]}";

  auto gc = ReadGeoJSON(rect2by2);

  auto &mp = gc.at(0);

  SUBCASE("Simple metrics"){
    CHECK(areaIncludingHoles(mp)==4);
    CHECK(perimExcludingHoles(mp)==8);
    CHECK(hullAreaPolygonOuterRings(mp)==4);
    CHECK(hullAreaOfHoles(mp)==0);
    CHECK(areaHoles(mp)==0);
    CHECK(perimHoles(mp)==0);
    CHECK(diameterOfEntireMultiPolygon(mp)==doctest::Approx(std::sqrt(2*2+2*2)));
  }

  SUBCASE("Bounding Circle"){
    const auto circle = GetBoundingCircle(mp);
    CHECK(areaIncludingHoles(circle)==doctest::Approx(2*M_PI));
  }
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

  CHECK(areaIncludingHoles(mp)==doctest::Approx(M_PI*33*33));
  CHECK(perimExcludingHoles(mp)==doctest::Approx(2*M_PI*33));
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
  for(auto &sn: getListOfUnboundedScores())
    CHECK(sn.size()<=10);
  for(auto &sn: getListOfBoundedScores())
    CHECK(sn.size()<=10);
}

TEST_CASE("Name Matching"){
  for(auto &sn: getListOfUnboundedScores())
    CHECK(unbounded_score_map.count(sn));
  for(auto &sn: getListOfBoundedScores())
    CHECK(bounded_score_map.count(sn));
}

TEST_CASE("Name uniqueness"){
  std::map<std::string, int> score_names;
  for(auto &sn: getListOfUnboundedScores())
    score_names[sn]++;
  for(auto &sn: getListOfBoundedScores())
    score_names[sn]++;
  for(const auto &kv: score_names)
    CHECK(kv.second==1);
}

TEST_CASE("Intersection area: big square and little square"){
  const std::string inita = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[0,0],[4,0],[4,4],[0,4],[0,0]]]}}]}";
  const std::string initb = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[1,1],[3,1],[3,3],[1,3],[1,1]]]}}]}";

  auto gca = ReadGeoJSON(inita);
  auto gcb = ReadGeoJSON(initb);

  SUBCASE("Area forward"){
    CHECK(IntersectionArea(gca[0],gcb[0])==4);
  }

  SUBCASE("Area backward"){
    gca.reverse();
    gcb.reverse();
    CHECK(IntersectionArea(gca[0],gcb[0])==4);
  }
}

TEST_CASE("Intersection areas big square with hole and little square"){
  const std::string inita = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[0,2],[2,0],[2,2],[0,2],[0,0]]]}}]}";
  const std::string initb = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[1,1],[3,1],[3,3],[1,3],[1,1]]]}}]}";

  auto gca = ReadGeoJSON(inita);
  auto gcb = ReadGeoJSON(initb);

  SUBCASE("Area forward"){
    CHECK(IntersectionArea(gca[0],gcb[0])==1);
  }

  //This should give the same answer as the foregoing since the underlying
  //clipper library will orientate things correctly for itself
  SUBCASE("Area backward"){
    gca.reverse();
    gcb.reverse();
    CHECK(IntersectionArea(gca[0],gcb[0])==1);
  }
}


TEST_CASE("Polygon with hole"){
  const std::string inita = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[0,0],[4,0],[4,4],[0,4],[0,0]],[[1,1],[2,1],[2,2],[1,2],[1,1]]]}}]}";
  const std::string initb = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[1,1],[3,1],[3,3],[1,3],[1,1]]]}}]}";

  auto gca = ReadGeoJSON(inita);
  auto gcb = ReadGeoJSON(initb);

  CHECK(areaIncludingHoles(gca[0])==16);
  CHECK(areaIncludingHoles(gcb[0])==4);
  CHECK(areaHoles(gca[0])==1);
  CHECK(IntersectionArea(gca[0],gcb[0])==3);
}

TEST_CASE("WKT output"){
  const std::string inita = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[0,0],[4,0],[4,4],[0,4],[0,0]],[[1,1],[2,1],[2,2],[1,2],[1,1]]]}}]}";
  const auto gca = ReadGeoJSON(inita);
  std::cout<<GetWKT(gca.at(0))<<std::endl;
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
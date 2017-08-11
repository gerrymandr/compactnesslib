#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "compactengine.hpp"
#include "geojson.hpp"

#include <iostream>
#include <fstream>
#include <streambuf>
#include <stdexcept>

namespace rj = rapidjson;

namespace complib {

template<class T>
Point2D ParsePoint(const T &p){
  auto v = p.GetArray();
  Point2D temp;
  temp.x = p[0].GetDouble();
  temp.y = p[1].GetDouble();
  return temp;
}

template<class T>
Ring ParseRing(const T &r){
  Ring temp;
  for(const auto &c: r.GetArray())
    temp.push_back(ParsePoint(c));
  return temp;
}


//rj::StringBuffer sb;
//rj::Writer<rj::StringBuffer> writer(sb);
//d["properties"].Accept(writer);
//std::string s = sb.GetString();

template<class T>
void PrintMembers(const T &d){
  std::cerr<<"Members: ";
  for (rj::Value::ConstMemberIterator itr = d.MemberBegin(); itr != d.MemberEnd(); ++itr)
    std::cerr<< itr->name.GetString() <<" ";
  std::cerr<<std::endl;
}

template<class T>
Polygon ParsePolygon(const T &coor){
  Polygon mp;
  mp.outer = ParseRing(coor[0]);
  for(rj::SizeType i=1;i<coor.Size();i++)
    mp.holes.push_back(ParseRing(coor[i]));

  return mp;
}

template<class T>
std::map<std::string, std::any> GetProperties(const T &d){
  std::map<std::string, std::any> props;
  for (rj::Value::ConstMemberIterator itr = d.MemberBegin(); itr != d.MemberEnd(); ++itr){
    std::any prop;
    switch(itr->value.GetType()){
      case 0:
        prop = "NULL";break;
      case 1:
        prop = "False";break;
      case 2:
        prop = "True";break;
      case 3:
        throw std::runtime_error("Object cannot be a property, yet!");
      case 4:
        throw std::runtime_error("Array cannot be a property, yet!");
      case 5:
        prop = itr->value.GetString();break;
      case 6:
        prop = itr->value.GetDouble();break;
      default:
        throw std::runtime_error("Unrecognized value!");
    }
    props[itr->name.GetString()] = prop;
  }
}

template<class T>
MultiPolygon ParseTopPolygon(const T &d){
  if(!d.HasMember("geometry"))
    throw std::runtime_error("Polygon has no 'geometry' member!");
  if(!d["geometry"].IsObject())
    throw std::runtime_error("Geometry not an object!");
  const auto geom = d["geometry"].GetObject();
  if(!geom.HasMember("coordinates"))
    throw std::runtime_error("Could not get coordinates!");
  const auto coor = geom["coordinates"].GetArray();

  MultiPolygon mps;
  mps.push_back(ParsePolygon(coor));
  return mps;
}

template<class T>
MultiPolygon ParseMultiPolygon(const T &d){
  MultiPolygon mps;

  if(!d.HasMember("geometry"))
    throw std::runtime_error("Polygon has no 'geometry' member!");
  if(!d["geometry"].IsObject())
    throw std::runtime_error("Geometry not an object!");
  const auto geom = d["geometry"].GetObject();
  if(!geom.HasMember("coordinates"))
    throw std::runtime_error("Could not get coordinates!");
  const auto coor = geom["coordinates"].GetArray();
  for(const auto &poly: coor)
    mps.emplace_back(ParsePolygon(poly.GetArray()));

  return mps;
}

template<class T>
MultiPolygon ParseFeature(const T &d){
  MultiPolygon mps;

  const std::string geotype = d["geometry"]["type"].GetString();
  if(geotype=="MultiPolygon")
    return ParseMultiPolygon(d);
  else if(geotype=="Polygon")
    return ParseTopPolygon(d);
  else
    throw std::runtime_error("Unexpected data type - skipping!");
}

GeoCollection ReadGeoJSON(std::string geojson){

  // 1. Parse a JSON string into DOM.
  rj::Document d;

  d.Parse(geojson.c_str());

  if(!d.IsObject())
    throw std::runtime_error("GeoJSON not an object!");

  if(!d.HasMember("type"))
    throw std::runtime_error("No type property!");
  if(!d["type"].IsString())
    throw std::runtime_error("Type not a string!");
  if(d["type"]!="FeatureCollection")
    throw std::runtime_error("Not a FeatureCollection!");

  GeoCollection mps;
  for(const auto &f: d["features"].GetArray())
    mps.push_back(ParseFeature(f.GetObject()));

  // // 3. Stringify the DOM
  // rj::StringBuffer buffer;
  // rj::Writer<rj::StringBuffer> writer(buffer);
  // d.Accept(writer);
  // // Output {"project":"rapidjson","stars":11}
  // std::cout << buffer.GetString() << std::endl;

  std::cout<<"Found features = "<<mps.size()<<std::endl;

  return mps;
}

GeoCollection ReadGeoJSONFile(std::string filename){
  std::ifstream fin(filename);

  std::string geojson((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());

  return ReadGeoJSON(geojson);
}

}
#include "json.hpp"
#include "compactengine.hpp"
#include "geojson.hpp"
#include <iostream>
#include <fstream>
#include <streambuf>
#include <stdexcept>
#include <sstream>

using json = nlohmann::json;

namespace complib {

//Use to deduce types
//template<typename T> struct TD;
//e.g. TD<decltype(WHAT_AM_I_VAR_NAME)> td;

template<class T>
Ring ParseRing(const T &r){
  Ring temp;
  for(const auto &c: r)
    temp.emplace_back(c[0],c[1]);
  return temp;
}


//rj::StringBuffer sb;
//rj::Writer<rj::StringBuffer> writer(sb);
//d["properties"].Accept(writer);
//std::string s = sb.GetString();

template<class T>
void PrintMembers(const T &d){
  std::cerr<<d<<std::endl;
}

template<class T>
Polygon ParsePolygon(const T &coor){
  Polygon mp;
  //First ring is the outer ring, all the others are holes
  for(const auto &r: coor)
    mp.push_back(ParseRing(r));

  return mp;
}

template<class T>
Props GetProperties(const T &d){
  Props props = d;
  return props;
}

template<class T>
const T& GetToCoordinates(const T &d){
  if(d.count("geometry")){
    return d.at("geometry").at("coordinates");
  } else if(d.count("coordinates")){
    return d.at("coordinates");
  } else {
    throw std::runtime_error("Could find neither geometry nor coordinates!");
  }
}

template<class T>
MultiPolygon ParseTopPolygon(const T &d){
  MultiPolygon mps;
  mps.push_back(ParsePolygon(GetToCoordinates(d)));
  return mps;
}

template<class T>
MultiPolygon ParseMultiPolygon(const T &d){
  MultiPolygon mps;
  for(const auto &poly: GetToCoordinates(d))
    mps.emplace_back(ParsePolygon(poly));

  return mps;
}

template<class T>
MultiPolygon ParseFeature(const T &d){
  MultiPolygon mps;

  const std::string geotype = d["geometry"]["type"];
  if(geotype=="MultiPolygon")
    return ParseMultiPolygon(d);
  else if(geotype=="Polygon")
    return ParseTopPolygon(d);
  else
    throw std::runtime_error("Unexpected data type - skipping!");
}

GeoCollection ReadGeoJSON(const std::string geojson){
  GeoCollection mps;
  auto d = json::parse(geojson);

  if(!d.is_object())
    throw std::runtime_error("GeoJSON not an object!");

  if(!d.count("type"))
    throw std::runtime_error("No type property!");
  if(!d["type"].is_string())
    throw std::runtime_error("Type not a string!");

  if(d["type"]=="MultiPolygon"){
    mps.push_back(ParseMultiPolygon(d));
  } else if(d["type"]=="Polygon"){
    mps.push_back(ParseTopPolygon(d));
  } else if(d["type"]=="FeatureCollection"){
    for(const auto &f: d["features"])
      mps.push_back(ParseFeature(f));
  } else {
    throw std::runtime_error("Not a FeatureCollection or MultiPolygon or Polygon!");
  }

  return mps;
}

GeoCollection ReadGeoJSONFile(std::string filename){
  std::ifstream fin(filename);

  std::string geojson((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());

  return ReadGeoJSON(geojson);
}

std::string OutScoreJSON(const GeoCollection &gc, const std::string id){
  std::ostringstream oss;

  const bool use_id = !id.empty();

  oss<<"{\n";
  for(unsigned int i=0;i<gc.size();i++){
    oss<<"\t\"";
    if(use_id)
      oss<<gc[i].props.at(id);
    else
      oss<<i;
    oss<<"\":{\n";

    int inserted = 0;
    for(unsigned int sn=0;sn<score_names.size();sn++){
      if(gc[i].scores.count(score_names[sn])){
        oss<<"\t\t\""<<score_names[sn]<<"\":"<<gc[i].scores.at(score_names[sn]);
        if(inserted++<gc[i].scores.size()-1)
          oss<<",\n";
      }
    }

    oss<<"\n\t}";
    if(i<gc.size()-1)
      oss<<",\n";
  }
  oss<<"\n}";

  return oss.str();
}

}
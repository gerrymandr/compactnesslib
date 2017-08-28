#include "wkt.hpp"
#include "geom.hpp"
#include <iostream>
#include <fstream>
#include <streambuf>
#include <stdexcept>
#include <sstream>
#include <string>

namespace complib {



void TrimStr(size_t &start, const std::string &str){
  for(;start<str.size();start++)
    if(str[start]!=' ')
      return;
  throw std::runtime_error("WKT ended with whitespace - expected something!");
}



double GetNumber(size_t &start, const std::string &wktstr){
  TrimStr(start,wktstr);
  size_t numend = start;
  for(;numend<wktstr.size();numend++)
    if(wktstr[numend]==' ' || wktstr[numend]==',' || wktstr[numend]==')')
      break;
  const double ret = std::stod(wktstr.substr(start,numend-start));
  start = numend;
  return ret;
}



Point2D ParsePoint(size_t &start, const std::string &wktstr){
  Point2D temp;
  temp.x = GetNumber(start,wktstr);
  temp.y = GetNumber(start,wktstr);
  for(;start<wktstr.size();start++)
    if(wktstr[start]==',' || wktstr[start]==')')
      return temp;
  throw std::runtime_error("Could not find the end of the point!");
}



Ring ParseRing(size_t &start, const std::string &wktstr){
  Ring ring;

  if(wktstr[start]!='(')
    throw std::runtime_error(std::string("Ring: Expected '(' found '")+wktstr[start]+"'!");
  TrimStr(++start, wktstr);

  while(true){
    TrimStr(start,wktstr);
    if(wktstr[start]==','){
      TrimStr(++start, wktstr); //Move forward from ',' and skip whitespace
    } else if(wktstr[start]==')'){
      start++;
      break;
    } else {
      ring.v.emplace_back(ParsePoint(start,wktstr));
    }
  }

  return ring;
}



Polygon ParsePolygon(size_t &start, const std::string &wktstr){
  Polygon poly;
  if(wktstr[start]!='(')
    throw std::runtime_error(std::string("Polygon: Expected '(' found '")+wktstr[start]+"'!");
  TrimStr(++start, wktstr);

  while(true){
    TrimStr(start,wktstr);
    if(wktstr[start]=='('){
      //First ring is the outer ring, all the others are holes
      poly.v.emplace_back(ParseRing(start,wktstr));
    } else if(wktstr[start]==','){
      TrimStr(++start, wktstr); //Move forward from ',' and skip whitespace
    } else if(wktstr[start]==')'){
      start++;
      break;
    }
  }

  return poly;
}



MultiPolygon ParseTopPolygon(size_t start, const std::string &wktstr){
  MultiPolygon mps;
  TrimStr(start,wktstr);
  mps.v.push_back(ParsePolygon(start,wktstr));
  return mps;
}



MultiPolygon ParseMultiPolygon(size_t start, const std::string &wktstr){
  MultiPolygon mps;
  TrimStr(start,wktstr);
  if(wktstr[start]!='(')
    throw std::runtime_error(std::string("MultiPolygon: Expected '(' found '")+wktstr[start]+"'!");
  TrimStr(++start, wktstr); //Move forward from '(' and skip whitespace

  while(true){
    TrimStr(start,wktstr);
    if(wktstr[start]=='('){
      mps.v.emplace_back(ParsePolygon(start,wktstr));
    } else if(wktstr[start]==','){
      TrimStr(++start, wktstr); //Move forward from ',' and skip whitespace
    } else if(wktstr[start]==')'){
      start++;
      break;
    }
  }

  return mps;
}



GeoCollection ReadWKT(std::string wktstr){
  GeoCollection gc;

  size_t start = 0;
  TrimStr(start,wktstr);

  //TODO: Trim string
  if(wktstr.compare(start,12,"MULTIPOLYGON")==0){
    gc.v.push_back(ParseMultiPolygon(start+12,wktstr));
  } else if(wktstr.compare(start,7,"POLYGON")==0){
    gc.v.push_back(ParseTopPolygon(start+7,wktstr));
  } else{
    throw std::runtime_error("Unrecognized geometry!");
  }

  gc.correctWindingDirection();
  
  return gc;
}

GeoCollection ReadWKTFile(std::string filename) {
  std::ifstream fin(filename);

  std::string wktstr((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());

  return ReadWKT(wktstr);
}

}
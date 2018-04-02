#include "wkt.hpp"
#include "geom.hpp"
#include <iostream>
#include <fstream>
#include <locale>         // std::locale, std::isspace
#include <streambuf>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <string>


namespace complib {



void TrimStr(size_t &start, const std::string &str){
  std::locale loc;

  for(;start<str.size();start++)
    if(!std::isspace(str[start],loc))
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
      ring.emplace_back(ParsePoint(start,wktstr));
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
      poly.emplace_back(ParseRing(start,wktstr));
    } else if(wktstr[start]==','){
      TrimStr(++start, wktstr); //Move forward from ',' and skip whitespace
    } else if(wktstr[start]==')'){
      start++;
      break;
    }
  }

  return poly;
}



MultiPolygon ParseTopPolygon(size_t &start, const std::string &wktstr){
  MultiPolygon mps;
  TrimStr(start,wktstr);
  mps.push_back(ParsePolygon(start,wktstr));
  return mps;
}



MultiPolygon ParseMultiPolygon(size_t &start, const std::string &wktstr){
  MultiPolygon mps;
  TrimStr(start,wktstr);
  if(wktstr[start]!='(')
    throw std::runtime_error(std::string("MultiPolygon: Expected '(' found '")+wktstr[start]+"'!");
  TrimStr(++start, wktstr); //Move forward from '(' and skip whitespace

  while(true){
    TrimStr(start,wktstr);
    if(wktstr[start]=='('){
      mps.emplace_back(ParsePolygon(start,wktstr));
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

  while(true){
    if(wktstr.compare(start,12,"MULTIPOLYGON")==0){
      start += 12; //Skip MULTIPOLYGON
      gc.push_back(ParseMultiPolygon(start,wktstr));
    } else if(wktstr.compare(start,7,"POLYGON")==0){
      start += 7; //Skip POLYGON
      gc.push_back(ParseTopPolygon(start,wktstr));
    } else{
      throw std::runtime_error("Unrecognized geometry!");
    }

    try {
      TrimStr(start,wktstr);
    } catch (...){
      break; //Okay, we've reached the end of the input data
    }
  }

  gc.correctWindingDirection();
  
  return gc;
}

GeoCollection ReadWKTFile(std::string filename) {
  std::ifstream fin(filename);

  std::string wktstr((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());

  return ReadWKT(wktstr);
}



std::string GetWKT(const MultiPolygon &mp){
  std::ostringstream ret;
  ret<<"MULTIPOLYGON (";
  for(unsigned int p=0;p<mp.size();p++){
    const auto &poly = mp.at(p);
    ret<<"(";
    for(unsigned int r=0;r<poly.size();r++){
      const auto &ring = poly.at(r);
      ret<<"(";
      for(unsigned int i=0;i<ring.size();i++){
        const auto &pt = ring.at(i);
        ret<<std::fixed<<std::setprecision(10)<<ring[i].x<<" "
           <<std::fixed<<std::setprecision(10)<<ring[i].y;
        if(i<ring.size()-1)
          ret<<",";
      }
      ret<<")";
      if(r<poly.size()-1)
        ret<<",";
    }
    ret<<")";
    if(p<mp.size()-1)
      ret<<",";
  }
  ret<<")";

  return ret.str();
}

}
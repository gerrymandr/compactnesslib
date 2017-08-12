#include "utility.hpp"
#include <ostream>
#include <string>
#include <stdexcept>
#include <iomanip>

namespace complib {

const int GEOJSON_PRECISION = 10;

std::ostream& operator<<(std::ostream &out, const std::any &propval){
  if(typeid(int) == propval.type()) 
    out<<std::any_cast<int>(propval);
  else if(typeid(long) == propval.type()) 
    out<<std::any_cast<long>(propval);
  else if(typeid(long long) == propval.type()) 
    out<<std::any_cast<long long>(propval);
  else if(typeid(double) == propval.type()) 
    out<<std::fixed<<std::setprecision(GEOJSON_PRECISION)<<std::any_cast<double>(propval);
  else if(typeid(std::string) == propval.type()) 
    out<<"\""<<std::any_cast<std::string>(propval)<<"\"";
  else 
    throw std::runtime_error("Unrecognized property type '" + std::string(propval.type().name()) + "'!");
  return out;
}

}
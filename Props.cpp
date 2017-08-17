#include "Props.hpp"
#include <iostream>

namespace complib {
  void PrintProps(const Props &ps){
    for(const auto &p: ps){
      // if(typeid(int) == p.second.type()) {
      //   std::cout<<p.first<<" = "<<std::any_cast<int>(p.second)<<std::endl;
      // } else if(typeid(long) == p.second.type()) {
      //   std::cout<<p.first<<" = "<<std::any_cast<long>(p.second)<<std::endl;
      // } else if(typeid(double) == p.second.type()) {
      //   std::cout<<p.first<<" = "<<std::any_cast<double>(p.second)<<std::endl;
      // } else if(typeid(std::string) == p.second.type()) {
      //   std::cout<<p.first<<" = "<<std::any_cast<std::string>(p.second).c_str()<<std::endl;
      // } else {
      //   std::cerr<<"Unrecognized property type '"<<p.second.type().name()<<"'!"<<std::endl;
      // }
      std::cout<<p.first<<" = "<<p.second<<std::endl;
    }
  }
}


// std::ostream& operator<<(std::ostream &out, const std::any &propval){
//   if(typeid(int) == propval.type()) 
//     out<<std::any_cast<int>(propval);
//   else if(typeid(long) == propval.type()) 
//     out<<std::any_cast<long>(propval);
//   else if(typeid(long long) == propval.type()) 
//     out<<std::any_cast<long long>(propval);
//   else if(typeid(double) == propval.type()) 
//     out<<std::fixed<<std::setprecision(GEOJSON_PRECISION)<<std::any_cast<double>(propval);
//   else if(typeid(std::string) == propval.type()) 
//     out<<"\""<<std::any_cast<std::string>(propval)<<"\"";
//   else 
//     throw std::runtime_error("Unrecognized property type '" + std::string(propval.type().name()) + "'!");
//   return out;
// }

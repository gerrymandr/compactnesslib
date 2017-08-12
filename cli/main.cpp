#include "engine/compactengine.hpp"
#include <iostream>
#include <fstream>

std::ostream& operator<<(std::ostream &out, const std::any &pval){
  if(typeid(int) == pval.type()) {
    out<<std::any_cast<int>(pval);
  } else if(typeid(long) == pval.type() || typeid(long long) == pval.type()) {
    out<<std::any_cast<long long>(pval);
  } else if(typeid(double) == pval.type()) {
    out<<std::any_cast<double>(pval);
  } else if(typeid(std::string) == pval.type()) {
    out<<std::any_cast<std::string>(pval).c_str();
  } else {
    throw std::runtime_error("Unrecognized field type (value write) '" + std::string(pval.type().name()) + "'!");
  }      
  return out;
}

int main(int argc, char **argv) {
  if(argc!=2){
    std::cerr<<"Syntax: "<<argv[0]<<" <INPUT>"<<std::endl;
    return -1;
  }

  auto mpolys = complib::ReadGeoJSONFile(argv[1]);

  complib::CalculateAllScores(mpolys);

  std::cout<<OutScoreJSON(mpolys, "")<<std::endl;

  return 0;
}
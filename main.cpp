#include "shapefile.hpp"
#include "compactengine.hpp"
#include <iostream>

int main(int argc, char **argv){
  if(argc!=3){
    std::cout<<"Syntax: "<<argv[0]<<" <DSN> <LAYER>"<<std::endl;
    return -1;
  }

  auto mpolys = ReadShapefile(argv[1], argv[2]);

  for(const auto &mp: mpolys){
    //mp.print();
    complib::ScorePolsbyPopper(mp);
    complib::ScoreSchwartzberg(mp);
    complib::ScoreConvexHull  (mp);
    complib::ScoreReock       (mp);
  }
}
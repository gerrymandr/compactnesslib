#include "../compactnesslib.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char **argv) {
  if(argc!=2){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input File>"<<std::endl;
    return -1;
  }

  std::string filename = argv[1];

  complib::GeoCollection gc;

  if(filename.find(".geojson")!=std::string::npos)
    gc = complib::ReadGeoJSONFile(filename);
  else if(filename.find(".shp")!=std::string::npos)
    gc = complib::ReadShapefile(filename);
  else
    throw std::runtime_error("Unrecognized file extension! Can use '.geojson' or '.shp'.");

  complib::CalculateAllScores(gc);

  std::cout<<OutScoreJSON(gc, "")<<std::endl;

  WriteShapefile(gc,"/z/outshape");

  return 0;
}
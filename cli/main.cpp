#include "../compactnesslib.hpp"
#include "Timer.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char **argv) {
  if(argc!=3){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input File> <Output File>"<<std::endl;
    std::cerr<<"\tIf <Output File> = '-' then GeoJSON is printed to stdout."<<std::endl;
    std::cerr<<"\tIf <Output File> = 'augment', the input shapefile is augmented to include scores."<<std::endl;
    return -1;
  }

  std::string in_filename  = argv[1];
  std::string out_filename = argv[2];

  std::cout<<"Processing '"<<in_filename<<"'..."<<std::endl;

  complib::GeoCollection gc;

  {
    Timer tmr;
    std::cout<<"Loading..."<<std::endl;
    if(in_filename.find(".geojson")!=std::string::npos)
      gc = complib::ReadGeoJSONFile(in_filename);
    else if(in_filename.find(".shp")!=std::string::npos)
      gc = complib::ReadShapefile(in_filename);
    else
      throw std::runtime_error("Unrecognized input file extension! Can use '.geojson' or '.shp'.");
    std::cout<<"Finished in = "<<tmr.elapsed()<<" s"<<std::endl;
  }

  {
    Timer tmr;
    std::cout<<"Scoring..."<<std::endl;
    complib::CalculateAllScores(gc);
    std::cout<<"Finished in = "<<tmr.elapsed()<<" s"<<std::endl;
  }

  {
    Timer tmr;
    std::cout<<"Writing..."<<std::endl;
    if(out_filename=="-"){
      std::cout<<OutScoreJSON(gc, "")<<std::endl;
    } else if(out_filename=="augment"){
      WriteShapeScores(gc, in_filename);
    } else if(out_filename.find(".geojson")!=std::string::npos){
      std::ofstream fout(out_filename);
      fout<<complib::OutScoreJSON(gc,"");
    } else if(out_filename.find(".shp")!=std::string::npos){
      WriteShapefile(gc, out_filename);
    } else if(out_filename.find(".csv")!=std::string::npos){
      std::ofstream fout(out_filename);
      fout<<complib::OutScoreCSV(gc,"");
    } else {
      throw std::runtime_error("Unrecognized output file directive! Can use '*.geojson' or '*.shp' or '-' or 'augment' .");
    }
    std::cout<<"Finished in = "<<tmr.elapsed()<<" s"<<std::endl;
  }

  return 0;
}
#include "../compactnesslib.hpp"
#include "Timer.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char **argv) {
  if(argc!=4){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input File> <Output File> <ID>"<<std::endl;
    std::cerr<<"\tIf <Output File> = '-json' then GeoJSON is printed to stdout."<<std::endl;
    std::cerr<<"\tIf <Output File> = '-csv' then CSV is printed to stdout."<<std::endl;
    std::cerr<<"\tIf <Output File> = 'augment', the input shapefile is augmented to include scores."<<std::endl;
    std::cerr<<"\tIf <ID>          = The attribute to which each output should be keyed."<<std::endl;
    return -1;
  }

  std::string in_filename  = argv[1];
  std::string out_filename = argv[2];
  std::string id           = argv[3];

  std::cout<<"Processing '"<<in_filename<<"'..."<<std::endl;

  complib::GeoCollection gc;

  {
    Timer tmr;
    std::cout<<"Loading..."<<std::endl;
    if(in_filename.compare(in_filename.size()-8,8,".geojson")==0)
      gc = complib::ReadGeoJSONFile(in_filename);
    else if(in_filename.compare(in_filename.size()-4,4,".shp")==0)
      gc = complib::ReadShapefile(in_filename);
    else if(in_filename.compare(in_filename.size()-4,4,".wkt")==0)
      gc = complib::ReadWKTFile(in_filename);
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
    if(out_filename=="-json"){
      std::cout<<OutScoreJSON(gc, id)<<std::endl;
    } else if(out_filename=="-csv"){
      std::cout<<OutScoreCSV(gc, id)<<std::endl;      
    } else if(out_filename=="augment"){
      WriteShapeScores(gc, in_filename);
    } else if(out_filename.compare(out_filename.size()-8,8,".geojson")==0){
      std::ofstream fout(out_filename);
      fout<<complib::OutScoreJSON(gc,id);
    } else if(out_filename.compare(out_filename.size()-4,4,".shp")==0){
      WriteShapefile(gc, out_filename);
    } else if(out_filename.compare(out_filename.size()-4,4,".csv")==0){
      std::ofstream fout(out_filename);
      fout<<complib::OutScoreCSV(gc,id);
    } else {
      throw std::runtime_error("Unrecognized output file directive! Can use '*.geojson' or '*.shp' or '-' or 'augment' .");
    }
    std::cout<<"Finished in = "<<tmr.elapsed()<<" s"<<std::endl;
  }

  return 0;
}
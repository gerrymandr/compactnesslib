#include "shapefile.hpp"
#include "compactengine.hpp"
#include <iostream>

//http://www.electiondataarchive.org/datacenter-gred.html

//ogr2ogr -f 'ESRI Shapefile' -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs' tl_2016_us_cd115.shp /z/conv.shp

//./dgfinder.exe /z/conv.shp conv /z/out

int main(int argc, char **argv){
  if(argc!=4){
    std::cout<<"Syntax: "<<argv[0]<<" <DSN> <LAYER> <OutName>"<<std::endl;
    return -1;
  }

  std::string dsn       = argv[1];
  std::string layername = argv[2];
  std::string outname   = argv[3];

  auto mpolys = ReadShapefile(dsn, layername);

  #pragma omp parallel for
  for(unsigned int i=0;i<mpolys.size();i++){
    mpolys[i].props["perim"]      = mpolys[i].perim();
    mpolys[i].props["area"]       = mpolys[i].area();
    mpolys[i].props["PolsbyPopp"] = complib::ScorePolsbyPopper(mpolys[i]);
    mpolys[i].props["Schwartzbe"] = complib::ScoreSchwartzberg(mpolys[i]);
    mpolys[i].props["ConvexHull"] = complib::ScoreConvexHull  (mpolys[i]);
    mpolys[i].props["Reock"]      = complib::ScoreReock       (mpolys[i]);
  }

  WriteShapefile(mpolys, outname, layername);
}
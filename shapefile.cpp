#include "shapefile.hpp"
#include "shapelib/shapefil.h"
#include <iostream>
#include <stdexcept>
#include <memory>
#include <map>
#include <algorithm>
#include <set>
#include <fstream>
#include <sstream>
#include "geom.hpp"

namespace complib {

static bool IsHole(const SHPObject *const psCShape, const int ringi){
  const int first_vtx = psCShape->panPartStart[ringi];

  int last_vtx;
  if(ringi>=psCShape->nParts-1)                   //Last ring in the set
    last_vtx = psCShape->nVertices;               //Choose last vertex
  else                                            //Not last right
    last_vtx = psCShape->panPartStart[ringi + 1]; //Last vertex is first of next ring

  const int size = last_vtx-first_vtx;

  const double *const x = psCShape->padfX + first_vtx;
  const double *const y = psCShape->padfY + first_vtx;

  double area = 0;

  //The "shoelace" algorithm
  int j = size-1;
  for(int i=0;i<size;i++){
    area += (x[j]*y[i])-(x[i]*y[j]);
    j = i;
  }

  return area>0;
}












static void ReadShapeAttributes(GeoCollection &gc, std::string filename){
  DBFHandle hDBF = DBFOpen( filename.c_str(), "rb" );
  if( hDBF == NULL )
    throw std::runtime_error("Failed to open file '"+filename+"'!");

  if( DBFGetFieldCount(hDBF) == 0 )
    throw std::runtime_error("No fields in the table file!");

  for(int iRecord = 0; iRecord < DBFGetRecordCount(hDBF); iRecord++ ){        
    for(int i = 0; i < DBFGetFieldCount(hDBF); i++ ){
      char szTitle[20]; //Code example showed 12, I extend to 20 for... safety?
      //char szFormat[32];
      int  nWidth;
      int  nDecimals;

      //TODO: Inefficient? We could do this once for each field above and then
      //refer to a vector of the stored field info.
      const DBFFieldType eType = DBFGetFieldInfo(hDBF, i, szTitle, &nWidth, &nDecimals);
      (void)eType;
            
      //USE IF YOU WISH TO DECODE ATTRIBUTES
      // if( DBFIsAttributeNULL( hDBF, iRecord, i ) ){
      //   if( eType == FTString )
      //     sprintf( szFormat, "%%-%ds", nWidth );
      //   else
      //     sprintf( szFormat, "%%%ds", nWidth );

      //   printf( szFormat, "(NULL)" );
      // } else {
      //   switch( eType ) {
      //     case FTString:
      //       sprintf( szFormat, "%%-%ds", nWidth );
      //       printf( szFormat, 
      //               DBFReadStringAttribute( hDBF, iRecord, i ) );
      //       break;
            
      //     case FTInteger:
      //       sprintf( szFormat, "%%%dd", nWidth );
      //       printf( szFormat, 
      //               DBFReadIntegerAttribute( hDBF, iRecord, i ) );
      //       break;
            
      //     case FTDouble:
      //       sprintf( szFormat, "%%%d.%dlf", nWidth, nDecimals );
      //       printf( szFormat, 
      //               DBFReadDoubleAttribute( hDBF, iRecord, i ) );
      //       break;
            
      //     default:
      //       break;
      //   }
      // }


      //GRAB ATTRIBUTES WITHOUT DECODING THEM INTO USABLE DATA
      const auto attrib = DBFReadStringAttribute( hDBF, iRecord, i );
      gc.at(iRecord).props[szTitle] = attrib;
    }
  }

  DBFClose( hDBF );
}






//TODO: Do we need to worry about layers?
void ReadShapes(GeoCollection &mgons, std::string filename){
  int nShapeType;
  int nEntities;
  double adfMinBound[4];
  double adfMaxBound[4];

  SHPHandle hSHP = SHPOpen(filename.c_str(), "rb");
  if( hSHP == NULL )
    throw std::runtime_error("Failed to open shapefile '" + filename + "'!");

  SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

  if(nShapeType!=SHPT_POLYGON)
    throw std::runtime_error("Can only work with polygon shapefiles!");

  SHPObject *psShape;
  for(int i=0;i<nEntities;i++){
    mgons.emplace_back();     //Add a new MultiPolygon 
    auto &mp = mgons.back();  //Get a reference to it

    psShape = SHPReadObject( hSHP, i );
    if(psShape==NULL)
      throw std::runtime_error("Couldn't load shape!");

    if( psShape->nParts > 0 && psShape->panPartStart[0] != 0 )
      throw std::runtime_error("panPartStart[0] should be 0, but is not!");

    //Loop through all the vertices of the multipolygon
    int ringi = 0; //Which ring we're considering
    for(int j=0; j < psShape->nVertices; j++ ){
      if( ringi < psShape->nParts && psShape->panPartStart[ringi] == j ){
        if(!IsHole(psShape,ringi))
          mp.emplace_back();
        mp.back().emplace_back();
        ringi++;
      }
      
      //if(psShape->bMeasureIsUsed){
      mp.back().back().emplace_back(psShape->padfX[j], psShape->padfY[j]);
    }

    SHPDestroyObject( psShape );
  }

  SHPClose( hSHP );
}

void ReadShapeProj(GeoCollection &gc, std::string filename){
  if(filename.size()>=4 && filename.substr(filename.size()-4)==".shp")
    filename = filename.substr(0,filename.size()-4);
  filename += ".prj";

  std::ifstream fin(filename);
  std::stringstream buffer;
  buffer << fin.rdbuf();
  gc.prj_str = buffer.str();
}

GeoCollection ReadShapefile(std::string filename){
  GeoCollection mgons;
  
  ReadShapes(mgons,filename);
  ReadShapeAttributes(mgons,filename);
  ReadShapeProj(mgons,filename);

  return mgons;
}






static void WriteShapes(const GeoCollection &gc, const std::string filename){
  SHPHandle hSHP = SHPCreate(filename.c_str(), SHPT_POLYGON);
  if(hSHP==NULL)
    throw std::runtime_error("Failed to create shapefile '" + filename + "'!");


  for(unsigned int id=0;id<gc.size();id++){
    std::vector<double> x;
    std::vector<double> y;
    std::vector<int> rings;

    for(const auto &poly: gc.at(id)){
      rings.push_back(x.size());
      for(const auto &p: poly.at(0)){
        x.push_back(p.x);
        y.push_back(p.y);
      }
      for(unsigned int r=1;r<poly.size();r++){
        rings.push_back(x.size());
        for(auto p=poly.at(r).rbegin();p!=poly.at(r).rend();p++){
          x.push_back(p->x);
          y.push_back(p->y);
        }
      }
    }

    SHPObject *psObject = SHPCreateObject(
      SHPT_POLYGON,         //Type
      id,                   //Id
      rings.size(),         //Parts/polygons
      rings.data(),         //Starts of rings
      NULL,                 //Assumes SHPP_RING
      x.size(),             //Number of vertices being passed
      x.data(),
      y.data(),
      NULL,                 //No z coordinates
      NULL                  //No measure vertices
    );

    SHPWriteObject( hSHP, -1, psObject );
    SHPDestroyObject( psObject );
  }
    
  SHPClose( hSHP );
}


DBFFieldType getType(const std::string &s){
  if(s.find(".")!=std::string::npos){ //Can be a double
    try {
      std::stod(s);
      return FTDouble;
    } catch(...) {
      return FTString;
    }
  } else { //Is an integer
    try {
      std::stoi(s);
      return FTInteger;
    } catch(...) {
      return FTString;
    }    
  }
}



void WriteShapeAttributes(const GeoCollection &gc, const std::string filename){
  DBFHandle hDBF = DBFCreate(filename.c_str());
  if(hDBF==NULL)
    throw std::runtime_error("Failed to create shapefile database '" + filename + "'!");

  class PropType {
   public:
    DBFFieldType type;
    int width;
    int field_id;
  };

  std::map<std::string, PropType> proptypes;

  for(const auto &poly: gc)
  for(const auto &prop: poly.props){
    if(!proptypes.count(prop.first)){
      auto &ptv = proptypes[prop.first];
      ptv.type  = getType(prop.second);
      ptv.width = (int)prop.second.size();
    } else {
      auto &ptv = proptypes[prop.first];
      if(ptv.type!=getType(prop.second))
        throw std::runtime_error("Property types for shapefile output don't match!");
      ptv.width = std::max(ptv.width,(int)prop.second.size());
    }
  }

  for(auto &p: proptypes){
    p.second.field_id = DBFAddField(hDBF, p.first.c_str(), p.second.type, p.second.width, (p.second.type==FTDouble)?10:0);
    if(p.second.field_id==-1)
      throw std::runtime_error("Failed to add field '"+p.first+"' to shapefile dbf!");
  }

  for(unsigned int id=0;id<gc.size();id++){
    for(const auto &prop: gc.at(id).props){
      const auto &atv = proptypes[prop.first];
      switch(atv.type){
        case FTDouble:  DBFWriteDoubleAttribute (hDBF, id, atv.field_id, std::stod(prop.second)); break;
        case FTString:  DBFWriteStringAttribute (hDBF, id, atv.field_id, prop.second.c_str()); break;
        case FTInteger: DBFWriteIntegerAttribute(hDBF, id, atv.field_id, std::stoi(prop.second)); break;
        default:
          throw std::runtime_error("Unknown field type encountered in shapefile output!");
      }
    }
  }

  DBFClose( hDBF );
}

void WriteShapeScores(const GeoCollection &gc, const std::string filename){
  DBFHandle hDBF = DBFOpen(filename.c_str(), "rb+");
  if(hDBF==NULL)
    throw std::runtime_error("Failed to create shapefile database '" + filename + "'!");

  const int fields = DBFGetFieldCount(hDBF);

  std::set<std::string> scoreset; //Gather all scores
  for(const auto &poly: gc)
  for(const auto &score: poly.scores)
    scoreset.insert(score.first);

  std::vector< std::pair<std::string, int> > scorefields; //Field name + id

  for(const auto &s: scoreset){
    int ret = DBFAddField(hDBF, s.c_str(), FTDouble, 40, 10);
    if(ret==-1)
      throw std::runtime_error("Failed to add field '"+s+"' to shapefile dbf!");
    scorefields.emplace_back(s, ret);
  }

  for(unsigned int id=0;id<gc.size();id++){
    for(unsigned int s=0;s<scorefields.size();s++){
      const auto &sf = scorefields.at(s);
      if(gc.at(id).scores.count(sf.first))
        DBFWriteDoubleAttribute(hDBF, id, sf.second, gc.at(id).scores.at(sf.first));
      else
        DBFWriteNULLAttribute(hDBF, id, sf.second);
    }
  }

  DBFClose( hDBF );
}

void WriteShapeProj(const GeoCollection &gc, std::string filename){
  if(gc.prj_str.empty())
    return;

  if(filename.size()>=4 && filename.substr(filename.size()-4)==".shp")
    filename = filename.substr(0,filename.size()-4);
  filename += ".prj";

  std::ofstream fout(filename);
  fout<<gc.prj_str;
}


void WriteShapefile(const GeoCollection &gc, const std::string filename){
  WriteShapes(gc,filename);
  WriteShapeAttributes(gc,filename);
  WriteShapeScores(gc,filename);
  WriteShapeProj(gc,filename);
}

}
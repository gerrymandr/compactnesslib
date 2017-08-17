#include "shapefile.hpp"
#include "shapelib/shapefil.h"
#include <iostream>
#include <stdexcept>
#include <memory>
#include <map>
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
  unsigned int j = size-1;
  for(unsigned int i=0;i<size;i++){
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
      char szFormat[32];
      int  nWidth;
      int  nDecimals;

      //TODO: Inefficient? We could do this once for each field above and then
      //refer to a vector of the stored field info.
      const DBFFieldType eType = DBFGetFieldInfo(hDBF, i, szTitle, &nWidth, &nDecimals);
            
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
  int bValidate     = false;
  const char *pszPlus;
  double adfMinBound[4];
  double adfMaxBound[4];
  int nPrecision = 15;

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


GeoCollection ReadShapefile(std::string filename){
  GeoCollection mgons;
  
  ReadShapes(mgons,filename);
  ReadShapeAttributes(mgons,filename);

  return mgons;
}





// void WriteShapefile(const GeoCollection &mps, std::string filename, std::string layername){
//   const char *pszDriverName = "ESRI Shapefile";
//   //const char *pszDriverName = "GeoJSON";

//   GDALAllRegister();

//   GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
//   if( poDriver == NULL )
//     throw std::runtime_error("Driver "+std::string(pszDriverName)+" not available!");

//   GDALDataset *poDS = poDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL );
//   if( poDS == NULL )
//     throw std::runtime_error("Creation of file '"+filename+"' failed!");

//   OGRLayer *poLayer;
//   poLayer = poDS->CreateLayer(layername.c_str(), NULL, wkbMultiPolygon, NULL );
//   if( poLayer == NULL )
//     throw std::runtime_error("Layer creation failed!");

//   for(const auto &pr: mps.front().props){
//     std::unique_ptr<OGRFieldDefn> oField;

//     std::string trunc_name = pr.first.substr(0,10);

//     if(typeid(int) == pr.second.type()) {
//       oField.reset( new OGRFieldDefn( trunc_name.c_str(), OFTInteger) );
//     } else if(typeid(long) == pr.second.type() || typeid(long long) == pr.second.type()) {
//       oField.reset( new OGRFieldDefn( trunc_name.c_str(), OFTInteger64) );
//     } else if(typeid(double) == pr.second.type()) {
//       oField.reset( new OGRFieldDefn( trunc_name.c_str(), OFTReal) );
//     } else if(typeid(std::string) == pr.second.type()) {
//       oField.reset( new OGRFieldDefn( trunc_name.c_str(), OFTString) );
//     } else {
//       throw std::runtime_error("Unrecognized field type (field creation) '" + std::string(pr.second.type().name()) + "'!");
//     }

//     oField->SetWidth(32);
//     oField->SetPrecision(16);

//     if( poLayer->CreateField( oField.get(), false ) != OGRERR_NONE )
//       throw std::runtime_error("Failed to create a field!");
//   }

//   for(const auto &mp: mps){
//     OGRFeature *poFeature;
//     poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );

//     if(poFeature == NULL)
//       throw std::runtime_error("Failed to create feature!");


//     //////////////////////
//     //This method is so stupid, but the only way I could find to get it to work
//     for(int i=0;i<poFeature->GetFieldCount();i++){
//       std::string search = poFeature->GetFieldDefnRef(i)->GetNameRef();
//       search             = search.substr(0,10);
      
//       auto found = mp.props.find(search);

//       if(found==mp.props.end()){
//         std::cerr<<"Cound not find key '"<<std::string(search)<<"'"<<std::endl;
//         continue;
//       }

//       if(typeid(int) == found->second.type()) {
//         poFeature->SetField(i, std::any_cast<int>(found->second) );
//       } else if(typeid(long) == found->second.type() || typeid(long long) == found->second.type()) {
//         poFeature->SetField(i, std::any_cast<GIntBig>(found->second) );
//       } else if(typeid(double) == found->second.type()) {
//         poFeature->SetField(i, std::any_cast<double>(found->second) );
//       } else if(typeid(std::string) == found->second.type()) {
//         poFeature->SetField(i, std::any_cast<std::string>(found->second).c_str() );
//       } else {
//         throw std::runtime_error("Unrecognized field type (value write) '" + std::string(found->second.type().name()) + "'!");
//       }      
//     }

//     ///////////////////////////
//     //This should work, but does not
//     // for(const auto &pr: mp.props){
//     //   if(poFeature->GetFieldIndex(pr.first.c_str())){
//     //     std::cerr<<"Warning! Could not find field '"<<pr.first<<"'"<<std::endl;
//     //     continue;
//     //   }

//     //   std::string trunc_name = pr.first.substr(0,10);

//     //   if(typeid(int) == pr.second.type()) {
//     //     poFeature->SetField(trunc_name.c_str(), std::any_cast<int>(pr.second) );
//     //   } else if(typeid(long) == pr.second.type()) {
//     //     poFeature->SetField(trunc_name.c_str(), std::any_cast<GIntBig>(pr.second) );
//     //   } else if(typeid(double) == pr.second.type()) {
//     //     poFeature->SetField(trunc_name.c_str(), std::any_cast<double>(pr.second) );
//     //   } else if(typeid(std::string) == pr.second.type()) {
//     //     poFeature->SetField(trunc_name.c_str(), std::any_cast<std::string>(pr.second).c_str() );
//     //   } else {
//     //     std::cerr<<"Unrecognized property type '"<<pr.second.type().name()<<"'!"<<std::endl;
//     //   }
//     // }

    
//     auto out_mp = MultiPolygonToOGR(mp);
//     poFeature->SetGeometryDirectly(out_mp);

//     if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
//       throw std::runtime_error("Failed to create feature in output!");

//     OGRFeature::DestroyFeature( poFeature );
//   }

//   GDALClose( poDS );
// }

}
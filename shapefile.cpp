#include "shapefile.hpp"
#include "shapelib/shapefil.h"
#include <iostream>
#include <stdexcept>
#include <memory>
#include <any>
#include <map>
#include "geom.hpp"

OGRPoint PointToOGR(const Point2D &p){
  OGRPoint pt;
  pt.setX(p.x);
  pt.setY(p.y);
  return pt;
}

OGRLinearRing* RingToOGR(const Ring &r){
  OGRLinearRing* olr = new OGRLinearRing();
  for(const auto &p: r)
    olr->addPoint(p.x,p.y);
  return olr;
}

OGRPolygon* PolygonToOGR(const Polygon &p){
  OGRPolygon* poly = new OGRPolygon();
  poly->addRingDirectly(RingToOGR(p.outer));
  for(const auto &h: p.holes)
    poly->addRingDirectly(RingToOGR(h));
  return poly;
}

OGRMultiPolygon* MultiPolygonToOGR(const MultiPolygon &mp){
  OGRMultiPolygon* omp = new OGRMultiPolygon();
  for(const auto &p: mp)
    omp->addGeometryDirectly(PolygonToOGR(p));
  return omp;
}

std::map<std::string,std::any> ReadOGRFields(OGRFeature &poFeature, OGRFeatureDefn &ofd){
  std::map<std::string,std::any> props;
  for(int iField = 0; iField < ofd.GetFieldCount(); iField++ ){
    
    OGRFieldDefn *poFieldDefn = ofd.GetFieldDefn( iField );

    auto &prop = props[poFieldDefn->GetNameRef()];

    if( poFieldDefn->GetType() == OFTInteger )
      prop = poFeature.GetFieldAsInteger( iField );
    else if( poFieldDefn->GetType() == OFTInteger64 )
      prop = poFeature.GetFieldAsInteger64( iField );
    else if( poFieldDefn->GetType() == OFTReal )
      prop = poFeature.GetFieldAsDouble(iField);
    else if( poFieldDefn->GetType() == OFTString )
      prop = std::string(poFeature.GetFieldAsString(iField));
    else
      prop = std::string(poFeature.GetFieldAsString(iField));
  }

  return props;
}

Ring ReadORGRing(const OGRLinearRing &oring){
  Ring ring;
  for(int i=0;i<oring.getNumPoints();i++)
    ring.emplace_back(oring.getX(i),oring.getY(i));
  return ring;
}

Polygon ReadOGRPolygon(const OGRPolygon &opoly){
  Polygon poly;
  poly.outer = ReadORGRing(*opoly.getExteriorRing());

  for(int i=0;i<opoly.getNumInteriorRings();i++)
    poly.holes.push_back(ReadORGRing(*opoly.getInteriorRing(i)));

  return poly;
}

MultiPolygon ReadOGRMultiPolygon(const OGRMultiPolygon &ompoly){
  MultiPolygon mpoly;

  for(int i=0;i<ompoly.getNumGeometries();i++){
    const auto g = ompoly.getGeometryRef(i);
    if(wkbFlatten(g->getGeometryType())!=wkbPolygon)
      throw std::runtime_error("MultiPolygon had a non-polygon subunit!");
    mpoly.push_back(ReadOGRPolygon(*static_cast<const OGRPolygon*>(g)));
  }

  return mpoly;
}

GeoCollection ReadShapefile(std::string filename, std::string layername){
  int nShapeType;
  int nEntities;
  int i;
  int iPart;
  int bValidate     = 0;
  int nInvalidCount = 0;
  int bHeaderOnly   = 0;
  const char *pszPlus;
  double adfMinBound[4];
  double adfMaxBound[4];
  int nPrecision = 15;


  SHPHandle *poDS = ShpOpen(filename.c_str(), "rb");
  if( poDS == NULL )
    throw std::runtime_error("Failed to open shapefile '" + filename + "'!");

  GeoCollection mgons;

  OGRFeature *poFeature;
  poLayer->ResetReading();
  while(int i=0;i<nEntities;i++){
    SHPObject *const psShape = SHPReadObject( hSHP, i );

    if(psShape==NULL)
      throw std::runtime_error("Couldn't load shape!");







      if( psShape->nParts > 0 && psShape->panPartStart[0] != 0 )
      {
          fprintf( stderr, "panPartStart[0] = %d, not zero as expected.\n",
                   psShape->panPartStart[0] );
      }

      for( j = 0, iPart = 1; j < psShape->nVertices; j++ )
      {
          const char  *pszPartType = "";

          if( j == 0 && psShape->nParts > 0 )
              pszPartType = SHPPartTypeName( psShape->panPartType[0] );
          
          if( iPart < psShape->nParts
              && psShape->panPartStart[iPart] == j )
          {
              pszPartType = SHPPartTypeName( psShape->panPartType[iPart] );
              iPart++;
              pszPlus = "+";
          }
          else
              pszPlus = " ";

          if( psShape->bMeasureIsUsed )
              printf("   %s (%.*g,%.*g, %.*g, %.*g) %s \n",
                     pszPlus,
                     nPrecision, psShape->padfX[j],
                     nPrecision, psShape->padfY[j],
                     nPrecision, psShape->padfZ[j],
                     nPrecision, psShape->padfM[j],
                     pszPartType );
          else
              printf("   %s (%.*g,%.*g, %.*g) %s \n",
                     pszPlus,
                     nPrecision, psShape->padfX[j],
                     nPrecision, psShape->padfY[j],
                     nPrecision, psShape->padfZ[j],
                     pszPartType );
      }

      if( bValidate )
      {
          int nAltered = SHPRewindObject( hSHP, psShape );

          if( nAltered > 0 )
          {
              printf( "  %d rings wound in the wrong direction.\n",
                      nAltered );
              nInvalidCount++;
          }
      }
      
      SHPDestroyObject( psShape );
    }

    SHPClose( hSHP );

    if( bValidate )
    {
        printf( "%d object has invalid ring orderings.\n", nInvalidCount );
    }








    auto geom_type = wkbFlatten(poGeometry->getGeometryType());

    if( geom_type == wkbPolygon ){
      //Create a new MultiPolygon
      mgons.emplace_back();
      //Load that with the polygon
      mgons.back().push_back(ReadOGRPolygon(*static_cast<OGRPolygon*>(poGeometry)));
    } else if( geom_type == wkbMultiPolygon){
      mgons.push_back(ReadOGRMultiPolygon(*static_cast<OGRMultiPolygon*>(poGeometry)));
    } else {
      throw std::runtime_error("Unrecognised geometry of type: "+std::to_string(wkbFlatten(poGeometry->getGeometryType())));
    }

    mgons.back().props = props;

    OGRFeature::DestroyFeature( poFeature );
  }
  SHPClose( poDS );

  return mgons;
}






void WriteShapefile(const GeoCollection &mps, std::string filename, std::string layername){
  const char *pszDriverName = "ESRI Shapefile";
  //const char *pszDriverName = "GeoJSON";

  GDALAllRegister();

  GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
  if( poDriver == NULL )
    throw std::runtime_error("Driver "+std::string(pszDriverName)+" not available!");

  GDALDataset *poDS = poDriver->Create(filename.c_str(), 0, 0, 0, GDT_Unknown, NULL );
  if( poDS == NULL )
    throw std::runtime_error("Creation of file '"+filename+"' failed!");

  OGRLayer *poLayer;
  poLayer = poDS->CreateLayer(layername.c_str(), NULL, wkbMultiPolygon, NULL );
  if( poLayer == NULL )
    throw std::runtime_error("Layer creation failed!");

  for(const auto &pr: mps.front().props){
    std::unique_ptr<OGRFieldDefn> oField;

    std::string trunc_name = pr.first.substr(0,10);

    if(typeid(int) == pr.second.type()) {
      oField.reset( new OGRFieldDefn( trunc_name.c_str(), OFTInteger) );
    } else if(typeid(long) == pr.second.type() || typeid(long long) == pr.second.type()) {
      oField.reset( new OGRFieldDefn( trunc_name.c_str(), OFTInteger64) );
    } else if(typeid(double) == pr.second.type()) {
      oField.reset( new OGRFieldDefn( trunc_name.c_str(), OFTReal) );
    } else if(typeid(std::string) == pr.second.type()) {
      oField.reset( new OGRFieldDefn( trunc_name.c_str(), OFTString) );
    } else {
      throw std::runtime_error("Unrecognized field type (field creation) '" + std::string(pr.second.type().name()) + "'!");
    }

    oField->SetWidth(32);
    oField->SetPrecision(16);

    if( poLayer->CreateField( oField.get(), false ) != OGRERR_NONE )
      throw std::runtime_error("Failed to create a field!");
  }

  for(const auto &mp: mps){
    OGRFeature *poFeature;
    poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );

    if(poFeature == NULL)
      throw std::runtime_error("Failed to create feature!");


    //////////////////////
    //This method is so stupid, but the only way I could find to get it to work
    for(int i=0;i<poFeature->GetFieldCount();i++){
      std::string search = poFeature->GetFieldDefnRef(i)->GetNameRef();
      search             = search.substr(0,10);
      
      auto found = mp.props.find(search);

      if(found==mp.props.end()){
        std::cerr<<"Cound not find key '"<<std::string(search)<<"'"<<std::endl;
        continue;
      }

      if(typeid(int) == found->second.type()) {
        poFeature->SetField(i, std::any_cast<int>(found->second) );
      } else if(typeid(long) == found->second.type() || typeid(long long) == found->second.type()) {
        poFeature->SetField(i, std::any_cast<GIntBig>(found->second) );
      } else if(typeid(double) == found->second.type()) {
        poFeature->SetField(i, std::any_cast<double>(found->second) );
      } else if(typeid(std::string) == found->second.type()) {
        poFeature->SetField(i, std::any_cast<std::string>(found->second).c_str() );
      } else {
        throw std::runtime_error("Unrecognized field type (value write) '" + std::string(found->second.type().name()) + "'!");
      }      
    }

    ///////////////////////////
    //This should work, but does not
    // for(const auto &pr: mp.props){
    //   if(poFeature->GetFieldIndex(pr.first.c_str())){
    //     std::cerr<<"Warning! Could not find field '"<<pr.first<<"'"<<std::endl;
    //     continue;
    //   }

    //   std::string trunc_name = pr.first.substr(0,10);

    //   if(typeid(int) == pr.second.type()) {
    //     poFeature->SetField(trunc_name.c_str(), std::any_cast<int>(pr.second) );
    //   } else if(typeid(long) == pr.second.type()) {
    //     poFeature->SetField(trunc_name.c_str(), std::any_cast<GIntBig>(pr.second) );
    //   } else if(typeid(double) == pr.second.type()) {
    //     poFeature->SetField(trunc_name.c_str(), std::any_cast<double>(pr.second) );
    //   } else if(typeid(std::string) == pr.second.type()) {
    //     poFeature->SetField(trunc_name.c_str(), std::any_cast<std::string>(pr.second).c_str() );
    //   } else {
    //     std::cerr<<"Unrecognized property type '"<<pr.second.type().name()<<"'!"<<std::endl;
    //   }
    // }

    
    auto out_mp = MultiPolygonToOGR(mp);
    poFeature->SetGeometryDirectly(out_mp);

    if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
      throw std::runtime_error("Failed to create feature in output!");

    OGRFeature::DestroyFeature( poFeature );
  }

  GDALClose( poDS );
}

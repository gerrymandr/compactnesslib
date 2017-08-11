#include "shapefile.hpp"
#include <ogrsf_frmts.h>
#include <iostream>
#include <stdexcept>
#include <memory>
#include "geom.hpp"

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

MultiPolygons ReadShapefile(std::string filename, std::string layername){
  GDALAllRegister();
  GDALDataset *poDS;
  poDS = (GDALDataset*) GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
  if( poDS == NULL )
    throw std::runtime_error("Failed to open shapefile '" + filename + "'!");

  MultiPolygons mgons;

  OGRLayer *poLayer;
  poLayer = poDS->GetLayerByName(layername.c_str());
  if(poLayer==NULL)
    throw std::runtime_error("Specified layer name not found!");

  OGRFeature *poFeature;
  poLayer->ResetReading();
  while( (poFeature = poLayer->GetNextFeature()) != NULL ){
    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
    int iField;
    for( iField = 0; iField < poFDefn->GetFieldCount(); iField++ ){
      OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( iField );
      if( poFieldDefn->GetType() == OFTInteger )
          printf( "%d,", poFeature->GetFieldAsInteger( iField ) );
      else if( poFieldDefn->GetType() == OFTInteger64 )
          printf( CPL_FRMT_GIB ",", poFeature->GetFieldAsInteger64( iField ) );
      else if( poFieldDefn->GetType() == OFTReal )
          printf( "%.3f,", poFeature->GetFieldAsDouble(iField) );
      else if( poFieldDefn->GetType() == OFTString )
          printf( "%s,", poFeature->GetFieldAsString(iField) );
      else
          printf( "%s,", poFeature->GetFieldAsString(iField) );
    }
    printf("\n");
    OGRGeometry *poGeometry;
    poGeometry = poFeature->GetGeometryRef();

    if(poGeometry==NULL)
      continue;

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
    OGRFeature::DestroyFeature( poFeature );
  }
  GDALClose( poDS );

  return mgons;
}

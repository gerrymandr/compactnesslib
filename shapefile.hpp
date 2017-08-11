#ifndef _shapefill_hpp_
#define _shapefill_hpp_

#include <string>
#include "geom.hpp"

GeoCollection ReadShapefile(std::string filename, std::string layername);
void WriteShapefile(const GeoCollection &mps, std::string filename, std::string layername);

#endif
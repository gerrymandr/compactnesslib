#ifndef _shapefill_hpp_
#define _shapefill_hpp_

#include <string>
#include "geom.hpp"

MultiPolygons ReadShapefile(std::string filename, std::string layername);
void WriteShapefile(const MultiPolygons &mps, std::string filename, std::string layername);

#endif
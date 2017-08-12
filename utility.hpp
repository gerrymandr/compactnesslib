#ifndef _utility_hpp_
#define _utility_hpp_

#include <ostream>
#include <any>

namespace complib {

std::ostream& operator<<(std::ostream &out, const std::any &propval);

}

#endif
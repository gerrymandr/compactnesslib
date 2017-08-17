#ifndef _props_hpp_
#define _props_hpp_

//#include <any>
#include <map>
#include <string>

namespace complib {
  //typedef std::map<std::string,std::any> Props;
  typedef std::map<std::string, std::string> Props;

  void PrintProps(const Props &ps);

  // std::ostream& operator<<(std::ostream &out, const std::any &propval);
}


#endif
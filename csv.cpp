#include "csv.hpp"
#include <string>
#include <set>
#include <sstream>

namespace complib {

std::string OutScoreCSV(const GeoCollection &gc, std::string id) {
  std::ostringstream oss;

  const bool use_id = !id.empty();

  std::set<std::string> score_names;
  for(const auto &mp: gc)
  for(const auto &s: mp.scores)
    score_names.insert(s.first);

  oss<<"id";
  for(const auto &sn: score_names)
    oss<<","<<sn;
  oss<<"\n";

  for(unsigned int i=0;i<gc.size();i++){
    if(use_id)
      oss<<gc[i].props.at(id);
    else
      oss<<i;

    for(const auto &sn: score_names){
      oss<<",";
      if(gc[i].scores.count(sn)){
        oss<<gc[i].scores.at(sn);
      } else {
        oss<<-9999;
      }
    }
    oss<<"\n";
  }

  return oss.str();
}

}
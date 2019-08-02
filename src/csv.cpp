#include <string>
#include <set>
#include <sstream>
#include <stdexcept>
#include <iomanip>

#include <compactnesslib/csv.hpp>

namespace complib {

std::string OutScoreCSV(const GeoCollection &gc, std::string id) {
  std::ostringstream oss;

  const bool use_id = !id.empty();

  std::set<std::string> scores_used;
  for(const auto &mp: gc)
  for(const auto &s: mp.scores)
    scores_used.insert(s.first);

  oss<<"id";
  for(const auto &sn: scores_used)
    oss<<","<<sn;
  oss<<"\n";

  for(unsigned int i=0;i<gc.size();i++){
    if(use_id){
      if(gc[i].props.count(id))
        oss<<gc[i].props.at(id);
      else
        throw std::runtime_error("Failed to find id property '"+id+"'");
    } else {
      oss<<i;
    }

    for(const auto &sn: scores_used){
      oss<<",";
      if(gc[i].scores.count(sn)){
        oss<<std::fixed<<std::setprecision(5)<<gc[i].scores.at(sn);
      } else {
        oss<<-9999;
      }
    }
    oss<<"\n";
  }

  return oss.str();
}

}
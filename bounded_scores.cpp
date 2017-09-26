#include "bounded_scores.hpp"
#include "geom.hpp"
#include <cmath>
#include <vector>
#include <stdexcept>
#include <unordered_map>
#include "lib/clipper.hpp"

#include <sstream>  //TODO
#include <iostream> //TODO

namespace complib {

namespace cl = ClipperLib;

double ScoreConvexHullPTB(const MultiPolygon &mp, const MultiPolygon &border){
  const double area      = areaIncludingHoles(mp);
  const double hull_area = IntersectionArea(mp.getHull(),border);
  double ratio = area/hull_area;
  if(ratio>1)
    ratio = 1;
  return ratio;
}



double ScoreReockPTB(const MultiPolygon &mp, const MultiPolygon &border){
  const auto   circle = GetBoundingCircle(mp);
  const auto   iarea  = IntersectionArea(circle, border);
  const double area   = areaIncludingHoles(mp);

  double ratio = area/iarea;
  if(ratio>1)
    ratio = 1;
  return ratio;
}



// std::string OutputPath(const cl::Path &path){
//   std::ostringstream out;
//   out << "POLYGON((";
//   for(unsigned int pi=0;pi<path.size();pi++){
//     out<<path[pi].X<<" "<<path[pi].Y;
//     if(pi<path.size()-1)
//       out<<",";
//   }
//   out<<"))";
//   return out.str();
// }

// std::string OutputPaths(const cl::Paths &paths){
//   std::ostringstream out;
//   out << "MULTIPOLYGON(";
//   for(unsigned int psi=0;psi<paths.size();psi++){
//     out<<"((";
//     for(unsigned int pi=0;pi<paths[psi].size();pi++){
//       out<<paths[psi][pi].X<<" "<<paths[psi][pi].Y;
//       if(pi<paths[psi].size()-1)
//         out<<",";
//     }
//     out<<"))";
//     if(psi<paths.size()-1)
//       out<<",";
//   }
//   out<<")";
//   return out.str();
// }


cl::Paths GetClipperRing(const cl::Paths &unit, const int pad_amount){
  //Grow the unit
  cl::Paths grown;
  {
    cl::ClipperOffset co;
    co.AddPaths(unit, cl::jtRound, cl::etClosedPolygon);
    co.Execute(grown, pad_amount);
  }

  //Shrink the unit
  cl::Paths shrunk;
  {
    cl::ClipperOffset co;
    co.AddPaths(unit, cl::jtRound, cl::etClosedPolygon);
    co.Execute(shrunk, -pad_amount);
  }

  //Get a unit ring
  cl::Paths ring;
  {
    cl::Clipper clpr;
    clpr.AddPaths(grown, cl::ptSubject, true);
    clpr.AddPaths(shrunk, cl::ptClip, true);
    cl::Paths solution;
    clpr.Execute(cl::ctDifference, ring, cl::pftEvenOdd, cl::pftEvenOdd);
  }

  return ring;
}

double ScoreBorderAreaUncertainty(const MultiPolygon &mp, const MultiPolygon &border){
  static int ringnum=0;
  ringnum++;
  //Amount by which we will grow the subunit
  const int pad_amount = 1000; //metres

  //If we've already determined the subunit is not an exterior child, then no
  //further calculation is necessary.
  if(mp.props.count("EXTCHILD") && mp.props.at("EXTCHILD")=="F")
    return 0;

  //Get paths for both the subunit and its superunit
  const auto paths_mp = ConvertToClipper(mp,false);
  const auto paths_bo = ConvertToClipper(border,false);

  //std::cerr<<ringnum<<",pathmp,\""<<OutputPaths(paths_mp)<<"\""<<std::endl;
  //std::cerr<<ringnum<<",pathbo,\""<<OutputPaths(paths_bo)<<"\""<<std::endl;

  //if(paths_mp.at(0).front()!=paths_mp.at(0).back())
    //std::cerr<<"NOT A RING!"<<std::endl;

  const auto mp_ring = GetClipperRing(paths_mp, pad_amount);
  const auto bo_ring = GetClipperRing(paths_bo, pad_amount);

  //std::cerr<<ringnum<<",mp_ring,\""<<OutputPaths(mp_ring)<<"\""<<std::endl;
  //std::cerr<<ringnum<<",bo_ring,\""<<OutputPaths(bo_ring)<<"\""<<std::endl;

  //Get the intersection of the border rings - uncertainty can only occur here
  cl::Paths ring_intersection;
  {
    cl::Clipper clpr;
    clpr.AddPaths(bo_ring, cl::ptSubject, true);
    clpr.AddPaths(mp_ring, cl::ptClip, true);
    cl::Paths solution;
    clpr.Execute(cl::ctIntersection, ring_intersection, cl::pftEvenOdd, cl::pftEvenOdd);
  }

  //std::cerr<<ringnum<<",ringint,\""<<OutputPaths(ring_intersection)<<"\""<<std::endl;


  //XOR the subunit and the superunit. This gives us the border uncertainty, but
  //also the entire rest of the superunit!
  cl::Paths xored;
  {
    cl::Clipper clpr;
    clpr.AddPaths(paths_mp, cl::ptSubject, true);
    clpr.AddPaths(paths_bo, cl::ptClip, true);
    cl::Paths solution;
    clpr.Execute(cl::ctXor, xored, cl::pftEvenOdd, cl::pftEvenOdd);
  }

  //std::cerr<<ringnum<<",xored,\""<<OutputPaths(xored)<<"\""<<std::endl;

  //Get the intersection of the border ring and the xored area - this is an upper
  //bound on the uncertain area
  cl::Paths isect;
  {
    cl::Clipper clpr;
    clpr.AddPaths(ring_intersection, cl::ptSubject, true);
    clpr.AddPaths(xored, cl::ptClip, true);
    cl::Paths solution;
    clpr.Execute(cl::ctIntersection, isect, cl::pftEvenOdd, cl::pftEvenOdd);
  }

  //std::cerr<<ringnum<<",isect,\""<<OutputPaths(isect)<<"\""<<std::endl;

  double area = 0;
  for(const auto &path: isect)
    area += cl::Area(path);

  return area;
}



void CalculateAllBoundedScores(
  GeoCollection &subunits,
  const GeoCollection &superunits,
  const std::string join_on
){
  CalculateListOfBoundedScores(subunits, superunits, join_on, getListOfBoundedScores());
}



void CalculateListOfBoundedScores(
  GeoCollection &subunits,
  const GeoCollection &superunits,
  const std::string join_on,
  std::vector<std::string> score_list
){
  if(score_list.empty())
    score_list = getListOfBoundedScores();
  else if(score_list.size()==1 && score_list.at(0)=="all")
    score_list = getListOfBoundedScores();


  if(join_on.empty() || superunits.size()==1){

    for(auto& sub: subunits){
      for(const auto &sn: score_list){
        if(bounded_score_map.count(sn))
          sub.scores[sn] = bounded_score_map.at(sn)(sub,superunits.at(0));
      }
    }

  } else {

    for(const auto &mp: subunits)
      if(!mp.props.count(join_on))
        throw std::runtime_error("At least one subunit was missing the joining attribute!");

    //A quick was to access superunits based on their key
    std::unordered_map<std::string, const MultiPolygon *> su_key;
    for(const auto &mp: superunits){
      if(!mp.props.count(join_on))
        throw std::runtime_error("At least one superunit was missing the joining attribute!");
      if(su_key.count(mp.props.at(join_on)))
        throw std::runtime_error("More than one superunit had the same key!");
      su_key[mp.props.at(join_on)] = &mp;
    }

    for(auto& sub: subunits){
      for(const auto &sn: score_list){
        if(bounded_score_map.count(sn))
          sub.scores[sn] = bounded_score_map.at(sn)(sub,*su_key.at(sub.props.at(join_on)));
      }
    }

  }
}



const std::vector<std::string>& getListOfBoundedScores(){
  static std::vector<std::string> score_names;
  if(!score_names.empty())
    return score_names;
  for(const auto &kv: bounded_score_map)
    score_names.push_back(kv.first);
  return score_names;
}



const bounded_score_map_t bounded_score_map({
  {"CvxHullPTB", ScoreConvexHullPTB},
  {"ReockPTB",   ScoreReockPTB},
  {"AreaUncert", ScoreBorderAreaUncertainty}
});

}

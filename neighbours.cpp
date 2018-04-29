#include "neighbours.hpp"
#include "geom.hpp"
#include "SpIndex.hpp"
#include "compactnesslib.hpp"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include "lib/doctest.h"

//TODO
#include <iostream>
#include <iomanip>

namespace complib {

const auto qnan = std::numeric_limits<double>::quiet_NaN();


typedef std::vector<double>               pointdimvec_t;
typedef std::vector<pointdimvec_t>        pointvec_t;
typedef std::vector<unsigned int>         ownervec_t;
typedef std::pair<ownervec_t, pointvec_t> owner_point_vec_t;


class SegmentGrid {
 private:
  std::vector<Segment> segments;
  typedef std::vector<int> segrefvec;

  //The grid. We assume that the data is sparse, since we are dealing with
  //polygons, and, therefore, use an unordered_map
  std::unordered_map<int, segrefvec> cells;

  double epsilon2;               //Threshold for two segments being neighbours

  double xmin;                   //Minimum x-coordinate of grid
  double ymin;                   //Minimum y-coordinate of grid
  double xmax;                   //Maximum x-coordinate of grid
  double ymax;                   //Maximum x-coordinate of grid
  double cellsize;               //Width & Height of a grid cell
  int    width;                  //Width of grid measured in cells
  int    height;                 //Height of grid measured in cells
  mutable int segcomp_count = 0; //Number of segment-segment comparisons

  bool doesThisVecHaveSegmentWithinEpsilon(
    const segrefvec &segvec,
    const Segment   &seg
  ) const {
    if(segvec.empty())
      return false;

    for(const auto &a: segvec){
      segcomp_count++;
      const auto dist = SegmentSegmentDistanceSquared(
        segments.at(a).first,
        segments.at(a).second,
        seg.first,
        seg.second
      );
      if(dist<=epsilon2)
        return true;
    }

    return false;
  }

 public:
  SegmentGrid() = default;

  SegmentGrid(
    const double xmin0,
    const double ymin0,
    const double xmax0,
    const double ymax0,
    const double cellsize0,
    const double epsilon0
  ){
    xmin     = xmin0;
    ymin     = ymin0;
    xmax     = xmax0;
    ymax     = ymax0;
    cellsize = cellsize0;
    epsilon2 = epsilon0*epsilon0;

    //Round up! Otherwise a segment end-point may fall off the grid.
    width    = std::ceil((xmax-xmin)/cellsize);
    height   = std::ceil((ymax-ymin)/cellsize);

    //Set hash table size to initially 10% of full size, based on some
    //experiments
    cells.reserve(width*height/10);
  }

  void printSegments() const {
    for(const auto &seg: segments)
      std::cerr<<"\t"
      <<"("<<seg.first.x<<","<<seg.first.y<<")"
      <<"-"
      <<"("<<seg.second.x<<","<<seg.second.y<<")"
      <<std::endl;
  }

  void addSegment(
    const Point2D &a,  //End point of segment
    const Point2D &b   //End point of segment
  ){
    //Get coordinates of cells
    const int ax = (a.x-xmin)/cellsize;
    const int ay = (a.y-ymin)/cellsize;
    const int bx = (b.x-xmin)/cellsize;
    const int by = (b.y-ymin)/cellsize;

    if(std::abs(ax-bx)>=2 || std::abs(ay-by)>=2)
      throw std::runtime_error("Segment endpoints are not in adjacent cells!");

    //Get cell addresses
    const int acelli = ay*width+ax;
    const int bcelli = by*width+bx;

    //Insert cells if not already present, get refernces to the cells
    segrefvec& acell = cells.emplace(acelli, segrefvec()).first->second; //.first is the map object (<int,Cell>), .second is the cell
    segrefvec& bcell = cells.emplace(bcelli, segrefvec()).first->second; //.first is the map object (<int,Cell>), .second is the cell

    segments.emplace_back(a,b);

    acell.emplace_back( segments.size()-1 );  
    bcell.emplace_back( segments.size()-1 );  
  }

  void addSegment(const Segment &seg){
    addSegment(seg.first, seg.second);
  }

  bool isThereANearbySegment(const Segment &seg){
    const int dx[9] = {0, -1, -1,  0,  1, 1, 1, 0, -1}; // x offsets of D8 neighbours, from a central cell
    const int dy[9] = {0,  0, -1, -1, -1, 0, 1, 1,  1}; // y offsets of D8 neighbours, from a central cell

    //Get coordinates of cells
    const int ax = (seg.first.x-xmin)/cellsize;
    const int ay = (seg.first.y-ymin)/cellsize;
    const int bx = (seg.second.x-xmin)/cellsize;
    const int by = (seg.second.y-ymin)/cellsize;

    if( ! (std::abs(ax-bx)<=1 && std::abs(ay-by)<=1) )
      throw std::runtime_error("Segment endpoints are not in adjacent cells!");

    //Get cell addresses
    const int acelli = ay*width+ax;
    const int bcelli = by*width+bx;

    //Since acell and bcell are neighbours, they share neighbouring cells. Here,
    //we keep track of the cells we've considered as neighbours to avoid
    //repeating searches. Since this is a small dataset, we use a vector to
    //avoid the overhead of using a set. Since each cell is surrounded by 8
    //neighbours plus itself, there are at most 18 entries in the list.
    // std::vector<int> explored;
    // explored.reserve(18);

    auto nloop = [&](const int cx, const int cy){
      //Loop through neighbours
      for(int n=0;n<=8;n++){
        const auto nx = cx+dx[n];                   //x-coordinate of neighbour
        const auto ny = cy+dy[n];                   //y-coordinate of neighbour

        if(nx<0 || ny<0 || nx==width || ny==height) //Is cell outside grid bounds?
          continue;                                 //Yep

        const int nidx = ny*width+nx;               //Index of neighbour cell

        if(!cells.count(nidx))                      //Does this cell exist?
          continue;                                 //Nope

        //Has this neighbouring cell already been explored?
        // if(std::find(explored.begin(),explored.end(),nidx)!=explored.end())
        //   continue;
        // explored.push_back(nidx);

        const auto& ncell = cells.at(nidx);

        if(doesThisVecHaveSegmentWithinEpsilon(ncell, seg))
          return true;
      }
      return false;
    };

    //Check if one subunit segment end point is near a parent segment
    if(nloop(ax,ay))
      return true;
    //Check if the other subunit segment end point is near a parent segment
    if(acelli!=bcelli && nloop(bx,by))
      return true;

    return false;    
  }

  int sizeCells() const {
    return cells.size();
  }

  int getWidth() const {
    return width;
  }

  int getHeight() const {
    return height;
  }

  int numSegSegComparisons() const {
    return segcomp_count;
  }

  void printCellCounts() const { //TODO: Cut
    for(int y=0;y<height;y++){
      for(int x=0;x<width;x++){
        const int c = y*width+x;
        if(cells.count(c)){
          const auto &cell = cells.at(c);
          std::cout<<std::setw(5)<<(cell.size())<<" ";
        } else {
          std::cout<<std::setw(5)<<"0"<<" ";
        }
      }
      std::cout<<"\n";
    }
    std::cout<<"\n";
  }
};


owner_point_vec_t GetDensifiedBorders(const GeoCollection &gc, const double maxdist){
  owner_point_vec_t temp;

  //Load all the unit's points into a vector that will be loaded into a kd- tree
  //for quick nearest-neighbour lookups
  for(unsigned int mpi=0;mpi<gc.size();mpi++){
    //Loop through all of the points of the superunit
    for(const auto &poly: gc.at(mpi))
    for(const auto &ring: poly)
    for(unsigned int i=0;i<ring.size();i++){
      const auto &a   = ring.at(i);
      const auto &b   = ring.at((i+1)%ring.size()); //Loop around to beginning
      const auto dist = EuclideanDistance(a,b);

      //Portion of distance between the two points taken up by maxdist - used as
      //step size for interpolation
      const auto step = maxdist/dist;

      //Calculate intermediate points using linear interpolation via method of
      //weighted averages
      int    si = 0; //Which portion of the weighted average we are on - prevents build up of floating errors
      double st = 0; //Portion of the average coming from start vs end point
      do {
        temp.first.emplace_back(mpi);
        temp.second.emplace_back(pointdimvec_t{
          (1-st)*a.x + st*b.x,
          (1-st)*a.y + st*b.y
        });
        si++;
        st = si*step;
      } while (st<1);

      temp.first.emplace_back(mpi);
    }
  }

  return temp;
}



std::vector<Segment> GetDensifiedBorderSegments(const MultiPolygon &mp, const double maxdist){
  std::vector<Segment> temp;

  //temp.resize(120e6);


  //Load all the unit's points into a vector
  for(const auto &poly: mp)
  for(const auto &ring: poly)
  for(unsigned int ri=0;ri<ring.size();ri++){
    const auto &a   = ring.at(ri);
    const auto &b   = ring.at((ri+1)%ring.size()); //Loop around to beginning
    const auto dist = EuclideanDistance(a,b);

    //Segment is already less than maximum distance, so accept it without
    //modification
    if(dist<=maxdist){
      temp.emplace_back(a,b);
      continue;
    }

    //Number of steps needed to get from a to b, round down since the last step
    //will be handled specially to prevent floating-point error
    const int steps = std::floor(dist/maxdist);

    //Size of the steps as a fraction of the total distance
    const double stepsize = maxdist/dist;

    //Previous point
    Point2D pta = a;

    //Start at step 1, since starting at 0 would give a segment from A to A. Use
    //inclusive bounds (<=) since we rounded down above
    for(int i=1;i<=steps;i++){
      const double st = i*stepsize;

      Point2D ptb(
        (1-st)*a.x + st*b.x,
        (1-st)*a.y + st*b.y
      );   

      temp.emplace_back(pta,ptb);
      pta = ptb;
    }

    temp.emplace_back(pta,b);

  }
  
  return temp;
}



TEST_CASE("Densified border segments"){
  const std::string rect2by2 = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[0,0],[2,0],[2,2],[0,2],[0,0]]]}}]}";
  
  const auto gc   = ReadGeoJSON(rect2by2);
  const auto segs = GetDensifiedBorderSegments(gc.at(0), 0.03);

  for(const auto &seg: segs)
    REQUIRE(EuclideanDistance(seg.first,seg.second)<=0.03000001);
}





///Our strategy here is to build an R*-tree (SpIndex) out of the bounding boxes
///of the units. We expand the bounding boxes by \p expand_bb_by to ensure that
///those belonging to neighbouring units overlap. We then discretize the
///perimeters of the units and use the discretized perimeters to check for
///proximity.
void FindNeighbouringUnits(
  GeoCollection &gc,  
  const double max_neighbour_dist,     ///< Distance within which a units are considered to be neighbours.
  const double expand_bb_by            ///< Distance by which units' bounding boxes are expanded. Only units with overlapping boxes are checked for neighbourness. Value should be >0.
){
  //Since a neighbour relationship is two-way, we can save time by making notes
  //of those units for which we have already identified neighbours and avoiding
  //repeating that work
  std::unordered_set<unsigned int> found_neighbours;

  //Add all of the units to the R*-tree so we can quickly find neighbours.
  //Expand the bounding boxes of the units so that they will overlap if they are
  //neighbours
  std::cerr<<"Creating R*-tree..."<<std::endl;
  SpIndex<double,int> gcidx;
  for(unsigned int gci=0;gci<gc.size();gci++)
    AddToSpIndex(gc.at(gci), gcidx, gci, expand_bb_by);
  gcidx.buildIndex();

  std::cerr<<"Finding neighbours..."<<std::endl;
  //Loop through all of the units
  #pragma omp parallel for default(none) shared(gc,gcidx)
  for(unsigned int gci=0;gci<gc.size();gci++){
    auto &unit = gc.at(gci);

    const auto this_bb = unit.bbox();

    //Segment grid for quick minimal-distance calculation
    SegmentGrid sg(
      this_bb.xmin(),
      this_bb.ymin(),
      this_bb.xmax(),
      this_bb.ymax(),
      1.05*max_neighbour_dist, //This padding accounts for floating-point issues
      max_neighbour_dist      
    );

    //Add points from one unit
    for(auto &seg: GetDensifiedBorderSegments(unit, max_neighbour_dist))
      sg.addSegment(seg);    

    //Find the neighbours of the unit by overlapping bounding boxes. (This is
    //thread-safe!)
    const auto neighbours = gcidx.query(unit);

    //std::cerr<<"Determining neighbourness..."<<std::endl;
    //Loop over the neighbouring units
    for(const auto &n: neighbours){
      auto &nunit = gc.at(n);

      //If we already know this unit is our neighbour, skip it (reduces workload
      //by about one-half).
      if(std::find(nunit.neighbours.begin(), nunit.neighbours.end(), gci)!=nunit.neighbours.end())
        continue;

      const auto n_segs = GetDensifiedBorderSegments(nunit, max_neighbour_dist);
      for(auto &seg: n_segs){
        if(sg.isThereANearbySegment(seg)){
          //Both of the following lins must be critical since this unit could be
          //accessed by a neighbour or access a neighbour at the same time that
          //memory is being accessed by another thread.
          #pragma omp critical 
          {
            unit.neighbours.push_back(n);
            nunit.neighbours.push_back(gci);
          }
          break;
        }
      }

      //std::cout<<"Blue="<<sg.sizeBlue()<<", Red="<<sg.sizeRed()<<", total="<<(sg.sizeBlue()+sg.sizeRed())<<", comparisons="<<sg.numSegSegComparisons()<<std::endl;
    }
  }

  for(auto &unit: gc){
    unit.props["NEIGHNUM"]   = std::to_string(unit.neighbours.size());
    unit.props["NEIGHBOURS"] = "";
    for(const auto &n: unit.neighbours)
      unit.props["NEIGHBOURS"] += std::to_string(n) + ",";
    if(unit.neighbours.size()>0)
      unit.props["NEIGHBOURS"].pop_back();
  }
}






//Finds external children who still have 100% overlap with their parent
void FindExternalChildren(
  GeoCollection &subunits,
  GeoCollection &superunits,
  const double edge_adjacency_dist        ///< Distance within which a subunit is considered to be on the border of a superunit.
){
  //External children are those who are within epsilon of the borders of their
  //parents. We will find the external children now.

  for(auto &sub: subunits)
    sub.props["EXTCHILD"] = "F";

  //First, we build segment grids for the parents
  std::vector<SegmentGrid> sgvec(superunits.size());

  #pragma omp parallel for default(none) shared(superunits,sgvec)
  for(unsigned int supi=0;supi<superunits.size();supi++){
    const auto &unit   = superunits.at(supi);
    const auto this_bb = unit.bbox();

    //Create a segment grid for quick minimal-distance calculation
    sgvec.at(supi) = SegmentGrid(
      this_bb.xmin(),
      this_bb.ymin(),
      this_bb.xmax(),
      this_bb.ymax(),
      1.05*edge_adjacency_dist, //This padding accounts for floating-point issues
      edge_adjacency_dist      
    );

    //Add points from the superunit into the segment grid
    for(auto &seg: GetDensifiedBorderSegments(unit, edge_adjacency_dist))
      sgvec.at(supi).addSegment(seg);    
  }



  //Now we loop through the subunits to see if they are near the edges of their parents
  #pragma omp parallel for default(none) shared(subunits,sgvec)
  for(unsigned int subi=0;subi<subunits.size();subi++){
    //Alias the current subunit
    auto &sub = subunits[subi];

    //We have already determined that this is an external child, probably
    //through area considerations
    if(sub.props["EXTCHILD"]=="T")
      continue;

    const auto sub_border = GetDensifiedBorderSegments(sub, edge_adjacency_dist);

    bool external = false;

    //Loop through all of the parents of the subunit and, for each parent, see
    //if any of the subunit's edges are near that parent's edges
    for(const auto supp: sub.parents)
    for(const auto &seg: sub_border){
      if(sgvec.at(supp.first).isThereANearbySegment(seg)){
        external = true;
        goto escape_extchild_search;
      }
    }

    escape_extchild_search:
    sub.props["EXTCHILD"] = external?"T":"F";
  }
}



///Our strategy here is to build an R*-tree (SpIndex) out of the bounding boxes
///of the superunits. We then check each potential subunit to find the
///superunits it overlaps. We then discretize the perimeters of the units and
///use the discretized perimeters to check for proximity. If a district is
///proximal to the edge of a parent, we can run a more expensive polygon-polygon
///intersection calculation. Otherwise, we choose an arbitrary point of the
///subunit and perform point-in-polygon checks until we find a parent. (There is
///more than one PIP check since the bounding box of the subunit may intersect
///multiple potential parents even if the subunit is only within one parent.)
void CalcParentOverlap(
  GeoCollection &subunits,
  GeoCollection &superunits,
  const double complete_inclusion_thresh, ///< A subunit with more fractional area than this in the parent are 100% included, all other potential parents are ignored
  const double not_included_thresh,       ///< A subunit with less fractional area than this in a parent disregards that parent
  const double edge_adjacency_dist,       ///< Distance within which a subunit is considered to be on the border of a superunit.  
  const bool   print_parent_columns       ///< If True, a column is printed for each possible parent indicating subunit inclusion fractions. Generally you should want this to be False.
){
  SpIndex<double,int> supidx;

  if(superunits.empty())
    throw std::runtime_error("No superunits provided!");

  //Clear previous data before continuing
  for(auto &sub: subunits){
    sub.parents.clear();
    sub.props["PARENTNUM"] = "";
    sub.props["PARENTS"]   = "";
    sub.props["PARENTPR"]  = "";
    sub.props["CENTROIDX"] = "";
    sub.props["CENTROIDY"] = "";
  }

  for(auto &sup: superunits){
    sup.children.clear();
    sup.props["CHILDNUM"]   = "";
    sup.props["CHILDREN"]   = "";
    sup.props["CHILDRENPR"] = "";
  }


  //Add all the superunits to an R*-tree so we can quickly find potential
  //children using minimum bounding boxes.
  for(unsigned int i=0;i<superunits.size();i++)
    AddToSpIndex(superunits.at(i), supidx, i, 10.);
  supidx.buildIndex();


  //Use index to find potential parents. Use potential parents to find external
  //children.
  std::cout<<"Finding potential parents..."<<std::endl;
  #pragma omp parallel for default(none) shared(subunits,superunits,supidx)
  for(unsigned int subi=0;subi<subunits.size();subi++){
    auto &sub = subunits.at(subi);

    //Mark potential parents
    sub.parents.clear();
    for(const auto pp: supidx.query(sub)) //Query the R* tree to find which superunits the subunit overlaps
      sub.parents.emplace_back(pp,qnan);
  }


  //Determine whether each child is an edge child
  std::cout<<"Finding external children..."<<std::endl;
  FindExternalChildren(
    subunits,
    superunits,
    edge_adjacency_dist
  );


  std::cout<<"Calculating parent overlaps..."<<std::endl;
  #pragma omp parallel for default(none) shared(subunits,superunits,std::cerr)
  for(unsigned int i=0;i<subunits.size();i++){
    auto &sub = subunits.at(i);

    if(sub.props.at("EXTCHILD")=="T"){
      const double area  = areaExcludingHoles(sub);

      const auto potential_parents = sub.parents;
      sub.parents.clear();

      //Loop over the parent units
      for(auto &pp: potential_parents){
        const auto p = pp.first; //Parent id
        const double iarea = IntersectionArea(sub, superunits.at(p));
        const double frac  = iarea/area;
        if(frac>complete_inclusion_thresh){
          sub.parents.clear();
          sub.parents.emplace_back(p, 1);
          break;
        } else if(frac>not_included_thresh){
          sub.parents.emplace_back(p, frac);
        }
      }
    } else {
      //The unit is not an exterior child, so it lies entirely within a
      //superunit and is not near that superunit's border. However, its bounding
      //box may still intersect more than one superunit. Therefore, we choose an
      //arbitrary point and check to see which superunit it lies in.

      const auto potential_parents = sub.parents;
      sub.parents.clear();

      //The easy case
      if(potential_parents.size()==1){
        sub.parents.emplace_back(potential_parents.front().first, 1);
        continue;
      }

      if(sub.v.empty())
        throw std::runtime_error("Subunit has no polygons!");
      if(sub.v.front().empty())
        throw std::runtime_error("Subunit has no rings!");
      if(sub.v.front().v.empty())
        throw std::runtime_error("Subunit has no rings!");

      //Get first point of first ring of first polygon
      const auto first_pt = sub.v.front().v.front().v.front();

      //Loop over the parent units
      for(auto &pp: potential_parents){
        const auto p = pp.first; //Parent id
        if(ContainsPoint(superunits.at(p), first_pt)){
          sub.parents.emplace_back(p, 1);
          break;
        }
      }      

      #ifdef COMPACTNESSLIB_WARNINGS
        if(sub.parents.size()==0){
          #pragma omp critical
          {
            std::cerr<<"Warning! Could not find a containing parent for a subunit! Dumping properties."<<std::endl;
            for(const auto &kv: sub.props)
              std::cerr<<"\t"<<kv.first<<" = "<<kv.second<<std::endl;
          }
        }
      #endif
    }
  }

  //Make sure that the parents all have the same information as the children so
  //we can access it from either direction. Don't parallelize this loop since
  //multiple subunits could be accessing parents' memory at the same time.
  for(unsigned int i=0;i<subunits.size();i++){
    for(auto &p: subunits[i].parents)
      superunits.at(p.first).children.emplace_back(i, p.second);
  }


  // for(auto &sub: subunits)
  //   sub.props["EXTCHILD"] = "F";



  //Stringify the information as properties suitable for output
  for(auto &sub: subunits){
    sub.props["PARENTNUM"] = std::to_string(sub.parents.size());
    for(const auto &p: sub.parents){
      sub.props["PARENTS"]  += std::to_string(p.first)  + ",";
      sub.props["PARENTPR"] += std::to_string(p.second) + ",";
    }

    //TODO: Nuclear option for testing
    if(print_parent_columns){
      for(unsigned int i=0;i<superunits.size();i++)
        sub.props["PAR" + std::to_string(i)] = std::to_string(0);
      for(const auto &p: sub.parents)
        sub.props["PAR" + std::to_string(p.first)] = std::to_string(p.second);
    }

    //Drop the trailing commas
    if(sub.parents.size()>0){
      sub.props["PARENTS"].pop_back();
      sub.props["PARENTPR"].pop_back();
    }
    const auto centroid = CentroidPTSH(sub);
    sub.props["CENTROIDX"] = std::to_string(centroid.x);
    sub.props["CENTROIDY"] = std::to_string(centroid.y);
  }

  for(auto &sup: superunits){
    sup.props["CHILDNUM"]   = std::to_string(sup.children.size());
    sup.props["CHILDREN"]   = "";
    sup.props["CHILDRENPR"] = "";
    for(const auto &c: sup.children){
      sup.props["CHILDREN"]   += std::to_string(c.first)  + ",";
      sup.props["CHILDRENPR"] += std::to_string(c.second) + ",";
    }

    //Drop the trailing commas
    if(sup.children.size()>0){
      sup.props["CHILDREN"].pop_back();
      sup.props["CHILDRENPR"].pop_back();
    }
  }
}



}

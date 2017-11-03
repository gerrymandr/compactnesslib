#include "neighbours.hpp"
#include "geom.hpp"
#include "lib/nanoflann.hpp"
#include "SpIndex.hpp"
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace complib {


typedef std::vector<double>               pointdimvec_t;
typedef std::vector<pointdimvec_t>        pointvec_t;
typedef std::vector<unsigned int>         ownervec_t;
typedef std::pair<ownervec_t, pointvec_t> owner_point_vec_t;


class SegmentGrid {
 private:
  std::vector<Segment> segments;
  typedef std::vector<Segment*> segrefvec;

  class Cell {
   public:
    segrefvec red;               //List of red segments which overlap cell
    segrefvec blue;              //List of blue segments which overlap cell
  };

  //The grid. We assume that the data is sparse, since we are dealing with
  //polygons, and, therefore, use an unordered_map
  std::unordered_map<int, Cell> cells;

  double xmin;                   //Minimum x-coordinate of grid
  double ymin;                   //Minimum y-coordinate of grid
  double xmax;                   //Maximum x-coordinate of grid
  double ymax;                   //Maximum x-coordinate of grid
  double cellsize;               //Width & Height of a grid cell
  int    width;                  //Width of grid measured in cells
  int    height;                 //Height of grid measured in cells
  mutable int segcomp_count = 0; //Number of segment-segment comparisons

  bool doTheseVecsHaveSegmentInEpsilon(
    const segrefvec &veca,
    const segrefvec &vecb,
    const double epsilon
  ) const {
    if(veca.size()==0 || vecb.size()==0)
      return false;

    for(const auto &a: veca)
    for(const auto &b: vecb){
      segcomp_count++;
      const auto dist = SegmentSegmentDistanceSquared(
        a->first,
        a->second,
        b->first,
        b->second
      );
      if(dist<=epsilon)
        return true;
    }

    return false;
  }

 public:
  SegmentGrid() = default;

  SegmentGrid(const double xmin0, const double ymin0, const double xmax0, const double ymax0, const double cellsize0){
    xmin     = xmin0;
    ymin     = ymin0;
    xmax     = xmax0;
    ymax     = ymax0;
    cellsize = cellsize0;

    //Round up! Otherwise a segment end-point may fall off the grid.
    width    = std::ceil((xmax-xmin)/cellsize);
    height   = std::ceil((ymax-ymin)/cellsize);
  }

  void addSegment(
    const Point2D &a,  //End point of segment
    const Point2D &b,  //End point of segment
    const bool red     //True => red; otherwise => blue
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
    Cell& acell = cells.emplace(acelli, Cell()).first->second; //.first is the map object (<int,Cell>), .second is the cell
    Cell& bcell = cells.emplace(bcelli, Cell()).first->second; //.first is the map object (<int,Cell>), .second is the cell

    segments.emplace_back(a,b);

    //Make note of which colours are present in the cell
    if(red){
      acell.red.emplace_back ( &segments[segments.size()-1] );  
      bcell.red.emplace_back ( &segments[segments.size()-1] );  
    } else {
      acell.blue.emplace_back( &segments[segments.size()-1] ); 
      bcell.blue.emplace_back( &segments[segments.size()-1] ); 
    }
  }

  bool areThereSegmentsWithinEpsilon(const double epsilon) const {
    const int dx[9] = {0, -1, -1,  0,  1, 1, 1, 0, -1}; // x offsets of D8 neighbours, from a central cell
    const int dy[9] = {0,  0, -1, -1, -1, 0, 1, 1,  1}; // y offsets of D8 neighbours, from a central cell

    //Iterate through cells
    for(const auto &kv: cells){
      const auto idx    = kv.first;                 //Index of current central cell
      const auto& ccell = kv.second;                //The current central cell

      const auto cy = idx/height;                   //y-coordinate of central cell
      const auto cx = idx%height;                   //x-coordinate of central cell

      //Loop through neighbours
      for(int n=0;n<=8;n++){
        const auto nx = cx+dx[n];                   //x-coordinate of neighbour
        const auto ny = cy+dy[n];                   //y-coordinate of neighbour

        if(nx<0 || ny<0 || nx==width || ny==height) //Is cell outside grid bounds?
          continue;                                 //Yep

        const int nidx = ny*width+nx;               //Index of neighbour cell

        if(!cells.count(nidx))                      //Does this cell exist?
          continue;                                 //Nope

        const auto& ncell = cells.at(nidx);

        if(doTheseVecsHaveSegmentInEpsilon(ccell.red,  ncell.blue, epsilon))
          return true;
        if(doTheseVecsHaveSegmentInEpsilon(ccell.blue, ncell.red,  epsilon))
          return true;
      }
    }

    return false;
  }

  int sizeBlue() const {
    int sum = 0;
    for(const auto &kv: cells)
      sum += kv.second.blue.size();
    return sum;
  }

  int sizeRed() const {
    int sum = 0;
    for(const auto &kv: cells)
      sum += kv.second.red.size();
    return sum;
  }

  int numSegSegComparisons() const {
    return segcomp_count;
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



std::vector<Segment> GetDensifiedBorderSegments(const GeoCollection &gc, const double maxdist){
  std::vector<Segment> temp;

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

      Point2D pta = a;
      do {
        Point2D ptb(
          (1-st)*a.x + st*b.x,
          (1-st)*a.y + st*b.y
        );
        Segment(pta,ptb)
        si++;
        st = si*step;
      } while (st<1);

      temp.first.emplace_back(mpi);
    }
  }

  return temp;
}



///////////////////////////////////////////////////////////
//The following is copied verbatim from a nanoflann example

/** A simple vector-of-vectors adaptor for nanoflann, without duplicating the storage.
  *  The i'th vector represents a point in the state space.
  *
  *  \tparam DIM If set to >0, it specifies a compile-time fixed dimensionality for the points in the data set, allowing more compiler optimizations.
  *  \tparam num_t The type of the point coordinates (typically, double or float).
  *  \tparam Distance The distance metric to use: nanoflann::metric_L1, nanoflann::metric_L2, nanoflann::metric_L2_Simple, etc.
  *  \tparam IndexType The type for indices in the KD-tree index (typically, size_t of int)
  */
template <class VectorOfVectorsType, typename num_t = double, int DIM = -1, class Distance = nanoflann::metric_L2, typename IndexType = size_t>
struct KDTreeVectorOfVectorsAdaptor {
  typedef KDTreeVectorOfVectorsAdaptor<VectorOfVectorsType,num_t,DIM,Distance> self_t;
  typedef typename Distance::template traits<num_t,self_t>::distance_t metric_t;
  typedef nanoflann::KDTreeSingleIndexAdaptor< metric_t,self_t,DIM,IndexType>  index_t;

  index_t* index; //! The kd-tree index for the user to call its methods as usual with any other FLANN index.

  /// Constructor: takes a const ref to the vector of vectors object with the data points
  KDTreeVectorOfVectorsAdaptor(const int dimensionality, const VectorOfVectorsType &mat, const int leaf_max_size = 10) : m_data(mat){
    assert(mat.size() != 0 && mat[0].size() != 0);
    const size_t dims = mat[0].size();
    if (DIM>0 && static_cast<int>(dims) != DIM)
      throw std::runtime_error("Data set dimensionality does not match the 'DIM' template argument");
    index = new index_t( dims, *this /* adaptor */, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size ) );
    index->buildIndex();
  }

  ~KDTreeVectorOfVectorsAdaptor() {
    delete index;
  }

  const VectorOfVectorsType &m_data;

  /** Query for the \a num_closest closest points to a given point (entered as query_point[0:dim-1]).
    *  Note that this is a short-cut method for index->findNeighbors().
    *  The user can also call index->... methods as desired.
    * \note nChecks_IGNORED is ignored but kept for compatibility with the original FLANN interface.
    */
  inline void query(const num_t *query_point, const size_t num_closest, IndexType *out_indices, num_t *out_distances_sq, const int nChecks_IGNORED = 10) const
  {
    nanoflann::KNNResultSet<num_t,IndexType> resultSet(num_closest);
    resultSet.init(out_indices, out_distances_sq);
    index->findNeighbors(resultSet, query_point, nanoflann::SearchParams());
  }

  /** @name Interface expected by KDTreeSingleIndexAdaptor
    * @{ */

  const self_t & derived() const {
    return *this;
  }
  self_t & derived()       {
    return *this;
  }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const {
    return m_data.size();
  }

  // Returns the dim'th component of the idx'th point in the class:
  inline num_t kdtree_get_pt(const size_t idx, int dim) const {
    return m_data[idx][dim];
  }

  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX & /*bb*/) const {
    return false;
  }
};



//////////////////////////////////////////////////////////////////
//Okay, back to my stuff

void FindNeighbouringDistricts(
  GeoCollection &gc,  
  const double max_neighbour_pt_dist,     ///< Distance within which a units are considered to be neighbours.
  const double max_boundary_pt_dist,      ///< Maximum distance between points on densified boundaries.
  const double expand_bb_by               ///< Distance by which units' bounding boxes are expanded. Only districts with overlapping boxes are checked for neighbourness. Value should be >0.
){
  //Since a neighbour relationship is two-way, we can save time by making notes
  //of those units for which we have already identified neighbours and avoiding
  //repeating that work
  std::unordered_set<unsigned int> found_neighbours;

  //Add all of the units to the R*-tree so we can quickly find neighbours.
  //Expand the bounding boxes of the units so that they will overlap if they are
  //neighbours
  std::cerr<<"Creating R*-tree..."<<std::endl;
  SpIndex gcidx;
  for(unsigned int i=0;i<gc.size();i++)
    AddToSpIndex(gc.at(i), gcidx, i, expand_bb_by);
  gcidx.buildIndex();

  //Loop through all of the units
  #pragma omp parallel for
  for(unsigned int gci=0;gci<gc.size();gci++){
    auto &unit = gc.at(gci);

    const auto this_bb = unit.bbox();

    // found_neighbours.insert(i); //TODO

    //Find the neighbours of the unit by overlapping bounding boxes
    const auto neighbours = gcidx.query(unit);

    std::cerr<<"Determining neighbourness..."<<std::endl;
    //Loop over the neighbouring units
    for(const auto &n: neighbours){
      auto &nunit = gc.at(n);

      //We've already determined the neighbour relationships for unit `n`, so
      //skip it.
      //TODO

      const auto bbn = nunit.bbox();

      const auto combined_bbox = this_bb+bbn;

      //Segment grid for quick minimal-distance calculation
      SegmentGrid sg(
        combined_bbox.xmin(),
        combined_bbox.ymin(),
        combined_bbox.xmax(),
        combined_bbox.ymax(),
        max_neighbour_pt_dist        
      );

      //Lambda to add points to the segment grid
      const auto add_points = [&](const MultiPolygon &mp, const bool red){
        //Loop over points of the superunit
        for(const auto &poly1: mp)
        for(const auto &ring1: poly1)
        for(unsigned int i1=0;i1<ring1.size();i1++)
          sg.addSegment( ring1.at(i1), ring1.at((i1+1)%ring1.size()), red);
      };

      //Add points from both MultiPolygons. Arbitrarily call one red and the
      //other blue
      add_points(unit,  true );
      add_points(nunit, false);

      if(sg.areThereSegmentsWithinEpsilon(max_neighbour_pt_dist))
        unit.neighbours.push_back(n);

      std::cout<<"Blue="<<sg.sizeBlue()<<", Red="<<sg.sizeRed()<<", total="<<(sg.sizeBlue()+sg.sizeRed())<<", comparisons="<<sg.numSegSegComparisons()<<std::endl;
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



void CalcParentOverlap(
  GeoCollection &subunits,
  GeoCollection &superunits,
  const double complete_inclusion_thresh, ///< A subunit with more fractional area than this in the parent are 100% included, all other potential parents are ignored
  const double not_included_thresh,       ///< A subunit with less fractional area than this in a parent disregards that parent
  const double max_boundary_pt_dist,      ///< Maximum distance between points on densified boundaries.
  const double edge_adjacency_dist        ///< Distance within which a subunit is considered to be on the border of a superunit.
){
  SpIndex supidx;

  //Add all the superunits to an R*-tree so we can quickly find potential
  //children using minimum bounding boxes.
  for(unsigned int i=0;i<superunits.size();i++)
    AddToSpIndex(superunits.at(i), supidx, i, 0.);
  supidx.buildIndex();

  #pragma omp parallel for
  for(unsigned int i=0;i<subunits.size();i++){
    auto &sub = subunits.at(i);

    const auto parents = supidx.query(sub);
    const double area  = areaExcludingHoles(sub);

    //Loop over the parent units
    for(auto &p: parents){
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
  }

  //Make sure that the parents all have the same information as the children so
  //we can access it from either direction
  for(unsigned int i=0;i<subunits.size();i++){
    for(auto &p: subunits[i].parents)
      superunits.at(p.first).children.emplace_back(i, p.second);
  }

  for(auto &sup: superunits)
    sup.props["CHILDREN"] = std::to_string(sup.children.size());

  for(auto &sub: subunits)
    sub.props["EXTCHILD"] = "F";


  //The algorithm relies on boundary points having a certain maximum spacing.
  //Ensure that the boundaries meet this requirement.
  const auto sup_densified_borders = GetDensifiedBorders(superunits, max_boundary_pt_dist);
  const auto sub_densified_borders = GetDensifiedBorders(subunits,   max_boundary_pt_dist);

  //kdtree that will be used to provide rapid nearest neighbour lookups among
  //the points which belong to the borders of the superunits
  typedef KDTreeVectorOfVectorsAdaptor< pointvec_t, double >  my_kd_tree_t;
  my_kd_tree_t border_idx(2 /*dim*/, sup_densified_borders.second, 10 /* max leaf */ );
  std::cerr<<"Constructing kd-tree..."<<std::endl;
  border_idx.index->buildIndex();
  std::cerr<<"done."<<std::endl;

  #pragma omp parallel for default(none) shared(subunits, border_idx)
  for(unsigned int subi=0;subi<subunits.size();subi++){
    //Alias the current subunit
    auto &sub = subunits[subi];

    //Loop through all of the points in the subunit
    for(const auto &pt: sub_densified_borders.second){
      //Find the nearest neighbour to the query point, irrespective of distance
      const size_t num_results = 1; //Number of nearest neighbours to find
      size_t nn_index;              //Index of the point that's been found
      double nn_dist_sqr;           //Squared distance to the point
      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&nn_index, &nn_dist_sqr);
      border_idx.index->findNeighbors(resultSet, pt.data(), nanoflann::SearchParams(10));

      if(nn_dist_sqr<edge_adjacency_dist*edge_adjacency_dist){
        sub.props["EXTCHILD"] = "T";
        goto escape_extchild_search;
      }
    }

    escape_extchild_search:
    (void)1;
  }



  for(auto &sub: subunits){
    sub.props["PARENTNUM"] = std::to_string(sub.parents.size());
    for(const auto &p: sub.parents){
      sub.props["PARENTS"]  += std::to_string(p.first)  + ",";
      sub.props["PARENTPR"] += std::to_string(p.second) + ",";
    }
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
    if(sup.children.size()>0){
      sup.props["CHILDREN"].pop_back();
      sup.props["CHILDRENPR"].pop_back();
    }
  }
}




}

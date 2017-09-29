#include "neighbours.hpp"
#include "geom.hpp"
#include "lib/nanoflann.hpp"
#include "SpIndex.hpp"
#include <vector>


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


namespace complib {

void FindNeighbouringDistricts(GeoCollection &gc){
  const double max_neighbour_pt_dist_squared = 1000*1000; //in metres

  //Point along the borders
  typedef std::vector< std::vector<double> > my_vector_of_vectors_t;
  std::vector< std::vector<double> > borders;
  //Which MultiPolygon each point came from
  std::vector<int> owner;

  //Load all the outer borders of the geometries into the borders point cloud
  for(unsigned int mpi=0;mpi<gc.size();mpi++){
    for(const auto &poly: gc.at(mpi))
    for(const auto &pt: poly.at(0)){
      borders.push_back(std::vector<double>({{pt.x,pt.y}}));
      owner.push_back(mpi);
    }
  }

  typedef KDTreeVectorOfVectorsAdaptor< my_vector_of_vectors_t, double >  my_kd_tree_t;
  my_kd_tree_t border_idx(2 /*dim*/, borders, 10 /* max leaf */ );
  border_idx.index->buildIndex();

  //Determine which districts border each other
  for(unsigned int i=0;i<borders.size();i++){
    nanoflann::SearchParams params;
    std::vector<std::pair<size_t, double> > temp;
    border_idx.index->radiusSearch(borders[i].data(), max_neighbour_pt_dist_squared, temp, params);
    for(const auto &sr: temp){
      if(owner[sr.first]==owner[i])
        continue;
      gc[owner[sr.first]].neighbours.insert(owner[i]);
      gc[owner[i]].neighbours.insert(owner[sr.first]);
    }
  }
}




void FindExteriorDistricts(
  GeoCollection &subunits,
  const GeoCollection &superunits,
  const int shrink,                //Amount by which to shrink superunits when trying to determine their children
  double border_dist_cutoff        //Subunits with at least one point within this distance of a superunit boundary are said to be part of the exterior set of that superunit
){
  //Square `border_dist_cutoff` because we are dealing with squared distances
  border_dist_cutoff *= border_dist_cutoff;

  //Data type for storing border information: see below for details.
  typedef std::vector< std::vector<double> > my_vector_of_vectors_t;
  //Vector which will contain all of the border points of the superunits. Later,
  //we'll dump them into a kd-tree. This allows us to rapidly build an index.
  std::vector< std::vector<double> > borders; //TODO: Data type
  //Which MultiPolygon each point came from
  std::vector<int> owner;

  //Shrink districts so that if two boundaries coalign there will be some
  //separation between the boundaries. This means that if a subunit is fully-
  //enclosed the nearest point to it will always be its superunit. If a subunit
  //is on the exterior, then the nearest points will (likely) be split between
  //two superunits. The exception is if a subunit shares a boundary with its
  //superunit. In this case the subunit may still have all its nearest points
  //inside of one superunit. However, we can say that any points on a subunit
  //sufficiently close to a superunit point indicate that the border is shared.
  #pragma omp parallel for
  for(unsigned int mpi=0;mpi<superunits.size();mpi++){
    //Shrink the superunit by `shrink` units
    const auto shrunk_superunit = BufferPath(superunits.at(mpi).clipper_paths, -shrink);

    //const auto shrunk_superunit = superunits.at(mpi).clipper_paths;

    //Loop through all of the points of the superunit
    #pragma omp critical
    for(const auto &path: shrunk_superunit)
    for(const auto &pt: path){
      //Add this point to the vector of border points
      borders.push_back(std::vector<double>({{pt.X,pt.Y}})); //TODO: Data type
      owner.push_back(mpi);
    }
  }

  //kdtree that will be used to provide rapid nearest neighbour lookups among
  //the points which belong to the borders of the superunits
  typedef KDTreeVectorOfVectorsAdaptor< my_vector_of_vectors_t, double >  my_kd_tree_t;
  my_kd_tree_t border_idx(2 /*dim*/, borders, 10 /* max leaf */ );
  border_idx.index->buildIndex();

  #pragma omp parallel for
  for(unsigned int subi=0;subi<subunits.size();subi++){
    //Alias the current subunit
    auto &sub = subunits[subi]; 

    //Holds parents so we can insert them multiple times
    std::set<int> parents; 
    //Holds onto exterior so we don't have to repeatedly access map
    bool exterior = false; 

    //Loop through all of the points in the subunit
    for(const auto &poly: sub.v)
    for(const auto &ring: poly)
    for(const auto &pt: ring){
      //One of the boundary points of a subunit
      const double qp[2] = {pt.x,pt.y};

      const size_t num_results = 1; //Number of nearest neighbours to find
      size_t ret_index;             //Index of the  point that's been found
      double out_dist_sqr;          //Squared distance to the point
      //Get the nearest neighbour
      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&ret_index, &out_dist_sqr);
      border_idx.index->findNeighbors(resultSet, &qp[0], nanoflann::SearchParams(10));

      //If the squared distance from the subunit point to its nearest neighbour
      //superunit point is less than the square distance cutoff, then we the
      //subunit as being exterior
      if(out_dist_sqr<border_dist_cutoff)
        exterior = true;

      parents.insert(owner.at(ret_index));

      //While it seems like we could break early if we've found an exterior
      //child, this isn't so: we still want to know what other parents the child
      //might have
    }

    //Note that a subunit might be an exterior subunit but have only one parent
    //because its border is coincident with a superunit
    if(exterior){
      sub.props["EXTCHILD"] = "T";
    } else {
      //However, if a subunit is not an exterior unit, then it can, by
      //definition, have only one parent.
      assert(sub.parents.size()==1);
      sub.props["EXTCHILD"] = "F";
    }

    for(const auto &p: parents)
      sub.parents.emplace_back(p, -1);

    // #pragma omp critical
    // for(const auto &p: parents)
    //   superunits.at(p).children.emplace_back(subi, -1);
  }
}



void CalcParentOverlap(GeoCollection &subunits, GeoCollection &superunits){
  const double not_parent_thresh = 0.997;

  for(unsigned int i=0;i<subunits.size();i++){
    auto &sub = subunits.at(i);
    //If the subunit is an interior unit (not an exterior unit) or has only one
    //parent, it must lie entirely in that parent
    if(sub.props["EXTCHILD"]=="F" || sub.parents.size()==1){
      assert(sub.parents.size()==1);
      sub.parents.at(0).second = 1;
      continue;
    }

    const double area = areaExcludingHoles(sub);

    //At this point, we suspect that some portion of the subunit lies outside of
    //the superunit parent. But, remember, we shrunk the superunit a bit, so
    //we'll need to be careful about attribution. Here, we intersect the
    //subunits with the whole area of the unshrunk superunits. Afterwards, we'll
    //check to see if we did any misattribution of parenthood.

    //Loop over the parent units
    for(auto &p: sub.parents){
      double iarea = IntersectionArea(sub, superunits.at(p.first));
      p.second     = iarea/area;
    }

    //Sort parents in order of decreasing percent area
    std::sort(
      sub.parents.begin(),
      sub.parents.end(),
      [](const MultiPolygon::parent_t &a, const MultiPolygon::parent_t &b){
        return a.second>b.second;
      }
    );

    //Running total of percent sums
    double contribution_sum = 0;
    //Pointer to the current superunit parent (starts with the parent within
    //which most of the subunit lies)
    auto p = sub.parents.begin();
    //Loop through all the parents
    for(;p!=sub.parents.end();p++){
      //Keep track of how much of the subunit's percent area is accounted for
      contribution_sum += p->second;
      //If enough of the subunit has been allocated, then we break, assuming the
      //rest can be chalked up to floating point errors
      if(contribution_sum>=not_parent_thresh)
        break;
    }
    //Eliminate the other parents
    sub.parents.erase(p,sub.parents.end());
    //Add the remainder to the final parent so that things add up closer to 1
    sub.parents.back().second += 1-contribution_sum;
  }

  //Make sure that the parents all have the same information as the children so
  //we can access it from either direction
  for(unsigned int i=0;i<subunits.size();i++)
  for(auto &p: subunits[i].parents)
    superunits.at(p.first).children.emplace_back(i, p.second);
}

}

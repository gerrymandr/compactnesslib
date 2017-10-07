#include "neighbours.hpp"
#include "geom.hpp"
#include "lib/nanoflann.hpp"
#include "SpIndex.hpp"
#include <vector>
#include <unordered_map>
#include <unordered_set>

//Loop through vector, remove elements not present in set
void VectorIntersect(std::vector<unsigned int> &vec, const std::unordered_set<unsigned int> &set){
  for(auto i=vec.begin();i!=vec.end();){
    if(set.count(*i)==0)
      i = vec.erase(i);
    else
      i++;
  }
}

//Loop through vector, remove elements not present in set
void VectorIntersect(std::vector<unsigned int> &vec, const unsigned int set){
  for(auto i=vec.begin();i!=vec.end();){
    if(*i!=set)
      i = vec.erase(i);
    else
      i++;
  }
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



void CalcParentOverlap(GeoCollection &subunits, GeoCollection &superunits){
  const double not_parent_thresh    = 0.997;
  const double max_boundary_pt_dist = 500;

  SpIndex<double, unsigned int> supidx;

  for(unsigned int i=0;i<superunits.size();i++)
    AddToSpIndex(superunits.at(i), supidx, i);
  supidx.buildIndex();


  for(unsigned int i=0;i<subunits.size();i++){
    auto &sub = subunits.at(i);

    const auto parents = supidx.query(sub);
    const double area  = areaExcludingHoles(sub);

    //Loop over the parent units
    for(auto &p: parents){
      const double iarea = IntersectionArea(sub, superunits.at(p));
      const double frac  = iarea/area;
      if(frac>0.997){
        sub.parents.clear();
        sub.parents.emplace_back(p, 1);
        break;
      } else if(frac>0.003){
        sub.parents.emplace_back(p, frac);
      }
    }
  }

  //Make sure that the parents all have the same information as the children so
  //we can access it from either direction
  for(unsigned int i=0;i<subunits.size();i++){
    auto &sub = subunits.at(i);
    sub.props["PARENTS"] = std::to_string(sub.parents.size());
    for(auto &p: subunits[i].parents)
      superunits.at(p.first).children.emplace_back(i, p.second);
  }

  for(auto &sup: superunits)
    sup.props["CHILDREN"] = std::to_string(sup.children.size());

  for(auto &sub: subunits)
    sub.props["EXTCHILD"] = "F";


  //The algorithm relies on boundary points having a certain maximum spacing.
  //Ensure that the boundaries meet this requirement.
  for(auto &sup: superunits)
    Densify(sup, max_boundary_pt_dist);

  for(auto &sub: subunits)
    Densify(sub, max_boundary_pt_dist);

  typedef std::vector< std::vector<double> > my_vector_of_vectors_t;

  //Vector which will contain all of the border points of the superunits. Later,
  //we'll dump them into a kd-tree. This allows us to rapidly build an index.
  std::vector< std::vector<double> > borders; //TODO: Data type
  //Which MultiPolygon each point came from
  std::vector<int> owner;

  //Load all the super unit points into a vector that will be loaded into a kd-
  //tree for quick nearest-neighbour lookups
  for(unsigned int mpi=0;mpi<superunits.size();mpi++){
    //Loop through all of the points of the superunit
    for(const auto &poly: superunits.at(mpi))
    for(const auto &ring: poly)
    for(const auto &pt: ring){
      //Add this point to the vector of border points
      borders.push_back(std::vector<double>({{pt.x,pt.y}}));
      owner.push_back(mpi);
    }
  }

  //kdtree that will be used to provide rapid nearest neighbour lookups among
  //the points which belong to the borders of the superunits
  typedef KDTreeVectorOfVectorsAdaptor< my_vector_of_vectors_t, double >  my_kd_tree_t;
  my_kd_tree_t border_idx(2 /*dim*/, borders, 10 /* max leaf */ );
  std::cerr<<"Constructing kd-tree..."<<std::endl;
  border_idx.index->buildIndex();
  std::cerr<<"done."<<std::endl;

  #pragma omp parallel for
  for(unsigned int subi=0;subi<subunits.size();subi++){
    //Alias the current subunit
    auto &sub = subunits[subi];

    //Loop through all of the points in the subunit
    for(const auto &poly: sub.v)
    for(const auto &ring: poly)
    for(const auto &pt: ring){
      //One of the boundary points of a subunit
      const double qp[2] = {pt.x,pt.y};

      //Find the nearest neighbour to the query point, irrespective of distance
      const size_t num_results = 1; //Number of nearest neighbours to find
      size_t nn_index;              //Index of the  point that's been found
      double nn_dist_sqr;           //Squared distance to the point
      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&nn_index, &nn_dist_sqr);
      border_idx.index->findNeighbors(resultSet, &qp[0], nanoflann::SearchParams(10));

      if(nn_dist_sqr<2000*2000){
        sub.props["EXTCHILD"] = "T";
        goto escape_extchild_search;
      }
    }

    escape_extchild_search:
    (void)1;
  }
}




}

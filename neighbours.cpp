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

void FindExteriorDistricts(GeoCollection &subunits, const GeoCollection &superunits){
  const double max_neighbour_pt_dist_squared = 1000*1000; //in metres

  //Point along the borders
  typedef std::vector< std::vector<double> > my_vector_of_vectors_t;
  std::vector< std::vector<double> > borders;
  //Which MultiPolygon each point came from
  std::vector<int> owner;

  //Load all the outer borders of the superunits into the point cloud
  for(unsigned int mpi=0;mpi<superunits.size();mpi++){
    for(const auto &poly: superunits.at(mpi))
    for(const auto &pt: poly.at(0)){
      borders.push_back(std::vector<double>({{pt.x,pt.y}}));
      owner.push_back(mpi);
    }
  }

  typedef KDTreeVectorOfVectorsAdaptor< my_vector_of_vectors_t, double >  my_kd_tree_t;
  my_kd_tree_t border_idx(2 /*dim*/, borders, 10 /* max leaf */ );
  border_idx.index->buildIndex();

  //Determine which districts are exterior districts
  for(auto &mp: subunits){
    for(const auto &poly: mp)
    for(const auto &ring: poly)
    for(const auto &pt: ring){
      const double qp[2] = {pt.x,pt.y};
      nanoflann::SearchParams params;
      std::vector<std::pair<size_t, double> > temp;
      border_idx.index->radiusSearch(qp, max_neighbour_pt_dist_squared, temp, params);
      if(temp.size()>0){
        //NOTE: Could make this all faster by promoting this to an actual
        //property, but then I would need to modify the output drivers.
        mp.props["EXTCHILD"] = "T"; 
        goto escape_the_multipolygon;
      }
    }
    escape_the_multipolygon:
    (void)1;
    if(!mp.props.count("EXTCHILD"))
      mp.props["EXTCHILD"] = "F";
  }
}



void CalcParentOverlap(GeoCollection &subunits, GeoCollection &superunits){
  //Make an Rtree!
  SpIndex<double, unsigned int> sp;

  //Add all superunits to the Rtree
  for(unsigned int sup=0;sup<superunits.size();sup++)
    AddToSpIndex(superunits.at(sup), sp, sup);
  sp.buildIndex();

  for(auto &sub: subunits){
    //Area of subunit
    const auto sub_area = areaIncludingHoles(sub);

    //All superunits which intersect the subunit
    const auto potential_parents = sp.query(sub);

    for(const auto &pp: potential_parents){
      auto ifrac = IntersectionArea(superunits.at(pp), sub)/sub_area;

      //Round to 1 if we're close enough
      if(ifrac>0.997)
        ifrac = 1;

      sub.parents.emplace_back(pp, ifrac);
    }
  }
}

}

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




void FindExteriorDistricts(
  GeoCollection &subunits,
  GeoCollection &superunits,
  const double max_boundary_pt_dist, //Maximum distance between boundary points
  const double border_dist_cutoff    //Points within at least this distance of a superunit boundary are potentially included in that superunit
){
  std::cerr<<"Subunit count   = "      <<subunits.size()       <<std::endl;
  std::cerr<<"Superunit count = "      <<superunits.size()     <<std::endl;
  std::cerr<<"Subunit point count   = "<<PointCount(subunits)  <<std::endl;
  std::cerr<<"Superunit point count = "<<PointCount(superunits)<<std::endl;

  //The algorithm relies on boundary points having a certain maximum spacing.
  //Ensure that the boundaries meet this requirement.
  for(auto &sup: superunits)
    Densify(sup, max_boundary_pt_dist);

  for(auto &sub: subunits)
    Densify(sub, max_boundary_pt_dist);

  //Data type for storing border information: see below for details.
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

    //Keep track of whether any parents are coincident on all boundary points.
    //To do so, we'll just store all the incoming parents in a long vector and
    //then get the vector's unique elements (if any) afterwards. This allows us
    //to take advantage of sequential memory and cuts computation in the loop.
    std::vector<unsigned int> bnd_parents_set;
    bool bnd_parents_init = false; //Indicates whether we've seen any boundary points yet
    //Keep track of which parents had how many child boundary points
    std::unordered_map<unsigned int, unsigned int> bnd_pts_per_parent; //<parent,point count> pairs
    //Holds parents so we can insert them multiple times
    std::unordered_set<unsigned int> parents;

    //Loop through all of the points in the subunit
    for(const auto &poly: sub.v)
    for(const auto &ring: poly)
    for(const auto &pt: ring){
      //One of the boundary points of a subunit
      const double qp[2] = {pt.x,pt.y};

      //Find the nearest neighbours with a radius of the query point
      std::vector< std::pair<size_t,double> > rs_matches; //<point index, distance pairs>
      nanoflann::SearchParams params;
      //Result is sorted by increasing distance. Distances are squared.
      border_idx.index->radiusSearch(&qp[0], border_dist_cutoff, rs_matches, params);

      //Find the nearest neighbour to the query point, irrespective of distance
      const size_t num_results = 1; //Number of nearest neighbours to find
      size_t nn_index;              //Index of the  point that's been found
      double nn_dist_Sqr;           //Squared distance to the point
      nanoflann::KNNResultSet<double> resultSet(num_results);
      resultSet.init(&nn_index, &nn_dist_Sqr);
      border_idx.index->findNeighbors(resultSet, &qp[0], nanoflann::SearchParams(10));

      //Parents as determined by radius search
      std::unordered_set<unsigned int> radius_parents;
      for(const auto &rn: rs_matches) //rn = radius neighbour
        radius_parents.insert(owner.at(rn.first));

      const unsigned int nn_parent = owner.at(nn_index);

      if(radius_parents.size()==0){
        //Point is unambiguously a member of one superunit

        parents.insert(nn_parent);

      } else if(radius_parents.size()==1){
        //Point is unambiguously a member of one superunit and on its border.

        //Since all members of the radius search belong to a single parent and
        //these points include the nearest neighbour, we can use `nn_parent` to
        //identify the parent.
        bnd_pts_per_parent[nn_parent]++;
        if(!bnd_parents_init){
          bnd_parents_set.push_back(nn_parent);
          bnd_parents_init=true;
        } else {
          VectorIntersect(bnd_parents_set, nn_parent);
        }


      } else if(radius_parents.size()>1) {
        //Point is on the border of more than one superunit. Make a note of
        //both.
        for(const auto &rn: radius_parents)
          bnd_pts_per_parent[rn]++;
        if(!bnd_parents_init){
          bnd_parents_set.insert(bnd_parents_set.end(), radius_parents.begin(), radius_parents.end());
          bnd_parents_init = true;
        } else {
          VectorIntersect(bnd_parents_set, radius_parents);
        }

      } else {
        throw std::logic_error("Execution should never reach this point!");
      }
    }

    //Initially mark subunit as not being on the exterior
    sub.props["EXTCHILD"] = "F";

    if(parents.size()==1){
      //Subunit unambiguously has only one parent. Even though the boundary
      //points might be associated with more than one parent, the ambiguity is
      //resolved by the interior points.
      sub.parents.emplace_back(*parents.begin(), 1); //Subunit is entirely in parent
      if(bnd_pts_per_parent.size()>0)
        sub.props["EXTCHILD"] = "T";
      std::cerr<<sub.props["GEOID"]<<" - One Parent"<<std::endl;
    } else if(parents.size()>1){
      //Subunit unambiguously has more than one parent. Even though the boundary
      //points might be associated with many parents, the ambiguity is resolved
      //by the interior points.
      for(const auto &p: parents)
        sub.parents.emplace_back(p, -1); //Fraction of area is unknown

      //To have more than one parent, the subunit must cross the boundaries of
      //the parents, so it is an exterior child
      sub.props["EXTCHILD"] = "T";

      std::cerr<<sub.props["GEOID"]<<" - Multiple parents"<<std::endl;
    } else if(parents.size()==0 && bnd_pts_per_parent.size()==0){
      //Uh oh, the subunit isn't associated with any superunit!
      throw std::runtime_error("Point found outside of superunit interior and boundary region!");
    } else if(parents.size()==0 && bnd_pts_per_parent.size()>0){

      //Subunit has no unambiguous parent, but has boundary points with many
      //parents. There are a couple of cases to distinguish between here:
      if(bnd_parents_set.size()==1){
        //1. The subunit boundary is everywhere coincident with the superunit
        //   boundary, but may contact multiple surrounding superunits.
        sub.parents.emplace_back(*bnd_parents_set.begin(), 1);
        sub.props["EXTCHILD"] = "T";
        std::cerr<<sub.props["GEOID"]<<" - Case 1: Everywhere coincident"<<std::endl;
      } else if(bnd_parents_set.size()==2){
        //2. The subunit boundary is everywhere coincident with the superunit
        //   boundary, but also with a single surrounding superunit (one 
        //   superunit contains another as a hole). We determine which the
        //   parent is by comparing the number of points of the superunits.
        //   Since the encompassing superunit has an outer boundary we haven't
        //   seen, we expect it to have many more points the the subunit
        //   boundary. This assumes that all boundaries have similar point
        //   densities - a property which should be guaranteed by the
        //   densification. This isn't foolproof, though!
        const auto subunit_pts = PointCount(sub);

        unsigned int min_pts    = std::numeric_limits<unsigned int>::max();
        unsigned int min_parent = std::numeric_limits<unsigned int>::max();

        //Find the potential parent whose point count most closely matches the
        //subunit point count
        for(const auto &sup: bnd_parents_set){
          const unsigned int sup_pts = PointCount(superunits.at(sup));
          const unsigned int pt_diff = std::abs((int)sup_pts-(int)subunit_pts);
          if(pt_diff<min_pts){
            min_pts    = pt_diff;
            min_parent = sup;
          }
        }
        sub.parents.emplace_back(min_parent, 1);
        sub.props["EXTCHILD"] = "T";
        std::cerr<<sub.props["GEOID"]<<" - Case 2: Hole"<<std::endl;
      } else if (bnd_parents_set.size()==0) { 
        //3. The subunit boundary is coincident to multiple superunit boundaries
        //   and no single superunit is present at all the subunit's boundary
        //   points. This can happen when border misalignment is large enough 
        //   that the subunit boundary meanders away from its actual superunit 
        //   but not into the interior of neighbouring superunits. In this case,
        //   we defer to the majority.
        unsigned int max_pts    = 0;
        unsigned int max_parent = std::numeric_limits<unsigned int>::max();
        for(const auto &kv: bnd_pts_per_parent){
          if(kv.second>max_pts){
            max_pts    = kv.second;
            max_parent = kv.first;
          }
        }
        sub.parents.emplace_back(max_parent, 1);
        sub.props["EXTCHILD"] = "T";
        std::cerr<<sub.props["GEOID"]<<" - Case 3: Multiple boundaries"<<std::endl;
      } else {
        std::cerr<<"parents size = "<<parents.size()<<std::endl;
        std::cerr<<"bnd_pts_per_parent size = "<<bnd_pts_per_parent.size()<<std::endl;
        std::cerr<<"bnd_parents_set size = "<<bnd_parents_set.size()<<std::endl;
        throw std::logic_error("Unexpected boundary point edge case! (Contact package maintainers, please.)");
      }
      //NOTE: Everything here is essentially an edge case, since it is expected
      //      that subunit-superunit coincidence will represent a minority of
      //      the cases encountered in actual data.
    } else {
      throw std::logic_error("Execution should never reach this point!");
    }

  }
}



// void FindExteriorDistricts(
//   GeoCollection &subunits,
//   const GeoCollection &superunits,
//   const int shrink,                //Amount by which to shrink superunits when trying to determine their children
//   double border_dist_cutoff        //Subunits with at least one point within this distance of a superunit boundary are said to be part of the exterior set of that superunit
// ){
//   //Square `border_dist_cutoff` because we are dealing with squared distances
//   border_dist_cutoff *= border_dist_cutoff;

//   //Data type for storing border information: see below for details.
//   typedef std::vector< std::vector<double> > my_vector_of_vectors_t;
//   //Vector which will contain all of the border points of the superunits. Later,
//   //we'll dump them into a kd-tree. This allows us to rapidly build an index.
//   std::vector< std::vector<double> > borders; //TODO: Data type
//   //Which MultiPolygon each point came from
//   std::vector<int> owner;

//   //Shrink districts so that if two boundaries coalign there will be some
//   //separation between the boundaries. This means that if a subunit is fully-
//   //enclosed the nearest point to it will always be its superunit. If a subunit
//   //is on the exterior, then the nearest points will (likely) be split between
//   //two superunits. The exception is if a subunit shares a boundary with its
//   //superunit. In this case the subunit may still have all its nearest points
//   //inside of one superunit. However, we can say that any points on a subunit
//   //sufficiently close to a superunit point indicate that the border is shared.
//   #pragma omp parallel for
//   for(unsigned int mpi=0;mpi<superunits.size();mpi++){
//     //Shrink the superunit by `shrink` units
//     const auto shrunk_superunit = BufferPath(superunits.at(mpi), -shrink);

//     //const auto shrunk_superunit = superunits.at(mpi).clipper_paths;

//     //Loop through all of the points of the superunit
//     #pragma omp critical
//     for(const auto &path: shrunk_superunit)
//     for(const auto &pt: path){
//       //Add this point to the vector of border points
//       borders.push_back(std::vector<double>({{pt.X,pt.Y}})); //TODO: Data type
//       owner.push_back(mpi);
//     }
//   }

//   //kdtree that will be used to provide rapid nearest neighbour lookups among
//   //the points which belong to the borders of the superunits
//   typedef KDTreeVectorOfVectorsAdaptor< my_vector_of_vectors_t, double >  my_kd_tree_t;
//   my_kd_tree_t border_idx(2 /*dim*/, borders, 10 /* max leaf */ );
//   border_idx.index->buildIndex();

//   #pragma omp parallel for
//   for(unsigned int subi=0;subi<subunits.size();subi++){
//     //Alias the current subunit
//     auto &sub = subunits[subi]; 

//     //Holds parents so we can insert them multiple times
//     std::set<int> parents; 
//     //Holds onto exterior so we don't have to repeatedly access map
//     bool exterior = false; 

//     //Loop through all of the points in the subunit
//     for(const auto &poly: sub.v)
//     for(const auto &ring: poly)
//     for(const auto &pt: ring){
//       //One of the boundary points of a subunit
//       const double qp[2] = {pt.x,pt.y};

//       const size_t num_results = 1; //Number of nearest neighbours to find
//       size_t ret_index;             //Index of the  point that's been found
//       double out_dist_sqr;          //Squared distance to the point
//       //Get the nearest neighbour
//       nanoflann::KNNResultSet<double> resultSet(num_results);
//       resultSet.init(&ret_index, &out_dist_sqr);
//       border_idx.index->findNeighbors(resultSet, &qp[0], nanoflann::SearchParams(10));

//       //If the squared distance from the subunit point to its nearest neighbour
//       //superunit point is less than the square distance cutoff, then we the
//       //subunit as being exterior
//       if(out_dist_sqr<border_dist_cutoff)
//         exterior = true;

//       parents.insert(owner.at(ret_index));

//       //While it seems like we could break early if we've found an exterior
//       //child, this isn't so: we still want to know what other parents the child
//       //might have
//     }

//     //Note that a subunit might be an exterior subunit but have only one parent
//     //because its border is coincident with a superunit
//     if(exterior){
//       sub.props["EXTCHILD"] = "T";
//     } else {
//       //However, if a subunit is not an exterior unit, then it can, by
//       //definition, have only one parent.
//       assert(sub.parents.size()==1);
//       sub.props["EXTCHILD"] = "F";
//     }

//     for(const auto &p: parents)
//       sub.parents.emplace_back(p, -1);

//     // #pragma omp critical
//     // for(const auto &p: parents)
//     //   superunits.at(p).children.emplace_back(subi, -1);
//   }
// }



// void CalcParentOverlap(GeoCollection &subunits, GeoCollection &superunits){
//   const double not_parent_thresh = 0.997;

//   for(unsigned int i=0;i<subunits.size();i++){
//     auto &sub = subunits.at(i);
//     //If the subunit is an interior unit (not an exterior unit) or has only one
//     //parent, it must lie entirely in that parent
//     if(sub.props["EXTCHILD"]=="F" || sub.parents.size()==1){
//       assert(sub.parents.size()==1);
//       sub.parents.at(0).second = 1;
//       continue;
//     }

//     const double area = areaExcludingHoles(sub);

//     //At this point, we suspect that some portion of the subunit lies outside of
//     //the superunit parent. But, remember, we shrunk the superunit a bit, so
//     //we'll need to be careful about attribution. Here, we intersect the
//     //subunits with the whole area of the unshrunk superunits. Afterwards, we'll
//     //check to see if we did any misattribution of parenthood.

//     //Loop over the parent units
//     for(auto &p: sub.parents){
//       double iarea = IntersectionArea(sub, superunits.at(p.first));
//       p.second     = iarea/area;
//     }

//     //Sort parents in order of decreasing percent area
//     std::sort(
//       sub.parents.begin(),
//       sub.parents.end(),
//       [](const MultiPolygon::parent_t &a, const MultiPolygon::parent_t &b){
//         return a.second>b.second;
//       }
//     );

//     //Running total of percent sums
//     double contribution_sum = 0;
//     //Pointer to the current superunit parent (starts with the parent within
//     //which most of the subunit lies)
//     auto p = sub.parents.begin();
//     //Loop through all the parents
//     for(;p!=sub.parents.end();p++){
//       //Keep track of how much of the subunit's percent area is accounted for
//       contribution_sum += p->second;
//       //If enough of the subunit has been allocated, then we break, assuming the
//       //rest can be chalked up to floating point errors
//       if(contribution_sum>=not_parent_thresh)
//         break;
//     }
//     //Eliminate the other parents
//     sub.parents.erase(p,sub.parents.end());
//     //Add the remainder to the final parent so that things add up closer to 1
//     sub.parents.back().second += 1-contribution_sum;
//   }

//   //Make sure that the parents all have the same information as the children so
//   //we can access it from either direction
//   for(unsigned int i=0;i<subunits.size();i++)
//   for(auto &p: subunits[i].parents)
//     superunits.at(p.first).children.emplace_back(i, p.second);
// }



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

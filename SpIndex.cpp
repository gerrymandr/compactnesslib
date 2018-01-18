#include "SpIndex.hpp"
#include <string>
#include <stdexcept>

namespace complib {

SpatialIndex::Region BoundingBoxToRegion( const BoundingBox &rect ){
  const double pt1[2] = { rect.xmin(), rect.ymin() };
  const double pt2[2] = { rect.xmax(), rect.ymax() };
  return SpatialIndex::Region( pt1, pt2, 2 );
}

class ObjVisitor : public SpatialIndex::IVisitor {
 private:
  std::vector<unsigned int>& mList;

 public:
  explicit ObjVisitor( std::vector<unsigned int> &list ) : mList( list ) {}

  void visitNode( const SpatialIndex::INode &n ) override {}

  void visitData( const SpatialIndex::IData &d ) override {
    mList.emplace_back( d.getIdentifier() );
}

  void visitData( std::vector<const SpatialIndex::IData *> &v ) override {}
};



class IteratorDataStream : public SpatialIndex::IDataStream {
 private:
  const idbb& mFi;
  idbb::const_iterator iter;
  bool started = false;
  SpatialIndex::RTree::Data *mNextData = nullptr;

  void readNextEntry() {
    if(iter==mFi.end()){
      mNextData = nullptr;
      return;
    }
    auto r = BoundingBoxToRegion(iter->second);
    mNextData = new SpatialIndex::RTree::Data(
      0,
      nullptr,
      r,
      iter->first
    );
    iter++;
  }

 public:
  explicit IteratorDataStream( const idbb& fi ) : mFi( fi ) {
    iter = mFi.begin();
    readNextEntry();
  }

  ~IteratorDataStream() {
    delete mNextData;
  }

  SpatialIndex::IData* getNext() override {
    SpatialIndex::RTree::Data *ret = mNextData;
    mNextData = nullptr;
    readNextEntry();
    return ret;
  }

  bool hasNext() override {
    return iter != mFi.end();
  }

  uint32_t size() override {
    throw std::runtime_error("IteratorDataStream::size() not available!");
  }

  void rewind() override { 
    throw std::runtime_error("IteratorDataStream::size() not available!");
  }
};


class SpIndexData {
 public:
  SpatialIndex::IStorageManager *mStorage = nullptr;
  SpatialIndex::ISpatialIndex *mRTree     = nullptr;

  SpIndexData() {
    initTree();
  }

  explicit SpIndexData( const idbb &fi){
    IteratorDataStream fids( fi );
    initTree( &fids );
  }

  ~SpIndexData() {
    delete mRTree;
    delete mStorage;
  }

  SpIndexData &operator=( const SpIndexData &rh ) = delete;

  void initTree( SpatialIndex::IDataStream *inputStream = nullptr ){
    // for now only memory manager
    mStorage = SpatialIndex::StorageManager::createNewMemoryStorageManager();

    // R-Tree parameters
    const double fillFactor                         = 0.7;
    const unsigned long indexCapacity               = 10;
    const unsigned long leafCapacity                = 10;
    const unsigned long dimension                   = 2;
    const SpatialIndex::RTree::RTreeVariant variant = SpatialIndex::RTree::RV_RSTAR;

    // create R-tree
    SpatialIndex::id_type indexId;

    if ( inputStream )
      mRTree = SpatialIndex::RTree::createAndBulkLoadNewRTree(
        SpatialIndex::RTree::BLM_STR, *inputStream, *mStorage, fillFactor, indexCapacity,
               leafCapacity, dimension, variant, indexId );
    else
      mRTree = SpatialIndex::RTree::createNewRTree( *mStorage, fillFactor, indexCapacity,
                                      leafCapacity, dimension, variant, indexId );
  }
};






SpIndex::SpIndex(){
  d.reset(new SpIndexData);
}

SpIndex::SpIndex( const idbb &fi){
  d.reset(new SpIndexData( fi ) );
}

void SpIndex::buildIndex(){
  d.reset(new SpIndexData(boxes_to_insert));
}

SpIndex& SpIndex::operator=( const SpIndex &other ){
  if ( this != &other )
    d = other.d;
  return *this;
}

void SpIndex::insert( unsigned int id, const BoundingBox &rect ){
  SpatialIndex::Region r( BoundingBoxToRegion( rect ) );

  try {
    d->mRTree->insertData( 0, nullptr, r, id );
  } catch ( Tools::Exception &e ) {
    throw std::runtime_error(std::string("libspatialindex tools exception: ")+e.what() );
  } catch ( const std::exception &e ) {
    throw std::runtime_error(std::string("std::exception in SpIndex: ")+e.what());
  } catch ( ... ) {
    throw std::runtime_error(std::string("Unknown libspatialindex exception caught!"));
}
}

void SpIndex::insertDeferred( const unsigned int id, const BoundingBox &bb ){
  boxes_to_insert.emplace_back(id,bb);
}

std::vector<unsigned int> SpIndex::query( const MultiPolygon &mp  ) const {
  return query(mp.bbox());
}

std::vector<unsigned int> SpIndex::query( const Point2D &pt       ) const {
  return query(BoundingBox(pt.x,pt.y,pt.x,pt.y));
}

std::vector<unsigned int> SpIndex::query( const BoundingBox &rect ) const {
  const SpatialIndex::Region r( BoundingBoxToRegion(rect) );

  //Object that will show us what results we've found
  std::vector<unsigned int> ret;
  ObjVisitor visitor(ret);

  #pragma omp critical
  d->mRTree->intersectsWithQuery( r, visitor );

  //Copy results    
  return ret;
}

/**
  Adds the bounding box of the multipolygon \p mp into SpIndex \sp associating
  it with an \p id. The bounding box is expanded by \p expandby units outward in 
  all directions.
*/
void AddToSpIndex(const MultiPolygon &mp, SpIndex &sp, const unsigned int id, const double expandby){
  auto bb = mp.bbox();
  bb.expand(expandby);
  sp.insertDeferred(id, bb);
}


}

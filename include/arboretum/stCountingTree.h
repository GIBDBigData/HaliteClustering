// -*-c++-*-

#ifndef __STCOUNTINGTREE_H
#define __STCOUNTINGTREE_H

#include "stCell.h"
#include "lmdbcxx/lmdb++.h"
#include <string>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "Normalization.h"
namespace Halite {

  template<typename D>
  class stCountingTree {
  public:
    stCountingTree(int H, const std::string& tmpdir, int DIM) {
      //empty tree
      this->H = H;
      sumOfPoints = 0;

      P.resize(DIM,0);

      lmdb_env=lmdb::env::create();
      
      //500 GiB
      lmdb_env.set_mapsize(500*1024*1024*1024L);
      lmdb_env.set_max_dbs(H-1);
      
      lmdb_env.open(tmpdir.c_str(), MDB_NOSYNC | MDB_WRITEMAP | MDB_NOLOCK);
      lmdb_txn = lmdb::txn::begin(lmdb_env);

       
      //create/open one dataset per tree level
      for (int i=0; i<H-1; i++) {
	//create dbName
	std::string dbName=str(boost::format("level_%1%.db") % i);
	//boost::filesystem::remove(dbName);
	

	//create and open dataset
	lmdb_levels.push_back(lmdb::dbi::open(lmdb_txn, dbName.c_str(),  MDB_CREATE));

      }
    }
    ~stCountingTree() {
      lmdb_txn.abort();
      
      /*      for (size_t i=0; i<H-1; i++) {
  	std::string dbName=str(boost::format("level_%1%.db") % i);
	boost::filesystem::remove(dbName);
	}*/

    }

    char insertPoint(const D *point, D *normPoint) {
      const size_t nPos = (P.size() + 7) / 8;

      //creates used arrays
      std::vector<D> min(P.size(), 0);
      std::vector<D> max(P.size(), 1);
      std::vector<unsigned char> fullId((H-1)*nPos, 0);

      std::vector<stCell> cellAndParents;
      for (size_t i = 0; i < H - 1; ++i)
	cellAndParents.push_back(stCell(P.size()));

      // normalizes the point
      char out = 0;
      normalization->normalize(point, normPoint);
      for (unsigned int i=0; i<P.size(); i++) {
	if (normPoint[i] < 0 || normPoint[i] > 1) {
	  out = 1; // invalid point
	}
      }

      if (!out) { // if the point is valid
	sumOfPoints++; // counts this point in the root
	// inserts this point in all tree levels
	deepInsert(0, min.data(), max.data(), normPoint, fullId, cellAndParents);
      }//end if

      // return the point "validity"
      return !out;
    }

    lmdb::dbi* getLMDB(size_t level) {
      return (level<H-1)?&lmdb_levels[level]:NULL;
    }

    lmdb::txn* getLMDBtxn() {
      return &lmdb_txn;
    }
 
    int getSumOfPoints() {
      return sumOfPoints;
    }
    int *getP() {
      return &P[0];
    }

    void setNormalization(std::shared_ptr<Normalization<D>> normalization) {
      this->normalization=normalization;
    }
    int findInNode(std::vector<stCell>& parents, stCell* cell, stCellId& id, int level) {
      const size_t nPos = (P.size() + 7) / 8;
      std::vector<unsigned char> fullId((level+1)*nPos, 0);

      //use id and parents to define fullId
      for (int l=0; l<level; l++) {
	parents[l].getId()->getIndex(&fullId[l*nPos]);
      }
      id.getIndex(&fullId[level*nPos]);

      //search the cell in current level
      lmdb::val value;
      bool lmdb_found = lmdb::dbi_get(lmdb_txn, lmdb_levels[level],
				      lmdb::val(fullId.data(), (level+1)*nPos),
				      value);
  
      if (lmdb_found) {
	*cell=stCell::deserialize(value.data<unsigned char>());
	cell->id = id;
      }

      return lmdb_found;
    }

    void commitCell(std::vector<stCell>& parents, stCell *cell, size_t level) {
      const size_t nPos = (P.size() + 7) / 8;
      std::vector<unsigned char> fullId((level+1)*nPos, 0);

      //use cell and parents to define fullId
      for (size_t l=0; l<level; l++) {
	parents[l].getId()->getIndex(&fullId[l*nPos]);
      }
      cell->getId()->getIndex(&fullId[level*nPos]);
    
      // insert/update the dataset
      unsigned char* const serialized = cell->serialize();

      lmdb::dbi_put(lmdb_txn, lmdb_levels[level],
		    lmdb::val(fullId.data(), (level+1)*nPos),
		    lmdb::val(serialized, stCell::size(P.size())));
      delete[] serialized;
    }
    void findParents(unsigned char *fullId, std::vector<stCell>& parentsVector, int level){
      const size_t nPos = (P.size() + 7) / 8;
      
      for (int l=0; l<level; l++) {
	//search the parent in level l
	lmdb::val key(fullId, (l+1)*nPos), value;
	
	bool lmdb_found=lmdb::dbi_get(lmdb_txn, lmdb_levels[l],
				     key, value);
	if(!lmdb_found) {
	  std::cerr<<"Error: NOT FOUND!\n";
	}
	parentsVector[l] = stCell::deserialize(value.data<unsigned char>());
	parentsVector[l].id.setIndex(fullId+(l*nPos));

      }
    }

  private:
    void deepInsert(size_t level, D* min, D *max, D *point, std::vector<unsigned char>& fullId, std::vector<stCell>& cellAndParents) {
      if (level < H-1) {
	// mounts the id of the cell that covers point in the current level
	// and stores in min/max this cell's lower and upper bounds
	D middle;
	stCellId cellId = stCellId(P.size());
	for (unsigned int i=0; i<P.size(); i++) {
	  middle = (max[i] - min[i]) / 2 + min[i];
	  if (point[i] > middle) { // if point[i] == middle, considers that the point
	    // belongs to the cell in the lower half regarding i
	    cellId.invertBit(i,P.size()); // turns the bit i to one
	    min[i] = middle;
	  } else {
	    // the bit i remains zero
	    max[i] = middle;
	  }
	}

	// updates the half space counts in level-1
	if (!level) { // updates in the root (level == 0)
	  for (unsigned int i=0; i<P.size(); i++) {
	    if (!cellId.getBitValue(i,P.size())) {
	      P[i]++; // counts the point, since it is in the root's lower half regarding i
	    }
	  }
	} else { // updates in the father (level > 0)
	  cellAndParents[level-1].insertPointPartial(&cellId, P.size());
	}

	// mounts the full id for the current cell
	int nPos = (P.size() - 1) / 8 + 1;
	cellId.getIndex(&fullId[level*nPos]);

	//searchs the cell in current level
	lmdb::val value;
	bool lmdb_found=lmdb::dbi_get(lmdb_txn, lmdb_levels[level],
				      lmdb::val(fullId.data(), (level+1)*nPos),
				      value);

	if (lmdb_found) {
	  cellAndParents[level]  = stCell::deserialize(value.data<unsigned char>());
	}

	//counts the new point
	cellAndParents[level].insertPoint();

	//inserts the point in the next level
	//and udates the partial counts for cellAndParents[level]
	deepInsert(level+1, min, max, point, fullId, cellAndParents);


	// insert/update the dataset
	unsigned char* serialized = cellAndParents[level].serialize();

	lmdb::dbi_put(lmdb_txn, lmdb_levels[level],
		      lmdb::val(fullId.data(), (level+1)*nPos),
		      lmdb::val(serialized, stCell::size(P.size())));

	delete[] serialized;
      }
    }

    std::vector<int> P;
    size_t sumOfPoints;
    size_t H;

    std::vector<lmdb::dbi> lmdb_levels;
    
    lmdb::env lmdb_env=NULL;
    lmdb::txn lmdb_txn=NULL;
    std::shared_ptr<Halite::Normalization<D>> normalization;
  };

}
#endif //__STCOUNTINGTREE_H

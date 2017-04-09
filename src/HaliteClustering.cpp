/**********************************************************************
 * GBDI Arboretum - Copyright (c) 2002-2009 GBDI-ICMC-USP
 *
 *                           Homepage: http://gbdi.icmc.usp.br/arboretum
 **********************************************************************/
/* ====================================================================
 * The GBDI-ICMC-USP Software License Version 1.0
 *
 * Copyright (c) 2009 Grupo de Bases de Dados e Imagens, Instituto de
 * Ciencias Matemáticas e de Computação, University of São Paulo -
 * Brazil (the Databases and Image Group - Intitute of Matematical and
 * Computer Sciences).  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *    if any, must include the following acknowledgment:
 *       "This product includes software developed by Grupo de Bases
 *        de Dados e Imagens, Instituto de Ciências Matemáticas e de
 *        Computação, University of São Paulo - Brazil (the Databases
 *        and Image Group - Intitute of Matematical and Computer
 *        Sciences)"
 *
 *    Alternately, this acknowledgment may appear in the software itself,
 *    if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names of the research group, institute, university, authors
 *    and collaborators must not be used to endorse or promote products
 *    derived from this software without prior written permission.
 *
 * 5. The names of products derived from this software may not contain
 *    the name of research group, institute or university, without prior
 *    written permission of the authors of this software.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OF THIS SOFTWARE OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * ====================================================================
 *                                            http://gbdi.icmc.usp.br/
 */
/**
 * @file
 * This file implements the class HaliteClustering.
 *
 * @version 1.0
 * @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
 * @author Agma Juci Machado Traina (agma@icmc.usp.br)
 * @author Christos Faloutsos (christos@cs.cmu.edu)
 * @author Caetano Traina Jr (caetano@icmc.usp.br)
 */
// Copyright (c) 2002-2009 GBDI-ICMC-USP

//------------------------------------------------------------------------------
// class HaliteClustering
//------------------------------------------------------------------------------

#include "HaliteClustering.h"
#include <vector>
#include <memory>
#include <boost/pending/disjoint_sets.hpp>
#include "BetaCluster.h"


namespace Halite {

  /* size_t HaliteClustering::getCenter(size_t level) {
    size_t center;

    //if (betaClusterCenter.getId()->getBitValue(0,DIM)) {
    //center = total - betaClusterCenterParents[level-1]->getP(0);
    //} else {
    center = betaClusterCenterParents[level-1].getP(0);
    //}
    return center;
    }*/
  template <typename D>
  HaliteClustering<D>::HaliteClustering (PointSource<D>& data, bool hardClustering, const std::string& tmpdir, NormalizationMode normalizationMode, D pThreshold, int H) {
    
    // stores DIM, H, hardClustering and initialLevel
    this->DIM=data.dimension();
    this->H = H;
    this->hardClustering = hardClustering;

    //Start at the 2nd level, since you need to look at your
    //parents for the statistical significance test
    this->initialLevel = 1;

    //"Laplacian" filter, using only direct neighbors
    this->centralConvolutionValue = 2 * DIM;
    this->neighbourhoodConvolutionValue = -1;

    // stores the pThreshold
    this->pThreshold = pThreshold;

    classifier=std::make_shared<Classifier<D> >();
    classifier->hardClustering=hardClustering;

    // builds the counting tree and inserts objects on it
    calcTree = new stCountingTree<D>(H, tmpdir, DIM);
        
    timeNormalization = clock(); //start normalization time
    readData(data, normalizationMode);
    timeNormalization = (clock()-timeNormalization); //total time spent in the normalization
    
  }//end HaliteClustering::HaliteClustering

  //---------------------------------------------------------------------------
  template <typename D>
  HaliteClustering<D>::~HaliteClustering() {

    // disposes the used structures
    delete calcTree;

  }//end HaliteClustering::~HaliteClustering

  //---------------------------------------------------------------------------
  template <typename D>
  void HaliteClustering<D>::findCorrelationClusters() {

     /**
     * Vector used when discovering relevant attributes.
     */
    std::vector<D> attributesRelevance(DIM);

    /**
     * Vectors to describe the positions of the beta-cluster center in the data space.
     */
    std::vector<D> minBetaClusterCenter(DIM,0.0);
    std::vector<D> maxBetaClusterCenter(DIM,0.0);

    // Vector used to indicate which direct neighbours of a central cell also belong to the beta-cluster.
    std::vector<char> neighbourhood(DIM);
    
    // Pointers to the central cell of a found beta-cluster and to a neighbour cell.
    stCell betaClusterCenter(DIM);
    // Vectors with pointers to the parents of a beta-cluster center and of a neighbour.
    std::vector<stCell> betaClusterCenterParents(H,stCell(DIM));

    std::vector<size_t> old_center(DIM, 0);
    std::vector<size_t> old_critical(DIM, 0);

   
    // defines when a new cluster is found
    size_t center;
    int ok, total;
    do { // looks for a cluster in each loop
      ok=0; // no new cluster was found
      // defines the initial grid level to analyze
      int level=initialLevel;
      do { // analyzes each level until a new cluster is found
	// apply the convolution matrix to each grid cell in the current level
	// finds the cell with the bigest convolution value
	if (walkThroughConvolution(level,betaClusterCenter, betaClusterCenterParents)) {
	  betaClusterCenter.useCell(); // visited cell
	  calcTree->commitCell(betaClusterCenterParents, &betaClusterCenter, level); // commit changes in the tree
	  if (level) {
	    // pointer to a neighbour of the father
	    stCell fatherNeighbour;
	    std::vector<stCell> cellP(level-1, stCell(DIM));

	    for (size_t i = 0; i < DIM; i++) {
	      // initiates total with the number of points in the father
	      total=betaClusterCenterParents[level-1].getSumOfPoints();

	      // discovers the number of points in the center
	      if (betaClusterCenter.getId()->getBitValue(i,DIM)) {
		center = total - betaClusterCenterParents[level-1].getP(i);
	      } else {
		center = betaClusterCenterParents[level-1].getP(i);
	      }//end if

	      // looks for the points in the direct neighbours of the father
	      if (internalNeighbour(i,betaClusterCenterParents[level-1],&fatherNeighbour,betaClusterCenterParents,level-1)) {
		total += fatherNeighbour.getSumOfPoints();
	      }//end if

	      std::vector<stCell> cellP = betaClusterCenterParents;

	      // looks for the external neighbour
	      if (externalNeighbour(i,betaClusterCenterParents[level-1],&fatherNeighbour,betaClusterCenterParents,cellP,level-1)) {
		total += fatherNeighbour.getSumOfPoints();
	      }//end if

	      // percentual of points in the center related to the average
	      attributesRelevance[i] = (100*center)/((D)total/6);
	      // right critical value for the statistical test
	      int criticalValue = GetCriticalValueBinomialRight2(total, (D)1/6, pThreshold);
	      size_t critVal=(size_t)criticalValue;
	      if (center > critVal && old_center[i] != center && old_critical[i] != critVal) {
		ok=1; // new cluster found
	      }//end if

	      old_center[i] = center;
	      old_critical[i] = criticalValue;
	    }//end for
	  } else { // analyzes each dimension based on the points distribution of the entire database
	    // initiate the total of points
	    total=calcTree->getSumOfPoints();
	    for (size_t i = 0; i < DIM; i++) {
	      // discovers the number of points in the center
	      if (betaClusterCenter.getId()->getBitValue(i,DIM)) {
		center = total - calcTree->getP()[i];
	      } else {
		center = calcTree->getP()[i];
	      }//end if
	      // percentual of points in the center related to the average
	      attributesRelevance[i] = (100*center)/((D)total/2);
	      // right critical value for the statistical test
	      int criticalValue = GetCriticalValueBinomialRight2(total, (D)1/2, pThreshold);
	      if (center > (size_t) criticalValue) {
		ok=1; // new cluster found
	      }//end if
	    }//end for
	  }//end if
	}//end if

	if (!ok) {
	  level++; // next level to be analyzed
	}//end if
      } while (!ok && level < H);//end do while

      if (!ok) {
	continue;
      }

      // discovers the cThreshold based on the minimum description length method
      const D cThreshold = calcCThreshold(attributesRelevance);
      
      // new cluster found
      classifier->betaClusters.push_back(BetaCluster<D>(level, DIM));
      printf("a beta-cluster was found at the Counting-tree level %d (%zu).\n",level,numBetaClusters()); // prints the level in which a new beta-cluster was found

      BetaCluster<D>& newCluster=classifier->betaClusters.back();

      // important dimensions
      for (size_t i = 0; i < DIM; i++) {
	newCluster.relevantDimension[i] = (attributesRelevance[i] >= cThreshold);
      }//end for
      
      // analyzes neighbours in important dimensions to verify which of them also belong to the found cluster
      for (size_t i = 0; i < DIM; i++) {
	neighbourhood[i] = 'N'; // no direct neighbour belongs to the cluster
      }//end for
      // center's position in the data space
      cellPosition(betaClusterCenter,betaClusterCenterParents,minBetaClusterCenter,maxBetaClusterCenter,level);
      // for each important dimension, analyzes internal and external neighbours to decide if they also
      // belong to the cluster
      for (size_t i = 0; i < DIM; i++) {
	if (newCluster.relevantDimension[i]) {
	  // looks for the internal neighbour
	  stCell neighbour;
	  if (internalNeighbour(i,betaClusterCenter,&neighbour,betaClusterCenterParents,level)) { // internal neighbour in important dimension always belongs to the cluster
	    //Position of a neighbour of the beta-cluster center in the data space, regarding a dimension e_j.
	    D minNeighbour, maxNeighbour;
	    cellPositionDimensionE_j(neighbour,betaClusterCenterParents,&minNeighbour,&maxNeighbour,level,i);
	    if (maxBetaClusterCenter[i] > maxNeighbour) {
	      neighbourhood[i]='I'; // inferior neighbour in i belongs to the cluster
	    } else {
	      neighbourhood[i]='S'; // superior neighbour in i belongs to the cluster
	    }//end if
	  }//end if

	  // looks for the external neighbour
	  std::vector<stCell> neighbourParents=betaClusterCenterParents;
	  
	  if (externalNeighbour(i,betaClusterCenter,&neighbour,betaClusterCenterParents,neighbourParents,level)) { // analyzes external neighbour to decide if it belongs to the cluster
	    if (neighbourhood[i] == 'N') {
	      //Position of a neighbour of the beta-cluster center in the data space, regarding a dimension e_j.
	      D minNeighbour,maxNeighbour;
	      cellPositionDimensionE_j(neighbour,neighbourParents,&minNeighbour,&maxNeighbour,level,i);
	      if (maxBetaClusterCenter[i] > maxNeighbour) {
		neighbourhood[i]='I'; // inferior neighbour in i belongs to the cluster
	      } else {
		neighbourhood[i]='S'; // superior neighbour in i belongs to the cluster
	      }//end if
	    } else {
	      neighbourhood[i]='B'; // both inferior and superior neighbours in i belong to the cluster
	    }//end if
	  }//end if
	}//end if
      }//end for

      // stores the description of the found cluster
      const D length = maxBetaClusterCenter[0] - minBetaClusterCenter[0];
      for (size_t i = 0; i < DIM; i++) {
	if (newCluster.relevantDimension[i]) { // dimension important to the cluster
	  // analyzes if the neighbours in i also belong to the cluster
	  switch (neighbourhood[i]) {
	  case 'B': // both inferior and superior neighbours in i belong to the cluster
	    minBetaClusterCenter[i]-=length;
	    maxBetaClusterCenter[i]+=length;
	    break;
	  case 'S': // superior neighbour in i belongs to the cluster
	    maxBetaClusterCenter[i]+=length;
	    break;
	  case 'I': // inferior neighbour in i belongs to the cluster
	    minBetaClusterCenter[i]-=length;
	  }//end switch
	  // new cluster description - relevant dimension
	  newCluster.min[i] = minBetaClusterCenter[i];
	  newCluster.max[i] = maxBetaClusterCenter[i];
	} else {
	  // new cluster description - irrelevant dimension
	  newCluster.min[i] = 0;
	  newCluster.max[i] = 1;
	}//end if
      }//end for
      // stops when no new cluster is found
    } while (ok);//end do while

    printf("\n%lu beta-clusters were found.\n",numBetaClusters()); // prints the number of beta-clusters found
    mergeBetaClusters(); // merges clusters that share some database space
    printf("\n%d correlation clusters left after the merging phase.\n",numCorrelationClusters()); // prints the number of correlation clusters found

  }//end HaliteClustering::findCorrelationClusters

  //---------------------------------------------------------------------------
  template<typename D>
  D HaliteClustering<D>::calcCThreshold(const std::vector<D>& attributesRelevance) {
    std::vector<D> sortedRelevance(DIM);

    for (size_t i = 0; i < DIM; i++) {
      sortedRelevance[i]=attributesRelevance[i];
    }//end for
    std::sort(sortedRelevance.begin(), sortedRelevance.end());

    D cThreshold = sortedRelevance[minimumDescriptionLength(sortedRelevance)];

    return cThreshold;

  }//end HaliteClustering::calcCThreshold

  //For positive numbers only
  template<typename D>
  inline D log2ceil(D d) {
    D cl = ceil(d);
    return cl ? log2(cl) : 0;
  }
  //---------------------------------------------------------------------------
  template<typename D>
  int HaliteClustering<D>::minimumDescriptionLength(const std::vector<D>& sortedRelevance) {
    
    int cutPoint=-1;
    D preAverage, postAverage, descriptionLength, minimumDescriptionLength;
    for (size_t i = 0; i < DIM; i++) {
      descriptionLength=0;
      // calculates the average of both sets
      preAverage=0;
      for (size_t j = 0;j<i;j++) {
	preAverage+=sortedRelevance[j];
      }//end for
      if (i) {
	preAverage/=i;
	descriptionLength += log2ceil(preAverage);	
      }//end if
      postAverage=0;
      for (size_t j = i;j<DIM;j++) {
	postAverage+=sortedRelevance[j];
      }//end for
      if (DIM-i) {
	postAverage/=(DIM-i);
	descriptionLength += log2ceil(postAverage);
      }//end if
      // calculates the description length
      for (size_t j = 0;j<i;j++) {
	descriptionLength += log2ceil(fabs(preAverage-sortedRelevance[j]));
      }//end for
      for (size_t j = i;j<DIM;j++) {
	descriptionLength += log2ceil(fabs(postAverage-sortedRelevance[j]));
      }//end for
      // verify if this is the best cut point
      if (cutPoint==-1 || descriptionLength < minimumDescriptionLength) {
	cutPoint=i;
	minimumDescriptionLength = descriptionLength;
      }//end if
    }//end for
    return cutPoint;

  }//end HaliteClustering::minimumDescriptionLength

  //---------------------------------------------------------------------------
  template<typename D>
  int HaliteClustering<D>::walkThroughConvolution(int level, stCell& betaClusterCenter, std::vector<stCell>& betaClusterCenterParents) {
    //try to get the db in the current level
    lmdb::dbi* lmdbi = calcTree->getLMDB(level);
    lmdb::txn* lmdbtxn = calcTree->getLMDBtxn();

    if (!lmdbi)
      return 0;

    // pointers to the parents of a cell
    std::vector<stCell> parentsVector(level, stCell(DIM));

    //prepare the fullId array
    int nPos = (int) ceil((D)DIM/8);

    std::vector<unsigned char> ccFullId((level+1)*nPos,0);

    //prepare the cell and the Dbts to receive
    //key/data pairs from the dataset
    lmdb::val lmdb_key, lmdb_val;
    
    // Get a cursor
    lmdb::cursor lmdb_cursor = lmdb::cursor::open(*lmdbtxn, *lmdbi);
      

    //prepare to walk through the level and find the
    //cell with the biggest convoluted value
    int biggestConvolutionValue = -MAXINT;
    int newConvolutionValue;
    char clusterFoundBefore;
    std::vector<D> maxCell(DIM, 0.0);
    std::vector<D> minCell(DIM, 0.0);

    // iterate over the database, retrieving each record in turn
    while (lmdb_cursor.get(lmdb_key, lmdb_val, MDB_NEXT)) {
     stCell cell = stCell::deserialize(lmdb_val.data<unsigned char>());

      
     //Does not analyze cells analyzed before
       
     if(cell.getUsedCell()) {
       continue;
     }
     // Does not analyze cells that can't be the biggest convolution center.
     // It speeds up the algorithm when neighbourhoodConvolutionValue <= 0

     if(neighbourhoodConvolutionValue<=0) {
       if(cell.getSumOfPoints()*centralConvolutionValue < biggestConvolutionValue
	  || (cell.getSumOfPoints()*centralConvolutionValue == biggestConvolutionValue
	      &&  memcmp(lmdb_key.data(), ccFullId.data(), lmdb_key.size())>=0)) {
	 continue;
       }
     }

      //set id for cell
      cell.id.setIndex(lmdb_key.data<unsigned char>()+ level*nPos);
      //finds the parents of cell
      calcTree->findParents(lmdb_key.data<unsigned char>(),parentsVector,level);

      // discovers the position of cell in the data space
      cellPosition(cell,parentsVector,minCell,maxCell,level);

      // verifies if this cell belongs to a cluster found before
      clusterFoundBefore=0;
      for (BetaCluster<D>& bCluster: classifier->betaClusters) {
	clusterFoundBefore = 1;
	for (size_t j = 0; j<DIM; j++) {
	  // Does not cut off cells in a level upper than the level where a cluster was found
	  if (!(maxCell[j] <= bCluster.max[j] && minCell[j] >= bCluster.min[j])) {
	    clusterFoundBefore = 0;
	    break;
	  }//end if
	}//end for
	if(clusterFoundBefore)
	  break;
      }//end for

      if (clusterFoundBefore) { // the cell doesn't belong to any found cluster
	// applies the convolution matrix to cell
	continue;
      }

      if (neighbourhoodConvolutionValue) {
	newConvolutionValue=applyConvolution(cell,parentsVector,level);
      } else {
	newConvolutionValue=centralConvolutionValue*(cell.getSumOfPoints()); // when the neighbourhood weight is zero
      }//end if

      if ( (newConvolutionValue > biggestConvolutionValue) ||
	   ((newConvolutionValue == biggestConvolutionValue) && (memcmp(lmdb_key.data(), ccFullId.data(), lmdb_key.size()) < 0))) {

	// until now, cell is the biggest convolution value, thus, set the new biggest value and copy
	// cell and its parents to betaClusterCenter and its parents
	biggestConvolutionValue = newConvolutionValue;
	std::copy(lmdb_key.data(), lmdb_key.data()+lmdb_key.size(), ccFullId.begin());
	
	betaClusterCenter=cell;
	for (int j = 0;j<level;j++) {
	  betaClusterCenterParents[j]=parentsVector[j];
	}//end for
      }//end if
    }//end while
 

    //closes the cursor
    lmdb_cursor.close();

    return 1; //Success
  }//end HaliteClustering::walkThroughConvolution

  //---------------------------------------------------------------------------
  template <typename D>
  int HaliteClustering<D>::applyConvolution(stCell& cell, std::vector<stCell>& cellParents, int level) {

    stCell neighbour;
 

    int newValue = cell.getSumOfPoints()*centralConvolutionValue;
    // looks for the neighbours
    for (size_t k=0;k<DIM;k++) {
      if (internalNeighbour(k,cell,&neighbour,cellParents,level)  ) {
	newValue+=(neighbour.getSumOfPoints()*neighbourhoodConvolutionValue);
      }//end if

      //copy parents
      std::vector<stCell> neighbourParents = cellParents;

      if (externalNeighbour(k,cell,&neighbour,cellParents,neighbourParents,level)) {
	newValue+=(neighbour.getSumOfPoints()*neighbourhoodConvolutionValue);
      }//end if
    }//end for
  
    // return the cell value after applying the convolution matrix
    return newValue;

  }//end HaliteClustering::applyConvolution

 
  //---------------------------------------------------------------------------
  template<typename D>
  void HaliteClustering<D>::cellPosition(stCell& cell, std::vector<stCell>& cellParents,
				      std::vector<D>& min, std::vector<D>& max, int level) {
      if (level) {
	cellPosition(cellParents[level-1],cellParents,min,max,level-1);
	for (size_t i = 0; i < DIM; i++) {
	  if (cell.id.getBitValue(i,DIM)) { // bit in the position i is 1
	    min[i] += ((max[i]-min[i])/2);
	  } else { // bit in the position i is 0
	    max[i] -= ((max[i]-min[i])/2);
	  }//end if
	}//end for
      } else { // level zero
	for (size_t i = 0; i < DIM; i++) {
	  if (cell.id.getBitValue(i,DIM)) { // bit in the position i is 1
	    min[i] = 0.5;
	    max[i] = 1;
	  } else { // bit in the position i is 0
	    min[i] = 0;
	    max[i] = 0.5;
	  }//end if
	}//end for
      }//end if
 
  }//end HaliteClustering::cellPosition

  //---------------------------------------------------------------------------
  template<typename D>
  void HaliteClustering<D>::cellPositionDimensionE_j(stCell& cell, std::vector<stCell>& cellParents,
						  D *min, D *max, int level, int j) {

      if (level) {
	cellPositionDimensionE_j(cellParents[level-1],cellParents,min,max,level-1,j);
	if (cell.id.getBitValue(j,DIM)) { // bit in the position j is 1
	  *min += ((*max-*min)/2);
	} else { // bit in the position j is 0
	  *max -= ((*max-*min)/2);
	}//end if
      } else { // level zero
	if (cell.id.getBitValue(j,DIM)) { // bit in the position j is 1
	  *min = 0.5;
	  *max = 1;
	} else { // bit in the position j is 0
	  *min = 0;
	  *max = 0.5;
	}//end if
      }//end if
 
  }//end HaliteClustering::cellPositionDimensionE_j

  //---------------------------------------------------------------------------
  template<typename D>
  int HaliteClustering<D>::externalNeighbour(int dimIndex, stCell& cell, stCell* neighbour,
					  std::vector<stCell>& cellParents, std::vector<stCell>& neighbourParents, int level) {
    if (level) {
      int found;
      stCell& father = cellParents[level-1];
      if (cell.getId()->getBitValue(dimIndex,DIM) ^ father.getId()->getBitValue(dimIndex,DIM)) { // XOR - different bit values -> starts going down in the tree
	//finds the internal neighbour of the father
	found = internalNeighbour(dimIndex, father, &neighbourParents[level-1], cellParents, level-1);
      } else {  // equal bit values -> continues going up in the tree
	// recursively, finds the external neighbour of the father
	found = externalNeighbour(dimIndex,father,&neighbourParents[level-1],cellParents,neighbourParents,level-1);
      }//end if
      if (found) { // father's neighbour was found
	// find the external neighbour of cell in dimension i
	return internalNeighbour(dimIndex, cell, neighbour, neighbourParents, level);
      }//end if
      return 0; // there is no father's neighbour
    }//end if
    return 0; // a cell in level zero never has an external neighbour

  }//end HaliteClustering::externalNeighbour

  //---------------------------------------------------------------------------
  template<typename D>
  int HaliteClustering<D>::internalNeighbour(int dimIndex, stCell& cell, stCell* neighbour,
					  std::vector<stCell>& cellParents, int level) {
    // creates the id that the neighbour should have
    stCellId neighboursId = cell.id;
    neighboursId.invertBit(dimIndex, DIM);

    int found = calcTree->findInNode(cellParents, neighbour, neighboursId, level);

    return found;

  }//end HaliteClustering::internalNeighbour

  //---------------------------------------------------------------------------
  template<typename D>
  void HaliteClustering<D>::readData(PointSource<D>& data, NormalizationMode mode) {
    normalization=std::make_shared<Normalization<D>>(data,mode);
    calcTree->setNormalization(normalization);
    classifier->normalization=normalization;

    std::vector<D> tmp(DIM);
    //We also calculate the number of data points here
    this->SIZE=0;
    //process each point
    for (data.restartIteration(); data.hasNext();) {
      const D* onePoint= data.readPoint();
      calcTree->insertPoint(onePoint,tmp.data()); //add to the grid structure
      this->SIZE++;
    }//end for


  }//end HaliteClustering::FastDistExponent

  template<typename D>
  void HaliteClustering<D>::mergeBetaClusters() {
    std::vector<size_t> rank(numBetaClusters());
    std::vector<size_t> parent(numBetaClusters());
    boost::disjoint_sets<size_t*,size_t* > ds(&rank[0], &parent[0]);

    for(size_t i=0; i<numBetaClusters(); i++) {
      ds.make_set(i);
    }
    for(size_t i=0; i<numBetaClusters(); i++) {
      BetaCluster<D>& iCl=classifier->betaClusters[i];
      for(size_t j=i+1; j<numBetaClusters(); j++) {
	BetaCluster<D>& jCl=classifier->betaClusters[j];
	if(shouldMerge(iCl,jCl)) {
	  ds.union_set(i,j);
	}
      }
    }
    int numCorrelationClusters=0;
    std::vector<int> relabeling(numBetaClusters(),-1);
    for(size_t i=0; i<numBetaClusters(); i++) {
      size_t rep = ds.find_set(i);
      if(relabeling[rep] == -1) {
	relabeling[rep] = numCorrelationClusters++;
      }
      classifier->betaClusters[i].correlationCluster = relabeling[rep];
    }
    
    correlationClusters.clear();
    correlationClusters.resize(numCorrelationClusters, CorrelationCluster(DIM));

    for (size_t i=0; i<numBetaClusters(); i++) {
      BetaCluster<D>& bCluster=classifier->betaClusters[i];
      for (size_t j = 0; j<DIM; j++) {
	CorrelationCluster& correlationCluster = correlationClusters[bCluster.correlationCluster];

	correlationCluster.relevantDimension[j] |=  bCluster.relevantDimension[j]; 
      }
    }

  }//end HaliteClustering::mergeBetaClusters

  template<typename D>
  int HaliteClustering<D>::shouldMerge(BetaCluster<D>& iCl,BetaCluster<D>& jCl) {
 
    // discovers if beta-cluster i shares database space with beta-cluster j
    int shareSpace=1;
    for(size_t k=0; shareSpace && k<DIM; k++) {
      if (!(iCl.max[k] > jCl.min[k] && iCl.min[k] < jCl.max[k])) {
	shareSpace=0; // beta-clusters i and j do not share database space
      }//end if
    }//end for

    if (shareSpace) {
      if (!hardClustering) {
	//does a PCA based analysis
	if (iCl.cost==-1) {
	  iCl.cost=cost(&iCl,NULL);
	}
	if (jCl.cost==-1) {
	  jCl.cost = cost(&jCl,NULL);
	}
	return (D)(iCl.cost+jCl.cost)>=cost(&iCl,&jCl); //merges if the merged cluster compacts best
      }
      return 1; //merge
    }
    return 0; //not merge
  }
  template<typename D>
  int HaliteClustering<D>::cost(BetaCluster<D>* iCl, BetaCluster<D>* jCl) {
    int cost=0;

    //prepare the input for PCA
    cv::Mat clusterMat = inputPCA(iCl, jCl);
    
    size_t clusterSize=clusterMat.rows;
    //applies PCA in the cluster
    if (clusterSize < DIM) {
      return 0; // not possible to apply PCA
    }
    cv::PCA princomp(clusterMat, // pass the data
		     cv::Mat(), // we do not have a pre-computed mean vector, so let the PCA engine to compute it
		     CV_PCA_DATA_AS_ROW, // indicate that the vectors are stored as matrix rows
		     (int)DIM // specify, how many principal components to retain
		     );
    clusterMat = princomp.project(clusterMat);

    //finds the minimum and maximum values of cluster in each PCA axis
    std::vector<D> min(DIM);
    std::vector<D> max(DIM);

    for (size_t d=0; d<DIM; d++) { // sets the values to the minimum/maximum possible here
      min[d] = std::numeric_limits<D>::max();
      max[d] = std::numeric_limits<D>::lowest();
    }// end for
    for (size_t p=0; p<clusterSize; p++) {
      for (size_t d=0; d<DIM; d++) {
	if (clusterMat.at<D>(p,d) > max[d]) {
	  max[d] = clusterMat.at<D>(p,d);
	}
	if (clusterMat.at<D>(p,d) < min[d]) {
	  min[d] = clusterMat.at<D>(p,d);
	}
      }
    }

    //cost for the points
    for (size_t p=0; p<clusterSize; p++) {
      for (size_t d=0; d<DIM; d++) {
	cost += indCost(clusterMat.at<D>(p,d)-(((max[d]-min[d])/2) + min[d])); //distance to the center in each dimension
      }
    }

    //cost for eigenvectors, min and max
    for (size_t d=0; d<DIM; d++) {
      cost += indCost(min[d]) + indCost(max[d]);
      for (size_t k=0; k<DIM; k++) {
	cost += indCost(princomp.eigenvectors.at<D>(d,k));
      }
    }

    //disposes the used memory
    clusterMat.release();

    return cost;
  }

  template<typename D>
  int HaliteClustering<D>::indCost(D n) {
    n = ceil(fabs(n*1000000)); //ignores the sign, since + and - cost the same, and
    //considers the cost of the smallest integer bigger than n
    if (n <= 1) {
      return 1; // zero and one cost one
    }
    return (int) ceil(log2(n)); // cost of n, when n > 1
  }
  template<typename D>
  cv::Mat HaliteClustering<D>::inputPCA(const BetaCluster<D>* iCl, const BetaCluster<D>* jCl) {
    //prepare the input for PCA
    std::vector<std::vector<D>> cluster;

    size_t iLevel=0, jLevel=0;
    if(iCl!=NULL) iLevel = iCl->level;
    if(jCl!=NULL) jLevel = jCl->level;
    size_t level=std::max(iLevel, jLevel);
  
    lmdb::dbi* lmdbi = calcTree->getLMDB(level);
    lmdb::txn* lmdbtxn = calcTree->getLMDBtxn();
    

    // pointers to the parents of a cell
    std::vector<stCell> parentsVector(level, stCell(DIM));

    //prepare the fullId array
    size_t nPos = (DIM + 7) / 8;

    lmdb::val lmdb_key, lmdb_val;
    
    // Get a cursor
    lmdb::cursor lmdb_cursor = lmdb::cursor::open(*lmdbtxn, *lmdbi);
   
    //prepare to walk through the level
    bool belongsTo;
 
    std::vector<D> maxCell(DIM, 0.0);
    std::vector<D> minCell(DIM, 0.0);

    // iterate over the database, retrieving each record in turn
    while (lmdb_cursor.get(lmdb_key, lmdb_val, MDB_NEXT)) {
      stCell cell = stCell::deserialize(lmdb_val.data<unsigned char>());
      
      //set id for cell
      cell.id.setIndex(lmdb_key.data<unsigned char>() + level*nPos);
      
      //finds the parents of cell
      calcTree->findParents(lmdb_key.data<unsigned char>(),parentsVector,level);
      // discovers the position of cell in the data space
      cellPosition(cell,parentsVector,minCell,maxCell,level);

      belongsTo = false;

      std::vector<D> cellCenter(DIM);
      for(size_t dim=0; dim<DIM; dim++) {
	cellCenter[dim]=(maxCell[dim]-minCell[dim])/2 + minCell[dim];
      }

      if(iCl != NULL && iCl->contains(cellCenter)) {
	belongsTo=true;
      } else if(jCl != NULL && jCl->contains(cellCenter)) {
	belongsTo=true;  
      }
     
      if (belongsTo) { // this cell belongs to the PCA input
	for (int p=0; p < cell.getSumOfPoints(); p++) {
	  cluster.emplace_back(DIM);
	  normalization->denormalize(cellCenter.begin(), cluster.back().begin());
	}
      }

    }
    
    //closes the cursor
    lmdb_cursor.close();

    //copy cluster to clusterMat (OpenCV format)
    cv::Mat clusterMat(cluster.size(), DIM, cv::DataType<D>::type);
    for (size_t p=0; p<cluster.size(); p++) {
      for (size_t d=0; d<DIM; d++) {
	clusterMat.at<D>(p,d) = cluster[p][d];
      }
    }

    return clusterMat;

  }

  template class HaliteClustering<float>;
  template class HaliteClustering<double>;
}



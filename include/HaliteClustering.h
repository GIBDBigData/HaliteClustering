// -*-c++-*-
/**********************************************************************
 * GBDI Arboretum - Copyright (c) 2002-2009 GBDI-ICMC-USP
 *
 *                           Homepage: http://gbdi.icmc.usp.br/arboretum
 **********************************************************************/
/* ====================================================================
 * The GBDI-ICMC-USP Software License Version 1.0
 *
 * Copyright (c) 2009 Grupo de Bases de Dados e Imagens, Instituto de
 * Ciências Matemáticas e de Computação, University of São Paulo -
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
 * This file defines the class HaliteClustering.
 *
 * @version 1.0
 * @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
 * @author Agma Juci Machado Traina (agma@icmc.usp.br)
 * @author Christos Faloutsos (christos@cs.cmu.edu)
 * @author Caetano Traina Jr (caetano@icmc.usp.br)
 *
 */
// Copyright (c) 2002-2009 GBDI-ICMC-USP

#ifndef __HALITECLUSTERING_H
#define __HALITECLUSTERING_H

#include <cv.h> //OpenCV
#include "arboretum/stCountingTree.h"
#include "PointSource.h"

#include "Utile.h"
#include "Classifier.h"
#include "CorrelationCluster.h"

#include <math.h>
#include <stdio.h>
#include <time.h>

namespace Halite {

  //----------------------------------------------------------------------------
  // class HaliteClustering
  //----------------------------------------------------------------------------
  /**
   * This class is used to find clusters in subspaces of the original data space.
   *
   * @version 1.0
   * @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
   * @author Agma Juci Machado Traina (agma@icmc.usp.br)
   * @author Christos Faloutsos (christos@cs.cmu.edu)
   * @author Caetano Traina Jr (caetano@icmc.usp.br)
   */
  //---------------------------------------------------------------------------
  template <typename D>
  class HaliteClustering {

  public:

    /**
     * Creates the needed structures to find correlation clusters.
     *
     * @param PointSource Source of the database objects.
     * @param pThreshold Threshold used to spot a beta-cluster, based on the binomial probability.
     * @param H The number of grid levels tobuild the counting tree.
     * @param hardClustering Choose between hard (1) and soft (0) clustering.
     *
     */
    HaliteClustering(PointSource<D>& data, bool hardClustering,
		     const std::string& tmpdir=".",
		     NormalizationMode NormalizationMode = NormalizationMode::Independent,
		     D pThreshold = 1e-10, int H = 4);

    /**
     * Disposes the allocated memory.
     */
    ~HaliteClustering();

    /**
     * Finds clusters in subspaces.
     */
    void findCorrelationClusters();

    /**
     * Gets the number of beta-clusters found.
     */
    size_t numBetaClusters() {
      return classifier->betaClusters.size();
    }

    /**
     * Gets the number of correlation clusters (merged).
     */
    int numCorrelationClusters() {
      return correlationClusters.size();
    }

    std::shared_ptr<Classifier<D>> getClassifier() {
      return classifier;
    }

    std::vector<CorrelationCluster>& getCorrelationClusters() {
      return correlationClusters;
    }

    /**
     * Gets the used counting tree.
     */
    stCountingTree<D> *getCalcTree() {
      return calcTree;
    }//end getCalcTree

    /**
     * Gets the time spent in the normalization.
     */
    clock_t getTimeNormalization() {
      return timeNormalization;
    }//end getTimeNormalization

  private:
    /**
     * Dimension of data
     */
    size_t DIM;

    /**
     * Size of data
     */
    size_t SIZE;

    /**
     * Time spent in the normalization.
     */
    clock_t timeNormalization;

    /**
     * Counting-tree pointer.
     */
    stCountingTree<D> *calcTree;


    /**
     * Number of grid levels in the counting tree.
     */
    int H;

    /**
     * Choose between hard and soft clustering.
     */
    bool hardClustering;

    /**
     * Defines the initial tree level to look for clusters.
     */
    int initialLevel;

    /**
     * Defines the convolution matrix (center and direct neighbours).
     */
    int centralConvolutionValue;
    int neighbourhoodConvolutionValue;

    /**
     * Defines the threshold used to spot a beta-cluster, based on the binomial probability.
     */
    D pThreshold;

    shared_ptr<Normalization<D>> normalization;
    shared_ptr<Classifier<D> > classifier;
    std::vector<CorrelationCluster> correlationClusters;


    /**
     * Merges beta-clusters that share some database space.
     */
    void mergeBetaClusters();

    int shouldMerge(BetaCluster<D>&  icl, BetaCluster<D>& j);

    int cost(BetaCluster<D>*  icl, BetaCluster<D>* j);

    int indCost(D n);

    cv::Mat inputPCA(const BetaCluster<D>* iCl, const BetaCluster<D>* jCl);

     /**
     * Calculates the cThreshold based on the Minimun Description Length (MDL) method.
     *
     * @param attributesRelevance Vector with the calculed relevances of each attribute.
     *
     */
    D calcCThreshold(const std::vector<D>& attributesRelevance);

    /**
     * Finds the best cut point position based on the MDL method.
     *
     * @param sortedRelevance Vector with the sorted relevances of each attribute.
     *
     */
    int minimumDescriptionLength(const std::vector<D>& sortedRelevance);

    /**
     * Walk through the counting tree applying the convolution matrix
     * to each cell in a defined level.
     *
     * @param level The counting tree level to be analyzed.
     *
     */
    int walkThroughConvolution(int level, stCell& betaClusterCenter, std::vector<stCell>& betaClusterCenterParents);

    /**
     * Applies the convolution matrix to a grid cell.
     *
     * @param cell The grid cell to apply the convolution matrix.
     * @param cellParents The cell parents.
     * @param level The cell level in the counting tree.
     *
     */
    int applyConvolution(stCell& cell, std::vector<stCell>& cellParents,
			 int level);

    /**
     * Finds the position of a cell in the data space.
     *
     * @param cell The analyzed cell.
     * @param cellParents The cell parents.
     * @param min Vector to receive the minimum position limits of cell in each dimension.
     * @param max Vector to receive the maximum position limits of cell in each dimension.
     * @param level The level of cell in the counting tree.
     *
     */
    void cellPosition(stCell& cell, std::vector<stCell>& cellParents,
		      std::vector<D>& min, std::vector<D>& max, int level);

    /**
     * Finds the position of a cell in the data space, regarding a dimension e_j.
     *
     * @param cell The analyzed cell.
     * @param cellParents The cell parents.
     * @param min The minimum position limit of cell in dimension e_j.
     * @param max The maximum position limit of cell in dimension e_j.
     * @param level The level of cell in the counting tree.
     * @param j The dimension to be considered.
     *
     */
    void cellPositionDimensionE_j(stCell& cell, std::vector<stCell>& cellParents,
				  D *min, D *max, int level, int j);

    /**
     * Finds the external face neighbour of cell in a determined dimension.
     *
     * @param dimIndex The dimension to look for the neighbour.
     * @param cell The analyzed cell.
     * @param neighbour Pointer to to the memory space where the found neighbour is returned.	
     * @param cellParents Vector with the cell parents.
     * @param neighbourParents Vector to receive the found neighbour parents.
     * @param level The level of cell in the counting tree.
     *
     * ps.: externalNeighbour assumes that cellParents and neighbourParents have the very
     *      same containt when the function is called.
     */
    int externalNeighbour(int dimIndex, stCell& cell, stCell* neighbour,
			  std::vector<stCell>& cellParents,
			  std::vector<stCell>& neighbourParents, int level);

    /**
     * Finds the internal face neighbour of cell in a determined dimension.
     *
     * @param dimIndex The dimension to look for the neighbour.
     * @param cell The analyzed cell.
     * @param neighbour Pointer to to the memory space where the found neighbour is returned.
     * @param cellParents Vector with the cell parents.	   
     * @param level The level of cell in the counting tree.
     *
     */
    int internalNeighbour(int dimIndex, stCell& cell, stCell* neighbour,
			  std::vector<stCell>& cellParents,
			  int level);

 
    /**
     * Normalize data points and insert them in the counting tree.
     * Method from the stFractalDimension class (adapted).
     *
     * @param PointSource Source of the database objects.
     * @param Normalization::Mode Determines how data will be normalized.
     *
     */
    void readData(PointSource<D>& data, NormalizationMode mode);

 
    size_t getCenter(size_t level);

  };//end HaliteClustering
  //----------------------------------------------------------------------------
} //end Halite namespace
#endif //__HALITECLUSTERING_H

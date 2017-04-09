/**********************************************************************
 * GBDI Arboretum - Copyright (c) 2002-2009 GBDI-ICMC-USP
 *
 *                           Homepage: http://gbdi.icmc.usp.br/arboretum
 **********************************************************************/
/* ====================================================================
 * The GBDI-ICMC-USP Software License Version 1.0
 *
 * Copyright (c) 2009 Grupo de Bases de Dados e Imagens, Instituto de
 * CiÍncias Matem·ticas e de ComputaÁ„o, University of S„o Paulo -
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
 *        de Dados e Imagens, Instituto de CiÍncias Matem·ticas e de
 *        ComputaÁ„o, University of S„o Paulo - Brazil (the Databases
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
 * This file demonstrates the usage of the HaliteClustering class.
 *
 * @version 1.0
 * @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
 * @author Agma Juci Machado Traina (agma@icmc.usp.br)
 * @author Christos Faloutsos (christos@cs.cmu.edu)
 * @author Caetano Traina Jr (caetano@icmc.usp.br)
 */
// Copyright (c) 2002-2009 GBDI-ICMC-USP

//------------------------------------------------------------------------------
// Halite
//------------------------------------------------------------------------------

#include <string.h>
#include <time.h>

#include "HaliteClustering.h"

#include "PointSource.h"

using namespace Halite;


// global variables
clock_t startTime;

/**
 * Initiates the measurement of run time.
 */
void initClock() {
	startTime=clock();
}//end initClock

/**
 * Prints the elapsed time.
 */
void printElapsed() { 
	printf("Elapsed time: %0.3lf sec.\n\n",(double)(clock()-startTime)/CLOCKS_PER_SEC);
}//end printElapsed

typedef double Dbl;
/**
 * Initiates the clustering process.
 *
 * @param argc Number of arguments received.
 * @param argv Array with the arguments received.
 * @return success (0), error (1).
 */
int main(int argc, char **argv) {
	initClock(); // initiates the meassurement of run time
	
	// first validations
	if (argc != 7) {
	  printf("Usage: Halite <hardClustering> <pThreshold> <H> <cache_points> <infile> <outfile>\n");
		return 1; //error
	}//end if
	
	if (atoi(argv[3]) < 2) {
		printf("Halite needs at least two resolution levels (H >= 2) to perform the clustering process.\n");
		return 1; //error
	}//end if
	
	char cachepoints=atoi(argv[4]);
	if (cachepoints<0 || cachepoints>1) {
		printf("Possible values for cache_points are 0 (do not cache) and 1 (cache data points in memory)\n");
		return 1;
	}//end if


	// opens/creates the used files
	FILE  *result;
	result=fopen(argv[6], "w");
	if (!result) {
		printf("Halite could not create the result file.\n");
		return 1; //error
	}//end if

	PointSource<Dbl>*  db;
	try {
	  db=new TextFilePointSource<Dbl>(argv[5]);
	} catch(std::exception& e) {
	  std::cout<<e.what()<<"\n";
		printf("'Halite could not open database file.\n");
		return 1;
	}
	size_t DIM=db->dimension();
	PointSource<Dbl>* memdb=NULL;
	PointSource<Dbl>* datasource=db;
	std::vector<Dbl*> objectsArray;
	if (cachepoints) { //unlimited memory
		for(db->restartIteration(); db->hasNext();) {
			Dbl* tmp = new Dbl[DIM];
			const Dbl* t=db->readPoint();
			for(unsigned int j=0; j<DIM; j++) {
				tmp[j]=t[j];
			}
			objectsArray.push_back(tmp);
		}//end for
		memdb=new ArrayOfPointersPointSource<Dbl>(&objectsArray[0],db->dimension(),objectsArray.size());
		datasource=memdb;
	}

	// creates an object of the class HaliteClustering
	HaliteClustering<Dbl> *sCluster = new HaliteClustering<Dbl>(*datasource, argv[1], ".", NormalizationMode::Independent, atof(argv[2]), atoi(argv[3]));		
	
	printf("The tree was built.\n");
	printElapsed(); // prints the elapsed time
	
	printf("Time spent in the normalization.\n");
	printf("Time: %0.3lf sec.\n\n",(Dbl)(sCluster->getTimeNormalization())/CLOCKS_PER_SEC);
	
	// looks for correlation clusters
	sCluster->findCorrelationClusters();
	
	printf("The clusters were found.\n");
	printElapsed(); // prints the elapsed time
	
	// mounts the result file


	unsigned int numCorrelationClusters=sCluster->numCorrelationClusters();

	shared_ptr<Classifier<Dbl>> classifier=sCluster->getClassifier();
	std::vector<CorrelationCluster>& correlationClusters=sCluster->getCorrelationClusters();
	// axes relevant to the found clusters
	for (unsigned int i=0; i<numCorrelationClusters; i++) {
	  CorrelationCluster& correlationCluster=correlationClusters[i];
	  fputs("ClusterResult",result);
	  fprintf(result,"%d",i+1);	
	  for (unsigned int j=0; j<DIM; j++) {
	    fputs(correlationCluster.relevantDimension[j] ? " 1" : " 0",result);
	  }//end for
	  fputs("\n",result); // writes the relevant axes to the current cluster in the result file
	}//end for
	
	// labels each point after the clusters found
	fputs("LABELING\n",result);
	const Dbl *onePoint;

	int point=0;
	for (datasource->restartIteration(); datasource->hasNext(); point++) {
		onePoint=datasource->readPoint();
		std::vector<int> clusters;
		classifier->assignToClusters(onePoint,std::back_inserter(clusters));
		for(unsigned int i=0; i<clusters.size(); i++) {
			fprintf(result, "%d %d\n", point+1, clusters[i]);
		}
		if(clusters.empty()) {
			fprintf(result, "%d %d\n", point+1, 0);
		}
	}//end for

	fclose(result); // the result file will not be used anymore
	delete db;

	printf("The result file was built.\n");
	printElapsed(); // prints the elapsed time

	// disposes the used structures
	if (cachepoints) { //unlimited memory
		// disposes objectsArray
		for (unsigned int i=0;i<objectsArray.size();i++) {
			delete [] objectsArray[i];
		}//end for	
		delete memdb;	
	}
	delete sCluster;
	
	printf("The used structures were deleted. Done!\n");
	printElapsed(); // prints the elapsed time
	//getc(stdin);
	return 0; // success
}//end main

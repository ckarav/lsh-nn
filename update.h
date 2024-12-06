#ifndef _UPDATE_
#define _UPDATE_

#include "enums.h"
#include "multi.h"

/************ Lloyd's update ***************/
//update a la Lloyd's (returns 1 if there was an actual update, else 0)
int Lloyd_update(METRIC_SPACES metricSpace, int *pointsClust, const multiMatrix *dataSet,
                  int rowNum, int colNum, int k, int *centroidsIndx, multiLongDist old_cost);

//finds the medoid of a cluster
int find_medoid(METRIC_SPACES metricSpace, int *pointsClust, const multiMatrix *dataSet,
                int rowNum, int colNum, int cluster);

//calculates local cost of a cluster relative to a point
multiLongDist local_cost(METRIC_SPACES metricSpace, int *pointsClust, const multiMatrix *dataSet,
                         int rowNum, int colNum, int point, int cluster);

/********************************************/
/***************** CLARANS ******************/
int Clarans(METRIC_SPACES metricSpace, int *pointClust, const multiMatrix *dataSet,
            int rowNum, int colNum, int k, int *centroidsIndx, multiLongDist old_cost, int claransIterations, int claransSetFraction);

//returns total cost of a configuration
multiLongDist config_cost(METRIC_SPACES metricSpace, int *pointClust, const multiMatrix *dataSet, int rowNum,
                          int colNum, int k, int *centroidsIndx);
/********************************************/
#endif // _UPDATE_

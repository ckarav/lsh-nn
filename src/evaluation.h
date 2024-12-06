#ifndef __EVALUATION__
#define __EVALUATION__

#include <stdio.h>
#include "enums.h"
#include "multi.h"

void silhouette(FILE* output, METRIC_SPACES metricSpace, const multiMatrix* datasetMatrix, int rowNum, int colNum, int clusterNum, int* pointsClust);
long double silhouette2(METRIC_SPACES metricSpace, const multiMatrix* datasetMatrix, int rowNum, int colNum, int clusterNum, int* pointsClust, int* centroidsIndx);

int pointsInClust(int *pointsClust, int clust, int rowNum); // Get the points in the specified cluster.

#endif

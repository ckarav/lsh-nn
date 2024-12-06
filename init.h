#ifndef __INIT__
#define __INIT__

#include "enums.h"
#include "multi.h"

void k_medoids_init(METRIC_SPACES metricSpace, const multiMatrix* datasetMatrix, int rowNum, int colNum, int clusterNum, int* centroidIndexes);

/*
 * Initialize the centroids for the clusters, using the the Concentration (Park-Jun) algorithm.
 * Instead of using a square 2D distance matrix, an upper triangular one is used to conserve memory.
 * Additionally, in case the metric space is the matrix space, then the distance matrix is not calculated again as it is already calculated in the dataset.
 */
void concentrate_init(METRIC_SPACES metricSpace, const multiMatrix* datasetMatrix, int rowNum, int colNum, int clusterNum, int* centroidIndexes);

#endif

#include <math.h>
#include "evaluation.h"
#include "update.h"
#include "distsNmore.h"
#include "ClusterElementsNum.h"

long double calcSi(METRIC_SPACES metricSpace, multiLongDist ai, multiLongDist bi, multiLongDist divz){
    long double Si;
    if(metricSpace == HAMMING_SPACE)
        Si = (long double)((bi.i - ai.i)/((double)(divz.i)));
    else
        Si = (long double)((bi.d - ai.d)/divz.d);
    return Si;
}


multiLongDist divByNum(METRIC_SPACES metricSpace, multiLongDist mlD, int num){
    if(metricSpace == HAMMING_SPACE)
        mlD.i /= num;
    else
        mlD.d /= (double)(num);
    return mlD;
}



int pointsInClust(int *pointsClust, int clust, int rowNum){
    int i, points;
    points = 0;
    for(i=0;i<rowNum;i++){
        if(pointsClust[i] == clust)
            points++;
    }
    points++;//for centroid
    return points;
}


void silhouette(FILE* output, METRIC_SPACES metricSpace, const multiMatrix* datasetMatrix, int rowNum, int colNum, int clusterNum, int* pointsClust) {
	long double a; // Average distance to other objects in the same cluster.
	long double b; // Min average distance to other cluster (next best or neighbor).
	long double s; // Silhouette metric for each object.
	long double sTotal = 0; // Average Silhouette metric for all the points in the dataset.
	ClusterElementsNumPtr array; // Stores the number of elements in each of the clusters, using the cluster's centroid.
	multiLongDist cost;
	int i, j, curIndex;

	// Allocate memory for the arrays.
	init(&array, pointsClust, rowNum, clusterNum); // Initialize the clusterElementsNum array.

	fprintf(output, "Silhouette: [");
	for (i = 0; i < rowNum; i++) {
		// Calculate the cost of the element i to the elements of its cluster.
		if (pointsClust[i] != -1) {
			cost = local_cost(metricSpace, pointsClust, datasetMatrix, rowNum, colNum, i, pointsClust[i]);
			if (metricSpace != HAMMING_SPACE)
				a = cost.d / (long double)(getClusterSize(array, pointsClust[i]) - 1);
			else
				a = (long double)cost.i / (long double)(getClusterSize(array, pointsClust[i]) - 1);
		} else {
			cost = local_cost(metricSpace, pointsClust, datasetMatrix, rowNum, colNum, i, i);
			if (metricSpace != HAMMING_SPACE)
				a = cost.d / (long double)(getClusterSize(array, i) - 1);
			else
				a = (long double)cost.i / (long double)(getClusterSize(array, i) - 1);
		}
		b = -1.0; // Set b[i] to a default value that will be overwritten by the min average distance to other clusters.

		for (j = 0; j < clusterNum; j++) {
			curIndex = getCentroidIndex(array, j); // Get the centroid of the j-th cluster of the ClusterElementsNum array.
			if ((curIndex == pointsClust[i]) || (curIndex == i)) continue; // Skip the average distance to the local cluster (it is already stored in a[i]).

			// Calculate the cost of the element i to the elements of the other clusters besides its own.
			cost = local_cost(metricSpace, pointsClust, datasetMatrix, rowNum, colNum, i, curIndex);

			// Calculate the average of the previous cost and override b[i] if the average is lower than it.
			if (metricSpace != HAMMING_SPACE) {
				cost.d = cost.d / (long double)getClusterSize(array, curIndex);
				if ((b == -1.0) || (cost.d < b)) b = cost.d;
			} else {
				cost.i = (long double)cost.i / (long double)getClusterSize(array, curIndex);
				if ((b == -1.0) || ((long double)cost.i < b)) b = (long double)cost.i;
			}
		}

		s = (b - a) / (long double)((a > b) ? a : b); // Calculate the Silhouette metric for the i-th element.
		fprintf(output, "%3.10Lf,", s);
		sTotal += s;
	}
	sTotal = sTotal / (long double)rowNum; // Calculate the average Silhouette metric for all the points in the dataset.

	fprintf(output, "%3.10Lf]\n", sTotal);

	// Deallocate memory.
	destroy(&array);
}


long double silhouette2(METRIC_SPACES metricSpace, const multiMatrix* datasetMatrix, int rowNum, int colNum, int clusterNum, int* pointsClust, int* centroidsIndx) {
	long double si, avgS, tmpS;
	multiLongDist ai, bi, minBi, divz;
	int iter;
	int i, j;

	si = 0;
	for(i=0;i<rowNum;i++){
		//calculate ai
		if(pointsClust[i] == (-1)){//if centroid
			ai = local_cost(metricSpace,pointsClust,datasetMatrix,rowNum,colNum,i,i);
			ai = divByNum(metricSpace,ai,pointsInClust(pointsClust,i,rowNum));//get avg
		}
		else{
			ai = local_cost(metricSpace,pointsClust,datasetMatrix,rowNum,colNum,i,pointsClust[i]);
			ai = divByNum(metricSpace,ai,pointsInClust(pointsClust,pointsClust[i],rowNum));
		}
		//calculate bi
		iter = 0;
		for(j=0;j<clusterNum;j++){
			if((centroidsIndx[j] == pointsClust[i]) || (centroidsIndx[j] == i))//skip same cluster
				continue;
			if(iter == 0){//first real iteration
				minBi = local_cost(metricSpace,pointsClust,datasetMatrix,rowNum,colNum,i,centroidsIndx[j]);
				minBi = divByNum(metricSpace,minBi,pointsInClust(pointsClust,centroidsIndx[j],rowNum));
				iter++;
				continue;
			}
			bi = local_cost(metricSpace,pointsClust,datasetMatrix,rowNum,colNum,i,centroidsIndx[j]);
			bi = divByNum(metricSpace,bi,pointsInClust(pointsClust,centroidsIndx[j],rowNum));
			if(lessThan2(metricSpace,bi,minBi))
				minBi = bi;
		}
		/////////////////
		//calculate s(i)
		if(lessThan2(metricSpace,ai,bi))
			divz = bi;
		else
			divz = ai;
		tmpS = calcSi(metricSpace,ai,bi,divz);
		si += calcSi(metricSpace,ai,bi,divz);
		//printf("%3.10Lf\n", tmpS);
	}
	avgS = si/(long double)rowNum;
	return avgS;
}

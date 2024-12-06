#ifndef _ASSIGNMENT_
#define _ASSIGNMENT_

#include "enums.h"
#include "multi.h"

//pointsClust is an array of size Nx2 that holds cluster index(and second nearest) of corresponding point,
//e.g. pointsClust[2][0] = 3 means 2nd point belongs to 3rd cluster

multiLongDist PAM_assignment(METRIC_SPACES metricSpace, int *pointsClust, const multiMatrix *dataSet,
		int rowNum, int colNum, int k, int *centroidsIndx);

//returns nearest centroid to point
multiDist nearestCentroid(METRIC_SPACES metricSpace, const multiMatrix *dataSet, int point,
                    int k, int *centroidsIndx, int colNum, int *nearest_cent);

/***********************************************/
/************* LSH/DBH functions ***************/
multiLongDist LSH_DBH_assign(METRIC_SPACES metricSpace, int *pointsClust, const multiMatrix *dataSet,
        int rowNum, int colNum, int L, int k, int *centroidsIndx, multiHashes *mHs, multiTables *mTs);
//returns min(dist between centers)/2 and stores max(dist between centers)/2 in threshold
double initRadius(METRIC_SPACES metricSpace, int k, int *centroidsIndx, const multiMatrix *dataSet,
                  int colNum, double *threshold);

//Range_Search around centroid and update points assigned accordingly
int Range_Search(METRIC_SPACES metricSpace, int L, multiHashes *mHs, multiTables *mTs, int iter,
                  double R, int colNum, int centroid, const multiMatrix *dataSet, int *state_array, int *pointsClust);
/////////////////////////////////////
//initializes hash tables of any metric space (Ris used only in case of vector space)
void initStructs(int L, int k, METRIC_SPACES metricSpace, multiHashes *mHs, multiTables *mTs,
                 const multiMatrix *dataSet, int colNum, int rowNum);

//destroys hash tables of any matric space
void destroyStructs(int L, METRIC_SPACES metricSpace, multiHashes *mHs, multiTables *mTs);


/***********************************************/

#endif // _ASSIGNMENT_


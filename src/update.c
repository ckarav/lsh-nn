#include "update.h"
#include "distsNmore.h"
#include "assignment.h"
#include "evaluation.h"
#include <math.h>
#include "rands.h"
#include <stdlib.h>

int Clarans(METRIC_SPACES metricSpace, int *pointClust, const multiMatrix *dataSet,
            int rowNum, int colNum, int k, int *centroidsIndx, multiLongDist old_cost, int claransIterations, int claransSetFraction){
    ////////////////////////////
    int i, old_clust, update, q, x, old_medoids_clust;
    int m, t; //pair (m,t) to be swapped
    multiLongDist new_cost;

    update = 0;
    q = claransSetFraction;
    q *= claransIterations; //s = 2 means PAM swap with twice the amount of pairs (m,t)
    ///////////////////////////
    i = 0;
    while(i < q){
        x = uniform_randInt(0,k*rowNum-1);//uniform int in range [0,k*N-1]
        m = x % k;
        t = x / k;
        if(pointClust[t] == (-1))//can't swap 2 centroids
            continue;
        //temporary change in configuration
        old_clust = centroidsIndx[m];
        centroidsIndx[m] = t;
        pointClust[old_clust] = 1;//anything >=0
        old_medoids_clust = pointClust[t];
        pointClust[t] = -1;
        ///////////////////
        new_cost = config_cost(metricSpace,pointClust,dataSet,rowNum,colNum,k,centroidsIndx);
        if(lessThan2(metricSpace,new_cost,old_cost)){//keep updated configuration
            update = 1;
            old_cost = new_cost;
        }
        else{//keep old configuration
            centroidsIndx[m] = old_clust;
            pointClust[old_clust] = -1;
            pointClust[t] = old_medoids_clust;
        }
        ///////////////////
        i++;
    }
    return update;
}

multiLongDist config_cost(METRIC_SPACES metricSpace, int *pointClust, const multiMatrix *dataSet, int rowNum,
                          int colNum, int k, int *centroidsIndx){
    /////////////////
    int i, nearest_cent;
    multiLongDist total_cost;
    multiDist tempDist;
    zeroInit(metricSpace,&total_cost);
    ////////////////////
    for(i=0;i<rowNum;i++){
        if(pointClust[i] == (-1))//skip centroids
            continue;
        tempDist = nearestCentroid(metricSpace,dataSet,i,k,centroidsIndx,colNum,&nearest_cent);
        total_cost = addDist(metricSpace,tempDist,total_cost);
    }
    return total_cost;
}

/*************************************************************/
int Lloyd_update(METRIC_SPACES metricSpace, int *pointsClust, const multiMatrix *dataSet,
                  int rowNum, int colNum, int k, int *centroidsIndx, multiLongDist old_cost){
    /////////////////////////
    int j, medoid, old_clust, update, old_medoids_clust;
    multiLongDist new_cost;
    update = 0;
    ////////////////////////
    for(j=0;j<k;j++){
    	if (pointsInClust(pointsClust, centroidsIndx[j], rowNum) == 1) continue;

        zeroInit(metricSpace,&new_cost);//init total cost with 0
        old_clust = centroidsIndx[j];
        medoid = find_medoid(metricSpace,pointsClust,dataSet,rowNum,colNum,old_clust);
        centroidsIndx[j] = medoid;
        ///////////////////////////////////
        //create new configuration
        //for(i=0;i<k;i++) PENDING
        old_medoids_clust = pointsClust[medoid];
        pointsClust[medoid] = -1;
        pointsClust[old_clust] = 1; //anything >=0
        new_cost = config_cost(metricSpace,pointsClust,dataSet,rowNum,colNum,k,centroidsIndx);
        if(lessThan2(metricSpace,new_cost,old_cost)){
            update = 1;
            old_cost = new_cost;
        }
        else{
            centroidsIndx[j] = old_clust;
            pointsClust[old_clust] = -1;
            pointsClust[medoid] = old_medoids_clust;//anything >=0
        }
    }
    ///////////////////////////
    return update;
}


int find_medoid(METRIC_SPACES metricSpace, int *pointsClust, const multiMatrix *dataSet,
                int rowNum, int colNum, int cluster){
    //////////////////////////////////
    int i, medoid;
    multiDist mD;
    multiLongDist bestDist, tempDist;
    i = 0;
    /////////////////////////////////
    while(pointsClust[i] != cluster)
        i++;
    //init medoid
    medoid = i; /*initialize medoid with first point of cluster*/
    bestDist = local_cost(metricSpace,pointsClust,dataSet,rowNum,colNum,i,cluster);
    ///////////
    i++;
    while(1){
        while((i<rowNum) && (pointsClust[i] != cluster))
            i++;
        if(i == rowNum) break;//we checked all points of cluster
        //find cluster cost with respect to point 'i'
        tempDist = local_cost(metricSpace,pointsClust,dataSet,rowNum,colNum,i,pointsClust[i]);
        if(lessThan2(metricSpace,tempDist,bestDist)){
            bestDist = tempDist;
            medoid = i;
        }
        i++;
    }
    return medoid; //return point index of cluster's medoid
}


multiLongDist local_cost(METRIC_SPACES metricSpace, int *pointsClust, const multiMatrix *dataSet,
                         int rowNum, int colNum, int point, int cluster){
    ////////////////////////
    int i;
    multiLongDist mlD;
    multiDist tempDist;

    zeroInit(metricSpace,&mlD);//init local cost with 0
    ////////////////////////
    for(i=0;i<rowNum;i++){
        if((i == point) || (pointsClust[i] != cluster))
            continue;
        tempDist = general_distance(metricSpace,point,i,dataSet,colNum);
        mlD = addDist(metricSpace,tempDist,mlD);
    }
    //calculate distance of point to cluster's centroid
    tempDist = general_distance(metricSpace,point,cluster,dataSet,colNum);
    mlD = addDist(metricSpace,tempDist,mlD);
    return mlD;
}

















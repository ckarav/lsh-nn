#include "assignment.h"
#include "distsNmore.h"
#include "hamming.h"
#include "matrix.h"
#include "vector.h"


multiLongDist PAM_assignment(METRIC_SPACES metricSpace, int *pointsClust, const multiMatrix *dataSet,
		int rowNum, int colNum, int k, int *centroidsIndx){
    ///////////////////////
    int i, nearest_cent;
    multiDist minDist;
    multiLongDist total_cost;
    zeroInit(metricSpace,&total_cost);//init total cost with 0
    ////////////////////
    for(i=0;i<rowNum;i++){
        if(pointsClust[i] < 0) //if (-1) then ith point is a centroid, so skip assignment
            continue;
        minDist = nearestCentroid(metricSpace,dataSet,i,k,centroidsIndx,colNum,&nearest_cent);
        pointsClust[i] = nearest_cent; //assign point to its nearest centroid
        total_cost = addDist(metricSpace,minDist,total_cost); //update total cost
    }
    return total_cost; //return clustering cost
}


multiDist nearestCentroid(METRIC_SPACES metricSpace, const multiMatrix *dataSet, int point,
                    int k, int *centroidsIndx, int colNum, int *nearest_cent){
    ////////////////////////////
    int j;
    multiDist minDist, tempDist;
    ///////
    for(j=0;j<k;j++){
        if(j == 0){ //case of first comparison with centroid
            minDist = general_distance(metricSpace,point,centroidsIndx[j],dataSet,colNum);
            (*nearest_cent) = j;
            continue;
        }
        tempDist = general_distance(metricSpace,point,centroidsIndx[j],dataSet,colNum);
        if(lessThan(metricSpace,tempDist,minDist)){
            minDist = tempDist; //updating best centroid
            (*nearest_cent) = j;
        }
    }
    (*nearest_cent) = centroidsIndx[(*nearest_cent)];
    return minDist;
}


/*******************************************************************/
/*********************** LSH/DBH functions *************************/

multiLongDist LSH_DBH_assign(METRIC_SPACES metricSpace, int *pointsClust, const multiMatrix *dataSet,
        int rowNum, int colNum, int L, int k, int *centroidsIndx, multiHashes *mHs, multiTables *mTs){
    ////////////////////////////////////////
    int i, j, flag, empty_ranges, *state_array, range_iter, nearest_cent;
    double R, threshold;
    multiDist minDist, tempDist;
    multiLongDist total_cost;

    range_iter = 0;//shows the number of range searches completed
    flag = 1;//indication for when to stop
    R = initRadius(metricSpace,k,centroidsIndx,dataSet,colNum,&threshold);
    //state_array tells us if a point is assigned and during which range search it got assigned
    state_array = malloc(rowNum*sizeof(int));//allocate space for state_array
    for(i=0;i<rowNum;i++){//initialize state_array
        if(pointsClust[i] == -1)//if centroid
            state_array[i] = -1;
        else//normal point
            state_array[i] = -2; //-2 encodes "not assigned"
    }
    zeroInit(metricSpace,&total_cost);//init total cost with 0
    //////////////////////
    while(flag){
        empty_ranges = 0;//condition to stop range searches if they're empty
        for(j=0;j<k;j++) //for each centroid
            empty_ranges += Range_Search(metricSpace,L,mHs,mTs,range_iter,R,colNum,centroidsIndx[j],dataSet,state_array,pointsClust);
        R *= 2; //increase radius for next search
        if((R > threshold) || (empty_ranges > k/2)) //condition to stop searches
            flag = 0;
        range_iter++;
    }
    /** Assign remaining(unassigned) points **/
    for(i=0;i<rowNum;i++){
        if(state_array[i] == (-1)) //if (-1) then ith point is a centroid, so skip assignment
            continue;
        //for each unassigned point, compare distances to all centroids
        if(state_array[i] == (-2)){
            for(j=0;j<k;j++){
                if(j == 0){ //case of firt comparison with centroid
                    minDist = general_distance(metricSpace,i,centroidsIndx[j],dataSet,colNum);
                    nearest_cent = j;
                    continue;
                }
                tempDist = general_distance(metricSpace,i,centroidsIndx[j],dataSet,colNum);
                if(lessThan(metricSpace,tempDist,minDist)){
                    minDist = tempDist; //updating best centroid
                    nearest_cent = j;
                }
            }
            pointsClust[i] = centroidsIndx[nearest_cent];
            total_cost = addDist(metricSpace,minDist,total_cost); //update total cost
        }
        else{
            tempDist = general_distance(metricSpace,i,centroidsIndx[pointsClust[i]],dataSet,colNum);
            total_cost = addDist(metricSpace,tempDist,total_cost);
        }
    }
    /////////////////////////
    free(state_array);
    return total_cost;
}

int Range_Search(METRIC_SPACES metricSpace, int L, multiHashes *mHs, multiTables *mTs, int iter,
                  double R, int colNum, int centroid, const multiMatrix *dataSet, int *state_array, int *pointsClust){
    /////////////////////
    int objID, empty_list;
    empty_list = 0;
    /////////////////////
    if(metricSpace == HAMMING_SPACE){
        ham_node *temp;
        ham_list range_results;
        hamList_init(&range_results);//init results list for hamming
        /////////////////////////////////
        range_ham2(L,mTs->HTables,mHs->Gh,R,(dataSet->c)[centroid],&range_results,state_array,iter);
        if(range_results.N == 0)
            empty_list = 1;
        temp = range_results.head;
        while(temp){
            objID = temp->item_index;
            if(state_array[objID] == (-2)){//previously unassigned point
                state_array[objID] = iter;//declare it was assigned during this iteration
                pointsClust[objID] = centroid;//assign point to centroid
            }
            else{//point assigned to previous cluster during this iteration so compare distances
                if(hamming_distance((dataSet->c)[objID],(dataSet->c)[centroid],colNum) < hamming_distance((dataSet->c)[objID],(dataSet->c)[pointsClust[objID]],colNum))
                    pointsClust[objID] = centroid;//assign point to current centroid
            }
            hamList_clear(temp);//remove same node from rest of list
            temp = temp->next;
        }
        /////////////////////////////////
        hamList_destroy(&range_results);
    }
    else if(metricSpace == VECTOR_SPACE){
        vec_list range_results;
        vec_node *temp;
        vecList_init(&range_results);
        /////////////////////
        range_vector2(L,mTs->VTables,mHs->Gv,R,(dataSet->d)[centroid],colNum,&range_results,state_array,iter);
        if(range_results.N == 0)
            empty_list = 1;
        temp = range_results.head;
        while(temp){
            objID = temp->item_index;
            if(state_array[objID] == (-2)){
                state_array[objID] = iter;
                pointsClust[objID] = centroid;
            }
            else{
                if(euclidean_distance((dataSet->d)[objID],(dataSet->d)[centroid],colNum)<euclidean_distance((dataSet->d)[objID],(dataSet->d)[pointsClust[objID]],colNum))
                    pointsClust[objID] = centroid;
            }
            vecList_clear(temp);
            temp = temp->next;
        }
        ///////////////////
        vecList_destroy(&range_results); //destroy results list
    }
    else{ //MATRIX_SPACE
        matrix_list range_results;
        matrix_node *temp;
        matList_init(&range_results);
        ///////////////////////
        range_matrix2(L,mTs->MTables,mHs->Gm,R,(dataSet->d)[centroid],&range_results,state_array,iter);
        if(range_results.N == 0)
            empty_list = 1;
        temp = range_results.head;
        while(temp){
            objID = temp->itemID;
            if(state_array[objID] == (-2)){
                state_array[objID] = iter;
                pointsClust[objID] = centroid;
            }
            else{
                if(((dataSet->d)[objID][centroid]) < ((dataSet->d)[objID][pointsClust[objID]]))
                    pointsClust[objID] = centroid;
            }
            matList_clear(temp);
            temp = temp->next;
        }
        ///////////////////////
        matList_destroy(&range_results);
    }
    return empty_list;
}

double initRadius(METRIC_SPACES metricSpace, int k, int *centroidsIndx, const multiMatrix *dataSet,
                  int colNum, double *threshold){
    int i, j;
    multiDist minDist, temp, maxDist;
    ////////////////
    for(i=0;i<k;i++){
        for(j=i;j<k;j++){ //find minimum distance between centtoids
            if((i==0) && (j==0)){ //case of first iteration, init min
                minDist = general_distance(metricSpace,centroidsIndx[i],centroidsIndx[j],dataSet,colNum);
                maxDist = minDist;
                continue;
            }
            temp = general_distance(metricSpace,centroidsIndx[i],centroidsIndx[j],dataSet,colNum);
            if(lessThan(metricSpace,temp,minDist))//
                minDist = temp;
            if(lessThan(metricSpace,maxDist,temp))
                maxDist = temp;
        }
    }
    ///////////////////////////////
    //divide distance by 2
    if(metricSpace == HAMMING_SPACE){
        (*threshold) = (double)(maxDist.i / 2);
        return ((double)(minDist.i / 2));
    }
    else{
        (*threshold) = maxDist.d / 2.0;
        return (minDist.d / 2.0);
    }
}

////////////////////
void initStructs(int L, int k, METRIC_SPACES metricSpace, multiHashes *mHs, multiTables *mTs,
                 const multiMatrix *dataSet, int colNum, int rowNum){
    //////////////////////////////////
    if(metricSpace == HAMMING_SPACE)
        InitHam(L,k,&(mTs->HTables),&(mHs->Gh),dataSet->c,colNum,rowNum);
    else if(metricSpace == VECTOR_SPACE)
        InitVec(L,k,&(mTs->VTables),&(mHs->Gv),dataSet->d,colNum,rowNum);
    else if(metricSpace == COSINE)
    	InitCosine(L,k,&(mTs->CTables),&(mHs->Gc),dataSet->d,rowNum,colNum);
    else //MATRIX_SPACE
        InitMatrix(L,k,&(mTs->MTables),&(mHs->Gm),dataSet->d,colNum);
}


void destroyStructs(int L, METRIC_SPACES metricSpace, multiHashes *mHs, multiTables *mTs){
    //////////////////////////////////////
    if(metricSpace == HAMMING_SPACE)
        FreeHam(L,mHs->Gh,mTs->HTables);
    else if(metricSpace == VECTOR_SPACE)
        FreeVec(L,mTs->VTables,mHs->Gv);
    else if(metricSpace == COSINE)
    	FreeCosine(L,mTs->CTables,mHs->Gc);
    else //MATRIX_SPACE
        FreeMatrix(L,mTs->MTables,mHs->Gm);
}

#include <stdio.h>
#include <stdlib.h>
#include "distsNmore.h"
#include "recommend.h"

void cluster_recommend(METRIC_SPACES metricSpace, int numOfUsers, int numOfItems,
		double **dataSet, int **five_best, double *avgRating, int *pointsClust){
	////////////////
	int i, j, k, my_clust;
	double z, sum, rating;
	/////
	for(i=0;i<numOfUsers;i++){
		//find the cluster this user is in
		if(pointsClust[i] == (-1))
			my_clust = i;
		else
			my_clust = pointsClust[i];
		z = calculateZ2(metricSpace,my_clust,i,numOfItems,numOfUsers,pointsClust,dataSet);
		for(k=0;k<numOfItems;k++){
			if(dataSet[i][k] != 0.0)
				continue;
			sum = 0.0;
			for(j=0;j<numOfUsers;j++){
				if((j == i) || ((pointsClust[j]!=(-1)) && (pointsClust[j]!=my_clust)))
					continue;
				if(((pointsClust[j]==(-1)) && (j==my_clust)) || ((pointsClust[j]!=(-1)) && (pointsClust[j]==my_clust))){
					if(dataSet[j][k] != 0.0){
						if(metricSpace == VECTOR_SPACE)
							sum += euclidean_distance(dataSet[i],dataSet[j],numOfItems) *(dataSet[j][k]-avgRating[j]);
						else
							sum += cosine_dist(dataSet[i],dataSet[j],numOfItems) *(dataSet[j][k]-avgRating[j]);
					}
				}
			}
			rating = z * sum;
			dataSet[i][k] = rating;
			addToTop5Matrix(i,dataSet,five_best[i],k,rating);
		}
	}

}


void NN_lsh_recommend(METRIC_SPACES metricSpace, int L, multiHashes *mHs, multiTables *mTs, int numOfUsers,
		int numOfItems, double **dataSet, int P, int **five_best, double *avgRating){
	/////////////////////////////
	int i, j, k;
	double z, sum, rating;
	neigh_struct **neighbours;
	init_neighMatrix(&neighbours,P,numOfUsers);
	////////
	/*Calculate P closest neighbors of each user*/
	for(i=0;i<numOfUsers;i++){
		if(metricSpace == VECTOR_SPACE)
			binary_repeated_search(L,mHs,mTs,numOfUsers,numOfItems,dataSet,neighbours,i,P);
		else
			binary_repeated_search2(L,mHs,mTs,numOfUsers,numOfItems,dataSet,neighbours,i,P);
		printf("COMPLETED USER %d\n", i);
	}
	/********************************************/
	for(i=0;i<numOfUsers;i++){
		z = calculateZ(neighbours[i],P);
		for(k=0;k<numOfItems;k++){
			if(dataSet[i][k] != 0.0)//item is already rated
				continue;
			sum = 0.0;
			for(j=0;j<P;j++){
				if(neighbours[i][j].neighbour_ID != (-1)){
					if(dataSet[neighbours[i][j].neighbour_ID][k] != 0.0)//neighbor has rated corresponding item
						sum += neighbours[i][j].similarity * (dataSet[neighbours[i][j].neighbour_ID][k] - avgRating[neighbours[i][j].neighbour_ID]);
				}
			}
			rating = z * sum;
			dataSet[i][k] = rating;
			addToTop5Matrix(i,dataSet,five_best[i],k,rating);
		}
	}
	////////////////////////////
	free_neighMatrix(neighbours,numOfUsers);
}

double maxUsrDist(METRIC_SPACES metricSpace, double **dataSet, int numOfUsers, int numOfItems, int userID){
	int i;
	double max, temp;
	max = 0.0;
	for(i=0;i<numOfUsers;i++){
		if(metricSpace == VECTOR_SPACE)
			temp = euclidean_distance(dataSet[userID],dataSet[i],numOfItems);
		else
			temp = cosine_dist(dataSet[userID],dataSet[i],numOfItems);
		if(temp > max )
			max = temp;
	}
	return max;
}


void binary_repeated_search(int L, multiHashes *mHs, multiTables *mTs, int numOfUsers,
		int numOfItems, double **dataSet, neigh_struct **neighbours, int userID, int P){
	///////////
	int i, temp_user;
	double R_min, R_max, R_mid, R0, max_dist, temp_dist, max_usr_dist;
	vec_list range_results;
	vec_node *temp;
	//////////////////////////////////
	/*initialize R0*/
	for(i=0;i<10;i++){
		temp_user = rand() % numOfUsers;
		if(i == 0){
			max_dist = euclidean_distance(dataSet[userID],dataSet[temp_user],numOfItems);
			continue;
		}
		temp_dist = euclidean_distance(dataSet[userID],dataSet[temp_user],numOfItems);
		if(temp_dist > max_dist)
			max_dist = temp_dist;
	}
	R0 = max_dist / 2.0;
	max_usr_dist = maxUsrDist(VECTOR_SPACE,dataSet,numOfUsers,numOfItems,userID);
	////////////////
	/*calculate R_min/max*/
	vecList_init(&range_results);
	range_vector(L,mTs->VTables,mHs->Gv,R0,dataSet[userID],numOfItems,&range_results);
	clearVecList(&range_results);//to remove duplicates
	if(range_results.N > (P+5)){
		R_min = R0;
		while(range_results.N > (P+5)){
			vecList_destroy(&range_results);
			vecList_init(&range_results);
			R_min /= 2.0;//halve the range
			range_vector(L,mTs->VTables,mHs->Gv,R_min,dataSet[userID],numOfItems,&range_results);
			clearVecList(&range_results);//to remove duplicates
		}
		R_max = R_min * 2.0;
	}
	else if(range_results.N < P){
		R_max = R0;
		while(range_results.N < P){
			vecList_destroy(&range_results);
			vecList_init(&range_results);
			R_max *= 2.0;//double the range
			if(R_max > max_usr_dist){
				range_vector_yolo(L,mTs->VTables,mHs->Gv,R_max,dataSet[userID],numOfItems,&range_results, 1);
				break;
			} else {
				range_vector_yolo(L,mTs->VTables,mHs->Gv,R_max,dataSet[userID],numOfItems,&range_results, 0);
			}
			clearVecList(&range_results);//to remove duplicates
		}
		R_min = R_max / 2.0;
	}
	else{//highly unlikely case of exactly P neighbors with randomized R0
		//simply fill neighbours array by copying from results list
		temp = range_results.head;
		for(i=0;i<P;i++){
			neighbours[userID][i].neighbour_ID = temp->item_index;
			neighbours[userID][i].similarity = euclidean_distance(dataSet[userID],dataSet[temp->item_index],numOfItems);
			temp = temp->next;
		}
		return;
	}
	///////////////
	/*Binary Search to find best range*/
	R_mid = (R_min + R_max) / 2.0;
	while((range_results.N < P) || (range_results.N > (P+5))){
		if(range_results.N < P){
			R_min = R_mid;
			R_mid = (R_min + R_max) / 2.0;
		}
		else if(range_results.N > (P+5)){
			R_max = R_mid;
			R_mid = (R_min + R_max) / 2.0;
		}
		if ((R_max - R_min) < 2.0){
			if(range_results.N > (P+5))
				break;
			else{
				vecList_destroy(&range_results);
				vecList_init(&range_results);
				R_mid *= 2.0;
				range_vector(L,mTs->VTables,mHs->Gv,R_mid,dataSet[userID],numOfItems,&range_results);
				clearVecList(&range_results);//to remove duplicates
				while(range_results.N < P){
					R_mid *= 2.0;
					vecList_destroy(&range_results);
					vecList_init(&range_results);
					if(R_mid > max_usr_dist){
						range_vector_yolo(L,mTs->VTables,mHs->Gv,R_max,dataSet[userID],numOfItems,&range_results, 1);
						break;
					} else {
						range_vector_yolo(L,mTs->VTables,mHs->Gv,R_max,dataSet[userID],numOfItems,&range_results, 0);
					}
					clearVecList(&range_results);//to remove duplicates
				}
				break;
			}
		}
		vecList_destroy(&range_results);
		vecList_init(&range_results);
		range_vector(L,mTs->VTables,mHs->Gv,R_mid,dataSet[userID],numOfItems,&range_results);
		clearVecList(&range_results);//to remove duplicates
	}
	//////////////
	/*Fill neighbours matrix*/
	temp = range_results.head;
	while(temp){
		addToNeighMatrix(P,userID,dataSet,neighbours[userID],numOfItems,temp);//will add P closest neighbors
		temp = temp->next;
	}
	////////////////
	vecList_destroy(&range_results);
	/////////////////////////////////
}


void binary_repeated_search2(int L, multiHashes *mHs, multiTables *mTs, int numOfUsers,
		int numOfItems, double **dataSet, neigh_struct **neighbours, int userID, int P){
	///////////
	int i, temp_user;
	double R_min, R_max, R_mid, R0, max_dist, temp_dist;
	cos_vec_list range_results;
	cos_vec_node *temp;
	//////////////////////////////////
	/*initialize R0*/
	for(i=0;i<10;i++){
		temp_user = rand() % numOfUsers;
		if(i == 0){
			max_dist = cosine_dist(dataSet[userID],dataSet[temp_user],numOfItems);
			continue;
		}
		temp_dist = cosine_dist(dataSet[userID],dataSet[temp_user],numOfItems);
		if(temp_dist > max_dist)
			max_dist = temp_dist;
	}
	R0 = max_dist / 2.0;
	////////////////
	/*calculate R_min/max*/
	cosList_init(&range_results);
	cos_vector(L,mTs->CTables,mHs->Gc,R0,dataSet[userID],numOfItems,&range_results);
	clearCosList(&range_results);//to remove duplicates
	if(range_results.N > (P+5)){
		R_min = R0;
		while(range_results.N > (P+5)){
			cosList_destroy(&range_results);
			cosList_init(&range_results);
			R_min /= 2.0;//halve the range
			cos_vector(L,mTs->CTables,mHs->Gc,R_min,dataSet[userID],numOfItems,&range_results);
			clearCosList(&range_results);//to remove duplicates
		}
		R_max = R_min * 2.0;
	}
	else if(range_results.N < P){
		R_max = R0;
		while(range_results.N < P){
			cosList_destroy(&range_results);
			cosList_init(&range_results);
			R_max *= 2.0;//double the range
			cos_vector(L,mTs->CTables,mHs->Gc,R_max,dataSet[userID],numOfItems,&range_results);
			clearCosList(&range_results);//to remove duplicates
		}
		R_min = R_max / 2.0;
	}
	else{//highly unlikely case of exactly P neighbors with randomized R0
		//simply fill neighbours array by copying from results list
		temp = range_results.head;
		for(i=0;i<P;i++){
			neighbours[userID][i].neighbour_ID = temp->item_index;
			neighbours[userID][i].similarity = cosine_dist(dataSet[userID],dataSet[temp->item_index],numOfItems);
			temp = temp->next;
		}
		return;
	}
	///////////////
	/*Binary Search to find best range*/
	R_mid = (R_min + R_max) / 2.0;
	while((range_results.N < P) || (range_results.N > (P+5))){
		if(range_results.N < P){
			R_min = R_mid;
			R_mid = (R_min + R_max) / 2.0;
		}
		else if(range_results.N > (P+5)){
			R_max = R_mid;
			R_mid = (R_min + R_max) / 2.0;
		}
		if((R_max - R_min) < 0.5){
			if(range_results.N > (P+5))
				break;
			else{
				cosList_destroy(&range_results);
				cosList_init(&range_results);
				R_mid *= 2.0;
				cos_vector(L,mTs->CTables,mHs->Gc,R_mid,dataSet[userID],numOfItems,&range_results);
				clearCosList(&range_results);//to remove duplicates
				while(range_results.N < P){
					R_mid *= 2.0;
					cosList_destroy(&range_results);
					cosList_init(&range_results);
					cos_vector(L,mTs->CTables,mHs->Gc,R_mid,dataSet[userID],numOfItems,&range_results);
					if(R_mid > 2.0)
						break;
					clearCosList(&range_results);//to remove duplicates
				}
				break;
			}
		}
		cosList_destroy(&range_results);
		cosList_init(&range_results);
		cos_vector(L,mTs->CTables,mHs->Gc,R_mid,dataSet[userID],numOfItems,&range_results);
		clearCosList(&range_results);//to remove duplicates
	}
	//////////////
	/*Fill neighbours matrix*/
	temp = range_results.head;
	while(temp){
		addToNeighMatrix2(P,userID,dataSet,neighbours[userID],numOfItems,temp);//will add P closest neighbors
		temp = temp->next;
	}
	////////////////
	cosList_destroy(&range_results);
	/////////////////////////////////
}


void addToTop5Matrix(int userID, double **dataSet, int *five_best, int itemID, double rating){
	///////////////
	int i, lowest_rated_item, lowest_rating;
	///
	for(i=0;i<5;i++){
		if(five_best[i] == (-1))
			break;
		if(i==0){
			lowest_rated_item = 0;
			lowest_rating = dataSet[userID][five_best[0]];
		}
		else{
			if(dataSet[userID][five_best[i]] < lowest_rating){
				lowest_rating = dataSet[userID][five_best[i]];
				lowest_rated_item = i;
			}
		}
	}
	if(i == 5){//no free spot, maybe replace with lowest rated item
		if(rating > lowest_rating)
			five_best[lowest_rated_item] = itemID;
	}
	else//assign free spot to newly rated item
		five_best[i] = itemID;
}

void addToNeighMatrix(int P, int userID, double **dataSet, neigh_struct *neighArray, int numOfItems,
		vec_node *neighbor){
	////////////////////
	int i, farthest_neigh;
	double max, dist;
	////
	for(i=0;i<P;i++){
		if(neighArray[i].neighbour_ID == (-1))
			break;
		if(i == 0){
			max = neighArray[i].similarity;
			farthest_neigh = 0;
		}
		else{
			dist = neighArray[i].similarity;
			if(dist > max){
				max = dist;
				farthest_neigh = i;
			}
		}
	}
	if(i == P){//no free spot, maybe replace with farthest neighbor
		dist = euclidean_distance(dataSet[userID],dataSet[neighbor->item_index],numOfItems);
		if(dist < max){
			neighArray[farthest_neigh].neighbour_ID = neighbor->item_index;
			neighArray[farthest_neigh].similarity = dist;
		}
	}
	else{//assign free spot to new neighbor
		neighArray[i].neighbour_ID = neighbor->item_index;
		neighArray[i].similarity = euclidean_distance(dataSet[userID],dataSet[neighbor->item_index],numOfItems);
	}
	//////////////////////
}

void addToNeighMatrix2(int P, int userID, double **dataSet, neigh_struct *neighArray, int numOfItems,
		cos_vec_node *neighbor){
	////////////////////
	int i, farthest_neigh;
	double max, dist;
	////
	for(i=0;i<P;i++){
		if(neighArray[i].neighbour_ID == (-1))
			break;
		if(i == 0){
			max = neighArray[i].similarity;
			farthest_neigh = 0;
		}
		else{
			dist = neighArray[i].similarity;
			if(dist > max){
				max = dist;
				farthest_neigh = i;
			}
		}
	}
	if(i == P){//no free spot, maybe replace with farthest neighbor
		dist = cosine_dist(dataSet[userID],dataSet[neighbor->item_index],numOfItems);
		if(dist < max){
			neighArray[farthest_neigh].neighbour_ID = neighbor->item_index;
			neighArray[farthest_neigh].similarity = dist;
		}
	}
	else{//assign free spot to new neighbor
		neighArray[i].neighbour_ID = neighbor->item_index;
		neighArray[i].similarity = cosine_dist(dataSet[userID],dataSet[neighbor->item_index],numOfItems);
	}
	//////////////////////
}

void init_neighMatrix(neigh_struct ***neighbours, int P, int numOfUsers){
	int i, j;
	(*neighbours) = malloc(numOfUsers*sizeof(neigh_struct *));
	for(i=0;i<numOfUsers;i++){
		(*neighbours)[i] = malloc(P*sizeof(neigh_struct));
		for(j=0;j<P;j++)
			(*neighbours)[i][j].neighbour_ID = -1;//init neighbors IDs with (-1) aka no neighbor yet
	}
}


void free_neighMatrix(neigh_struct **neighbours, int numOfUsers){
	int i;
	for(i=0;i<numOfUsers;i++)
		free(neighbours[i]);
	free(neighbours);
}

double calculateZ(neigh_struct *neighbours, int P){
	int i;
	double sum;
	sum = 0.0;
	for(i=0;i<P;i++){
		if(neighbours[i].neighbour_ID != (-1))
			sum += neighbours[i].similarity;
	}
	if(sum == 0.0){
		printf("---ERROR---\n");
		return 0.0;
	}
	return (1.0/sum);
}

double calculateZ2(METRIC_SPACES metricSpace, int my_clust, int userID, int numOfItems, int numOfUsers, int *pointsClust, double **dataSet){
	int i;
	double sum;
	sum = 0.0;
	for(i=0;i<numOfUsers;i++){
		if(i == userID)
			continue;
		if(((pointsClust[i] == (-1)) && (i == my_clust)) || ((pointsClust[i]!=(-1)) && (pointsClust[i]==my_clust))){
			if(metricSpace == VECTOR_SPACE)
				sum += euclidean_distance(dataSet[userID],dataSet[i],numOfItems);
			else
				sum += cosine_dist(dataSet[userID],dataSet[i],numOfItems);
		}
	}
	return (1.0/sum);
}




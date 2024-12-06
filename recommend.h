#ifndef __RECOMMEND__
#define __RECOMMEND__

#include "vector.h"
#include "enums.h"
#include "multi.h"
#include "cosine.h"

typedef struct neighbours_struct{
	int neighbour_ID;
	double similarity;//similarity of neighbor with corresponding user
}neigh_struct;


///////////////////////
//clustering
void cluster_recommend(METRIC_SPACES metricSpace, int numOfUsers, int numOfItems,
		double **dataSet, int **five_best, double *avgRating, int *pointsClust);
///////////////////////

void NN_lsh_recommend(METRIC_SPACES metricSpace, int L, multiHashes *mHs, multiTables *mTs, int numOfUsers,
		int numOfItems, double **dataSet, int P, int **five_best, double *avgRating);


//stores P closest neighbors of user with userID in "neighbours" array
void binary_repeated_search(int L, multiHashes *mHs, multiTables *mTs, int numOfUsers,
		int numOfItems, double **dataSet, neigh_struct **neighbours, int userID, int P);

//binary2 only for cosine
void binary_repeated_search2(int L, multiHashes *mHs, multiTables *mTs, int numOfUsers,
		int numOfItems, double **dataSet, neigh_struct **neighbours, int userID, int P);
/////
void addToNeighMatrix(int P, int userID, double **dataSet, neigh_struct *neighArray, int numOfItems,
		vec_node *neighbor);
//add2 for cosine
void addToNeighMatrix2(int P, int userID, double **dataSet, neigh_struct *neighArray, int numOfItems,
		cos_vec_node *neighbor);

void addToTop5Matrix(int userID, double **dataSet, int *five_best, int itemID, double rating);

void init_neighMatrix(neigh_struct ***neighbours, int P, int numOfUsers);
void free_neighMatrix(neigh_struct **neighbours, int numOfUsers);
double calculateZ(neigh_struct *neighbours, int P);//returns normalizing factor 'z' of a user
double calculateZ2(METRIC_SPACES metricSpace, int my_clust, int userID, int numOfItems, int numOfUsers, int *pointsClust, double **dataSet);
/////////
#endif

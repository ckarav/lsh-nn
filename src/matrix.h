#ifndef _MATRIX_
#define _MATRIX_

#include <stdio.h>

/************** MATRIX ****************/


typedef struct h_x1x2t1t2{
    int x1, x2; //x1, x2 representing the itemID of corresponding points
    double t1, t2; //tresholds
    double d_x1x2;//distance(x1,x2)
}h_x1x2t1t2;

typedef struct matrix_hashF{
    h_x1x2t1t2 *Hi; // 'k' H functions that define the main hash_function
    int k;
}matrix_hashF;

/** BUCKET LIST **/

typedef struct matrix_node{
    int itemID;
    struct matrix_node *next;
}matrix_node;

typedef struct matrix_list{
    matrix_node *head;
    int N;
}matrix_list;

typedef struct matrix_table{
    matrix_list *buckets;
    int size;
}matrix_table;
////////////////////////////////////////////////
/*matrix hash functions*/
void init_matHash(matrix_hashF *MH, int k, int n, const double **m_array); //init matrix hash function
int hash_item(matrix_hashF *MH, double *item_dist); //returns result of hash function (after binary to decimal conversion)
void destroy_matrixHF(matrix_hashF *MH);
int isDuplicateMat(matrix_hashF *H, int n, matrix_hashF *G); /*checks that G is differnt from each of the 'n' funcs in the H array*/

double calculate_median(int x1, int x2, int n, const double **d);//calculates and returns "t1"
////////////////////////////
/*matrix tables functions*/
void matList_init(matrix_list *MB);
void matList_add(matrix_list *MB, int index);//add new node at head of list
void matList_destroy(matrix_list *MB);
int matList_size(matrix_list *MB);
void matList_clear(matrix_node *cur);
double minMatDist(matrix_list *MB, double *item_dist, matrix_node **min_node, int L);
void matrixRange(matrix_list *MB, double *item_dist, double R, matrix_list *results);
///// for clustering /////
void matrixRange2(matrix_list *MB, double *item_dist, double R, matrix_list *results, int *state_array, int iter);
//////////////////////////

void mat_table_init(int k, matrix_table *MT); //create matrix table with 2^k cells
void mat_table_destroy(matrix_table *MT);
void mat_table_insert(matrix_table *MT, matrix_hashF *G, int item_index, double *item_dist);
//////////////////////////////
/****************************/
/******** ALGORITHMS ********/
void Matrix_Search(int L, int k, const double **dists, int n, double R, double **queries, int Qnum, const char **itemNames, char **QNames, FILE *out);//main matrix_space function
int approx_NN_matrix(int L, matrix_table *tables, matrix_hashF *G, double *query, double *aDist);//returns ID of NN
void range_matrix(int L, matrix_table *tables, matrix_hashF *G, double R, double *query, matrix_list *results);//stores neighbors within R, in results list
/////////////////////////////////////////////////////////////////
/** functions for clustering **/
void InitMatrix(int L, int k, matrix_table **Tables, matrix_hashF **Gs, const double **dists, int n);
void FreeMatrix(int L, matrix_table *Tables, matrix_hashF *Gs);

void range_matrix2(int L, matrix_table *tables, matrix_hashF *G, double R, double *query, matrix_list *results, int *state_array, int iter);
//////////////////////////////////////////

/****************************************/

#endif // _MATRIX_

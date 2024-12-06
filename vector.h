#ifndef _VECTOR_
#define _VECTOR_

#include <stdio.h>

/**************** VECTOR *****************/

typedef struct h_func{
    double *vec; //d-dimensional vector with coords i.i.d. by the standard normal
    double t;
    double w;
}h_func;

typedef struct vector_hashFun{
    h_func *Hi; // 'k' functions that define the hashing function
    int *Ri; //pointer to k-array with Ris allocated in heap
    int k;
    int d;
}vector_hashFun;

/** BUCKET LIST **/

typedef struct vec_node{
    int item_index;
    long long objID; //remember object ID to maybe skip calculating the distance (long/unsigned int PENDING)
    double *vec; /*stored d-dimensional vector array*/
    struct vec_node *next; //next node
}vec_node;

typedef struct vec_list{
    vec_node *head;
    int N; //number of items in list
}vec_list;

typedef struct vec_table{
    vec_list *buckets;
    int size; //#cells of table
}vec_table;

/////////////////////////////////////////////
/*vector hash functions*/
void init_vecHash(vector_hashFun *VH, int d, int k, int w, int *Ri); /*initializes vector hash function structure*/
int hash_vec(vector_hashFun *VH, double *input, long long *objID, int table_size); /*hashes input vector, returns corresponding value and stores object ID in objID*/
void destroy_vectorHF(vector_hashFun *VH); /*free vector hash function resources*/
int isDuplicateVec(vector_hashFun *H, int n, vector_hashFun *G); /*checks that G is differnt from each of the 'n' funcs in the H array*/

int difVectors(double *vA, double *vB, int d);//returns 1 if vA != vB
////////////////////////////////
/*vector tables functions*/
void vecList_init(vec_list *VB); /*initializes list*/
void vecList_add(vec_list *VB, double *v, int index, long long obj); /* adds entry (at Head of list) to list with vector "v" */
void vecList_destroy(vec_list *VB); /*destroy entire list*/
int vecList_size(vec_list *VB); /*returns list's size*/
void vecList_clear(vec_node *cur);
void vecList_clear2(vec_list *VB, vec_node *cur);
void clearVecList(vec_list *VB);
double minVecDist(vec_list *VB, double *query, long long QobjID, int d, vec_node **min_node, int L); /*traverses list, returns min_dist and stores pointer of respective item to min_node*/
void vecRange(vec_list *VB, double *query, long long QobjID, int d, double R, vec_list *results); /*traverse list and adds to results all entries within range R of query*/
void vecRangeYolo(vec_list *VB, double *query, long long QobjID, int d, double R, vec_list *results, int yolo); /*traverse list and adds to results all entries within range R of query*/
///// for clustering //////
void vecRange2(vec_list *VB, double *query, int QobjID, int d, double R, vec_list *results, int *state_array, int iter);
///////////////////////////

void vec_table_init(int n, vec_table *VT); /*create vector table with n cells*/
void vec_table_destroy(vec_table *VT); /*ANNIHILATES  a vector table*/
void vec_table_insert(vec_table *VT, vector_hashFun *G, double *item, int item_id); /*inserts new entry in table*/
///////////////////////
/********************************/
/********** ALGORITHMS **********/
void Vector_Search(int L, int k, const double **vectors, int colNum, int rowNum, double R, double **queries, int Qnum, const char **itemNames, char **QNames, FILE *out); /*Main Vector Search function*/
int approx_NN_vector(int L, vec_table *tables, vector_hashFun *Gs, double *query, int d, double *aDist); /*approximate NN in vector space*/
void range_vector(int L, vec_table *tables, vector_hashFun *Gs, double R, double *query, int d, vec_list *results); /*search for neighbors within range R*/
void range_vector_yolo(int L, vec_table *tables, vector_hashFun *Gs, double R, double *query, int d, vec_list *results, int yolo);
/////////////////////////////////////////////////////////////////
/** functions for clustering **/
void InitVec(int L, int k, vec_table **Tables, vector_hashFun **Gs, const double **vectors, int colNum, int rowNum);
void FreeVec(int L, vec_table *Tables, vector_hashFun *Gs);

void range_vector2(int L, vec_table *tables, vector_hashFun *Gs, double R, double *query, int d, vec_list *results, int *state_array, int iter);
/////////////////////////////////////

/*****************************************/
#endif // _VECTOR_

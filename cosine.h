#ifndef _COSINE_
#define _COSINE_

typedef struct hi_func{
	double **ri;
	int k;
	int d;
}cos_hashFun;

/** BUCKET LIST **/

typedef struct cos_vec_node{
    int item_index;
    double *vec; /*stored d-dimensional vector array*/
    struct cos_vec_node *next; //next node
}cos_vec_node;

typedef struct cos_vec_list{
    cos_vec_node *head;
    int N; //number of items in list
}cos_vec_list;

typedef struct cos_vec_table{
    cos_vec_list *buckets;
    int size; //#cells of table
}cos_vec_table;

///////////////////////////////////////////////
void init_cos_vecHash(cos_hashFun *CH, int d, int k); /*initializes cosine hash function structure*/
int hash_cos_vec(cos_hashFun *CH, double *input); /*hashes input vector and returns corresponding value*/
void destroy_cos_vectorHF(cos_hashFun *CH); /*free cosine hash function resources*/
/////////////////
void cosList_init(cos_vec_list *CB); /*initializes list*/
void cosList_add(cos_vec_list *CB, double *v, int index); /* adds entry (at Head of list) to list with vector "v" */
void cosList_destroy(cos_vec_list *CB); /*destroy entire list*/
int cosList_size(cos_vec_list *CB); /*returns list's size*/
void cosList_clear2(cos_vec_list *CB, cos_vec_node *cur);
void clearCosList(cos_vec_list *CB);
void cosRange(cos_vec_list *CB, double *query, int d, double R, cos_vec_list *results); /*traverse list and adds to results all entries within range R of query*/
///////////////////
void cos_vec_table_init(int k, cos_vec_table *CT); /*create cosine table with n cells*/
void cos_vec_table_destroy(cos_vec_table *CT); /*ANNIHILATES  a cosine table*/
void cos_vec_table_insert(cos_vec_table *CT, cos_hashFun *G, double *item, int item_id); /*inserts new entry in table*/

//////////////////////////////////
void cos_vector(int L, cos_vec_table *tables, cos_hashFun *Gs, double R, double *query, int d, cos_vec_list *results); /*search for neighbors within range R*/

void InitCosine(int L, int k, cos_vec_table **Tables, cos_hashFun **Gs,
		const double **ratings, int numOfUsers, int numOfItems);
void FreeCosine(int L, cos_vec_table *Tables, cos_hashFun *Gs);

/////////////////////////////////

#endif

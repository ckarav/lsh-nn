#ifndef _HAM_
#define _HAM_

#include <stdio.h>

/************ HAMMING *************/

typedef struct ham_hash{
    int L, k; /*parameters*/
    int *g; /*hash function g is an array of k integers(in range [1..d]), representing the index of each h*/
}ham_hash;

/*** BUCKET LIST ***/

typedef struct ham_node{
    int item_index;
    char *item; /*stored string from input*/
    struct ham_node *next; /*next node in list*/
}ham_node;

typedef struct ham_list{
    ham_node *head;
    int N; /*number of items in list*/
}ham_list;

typedef struct ham_table{
    ham_list *buckets;
    int size; // #cells of table
}ham_table;

///////////////////////////////
/*hamming hash functions*/
void init_hamHash(ham_hash *HH, int L_val, int k_val, int d); /*initializes hamming hash structure*/
int hash_string(ham_hash *HH, const char *input); /*takes the input strings and returns result of hash function (in decimal)*/
void destroy_ham(ham_hash *HH); /*frees resources of ham_hash struct*/
int isDuplicateHam(ham_hash *H, int n, int k, ham_hash *G); /*checks that G is differnt from each of the 'n' funcs in the H array*/

//////////////////////////////
/*hamming tables functions*/
void hamList_init(ham_list *HB); /*initializes list*/
void hamList_add(ham_list *HB, char *s, int index); /* adds entry (at Head of list) to list with string "s" */
void hamList_destroy(ham_list *HB); /*destroy entire list*/
int hamList_size(ham_list *HB); /*returns list's size*/
void hamList_clear(ham_node *cur);/*removes similar nodes after "cur" node*/
int minHamDist(ham_list *HB, char *query, ham_node **min_node, int L); /*traverses list, returns min_dist and stores pointer of respective item to min_node*/
void hamRange(ham_list *HB, char *query, double R, ham_list *results); /*traverse list and adds to results all entries within range R of query*/
///// for clustering //////
void hamRange2(ham_list *HB, char *query, double R, ham_list *results, int *state_array, int iter);
//////////////////////////

void ham_table_init(int k, ham_table *HT); /*create hamming table with 2^k cells*/
void ham_table_destroy(ham_table *HT); /*can you guess what this function does?*/
void ham_table_insert(ham_table *HT, ham_hash *G, const char *item, int item_id); /*inserts new entry in table*/
//////////////////////////////
/**********************************/
/********** ALGORITHMS ************/
void Hamming_Search(int L, int k, const char **bitstrings, int colNum, int rowNum, double R, char **queries, int Qnum, const char **itemNames, char **QNames, FILE *out); /*Main Hamming function*/
int approx_NN_ham(int L, ham_table *tables, ham_hash *G, char *query, int *aDist); /*executes approximate NN algorithm in hamming space (returns -1 if nothing found)*/
void range_ham(int L, ham_table *tables, ham_hash *G, double R, char *query, ham_list *results); /*search for neighbors in range R and stores them in results*/

/////////////////////////////////////////////////////////////////
/** functions for clustering **/
void InitHam(int L, int k, ham_table **Tables, ham_hash **Gs, const char **bitstrings, int colNum, int rowNum);
void FreeHam(int L, ham_hash *Gs, ham_table *Tables);

void range_ham2(int L, ham_table *tables, ham_hash *G, double R, char *query, ham_list *results, int *state_array, int iter);
//////////////////////////////////////////////////////////
/**********************************/

#endif

#include "cosine.h"
#include "distsNmore.h"
#include <stdlib.h>
#include <stdio.h>
#include "rands.h"
#include <math.h>


void init_cos_vecHash(cos_hashFun *CH, int d, int k){
	int i;
	CH->d = d;
	CH->k = k;
	CH->ri = malloc(k*sizeof(double *));
	for(i=0;i<k;i++){
		(CH->ri)[i] = malloc(d*sizeof(double));
		gauss_rand2((CH->ri)[i],d);
	}
}

int hash_cos_vec(cos_hashFun *CH, double *input){
	char *hash_output;
	int i, dec_val;
	//////
	hash_output = malloc((CH->k)*sizeof(char));
	////
	for(i=0;i<(CH->k);i++){
		if(inner_product((CH->ri)[i],input,CH->d) >= 0.0)
			hash_output[i] = '1';
		else
			hash_output[i] = '0';
	}
	dec_val = bin_to_dec(hash_output,CH->k);
	////////
	free(hash_output);
	return dec_val;
}

void destroy_cos_vectorHF(cos_hashFun *CH){
	int i;
	for(i=0;i<(CH->k);i++)
		free((CH->ri)[i]);
	free(CH->ri);
}

/////////////////

void cosList_init(cos_vec_list *CB){
	CB->head = NULL;
	CB->N = 0;
}

void cosList_add(cos_vec_list *CB, double *v, int index){
	cos_vec_node *temp;
	temp = malloc(sizeof(cos_vec_node));
	temp->item_index = index;
	temp->vec = v;
	temp->next = CB->head;
	CB->head = temp;
	CB->N++;
}

void cosList_destroy(cos_vec_list *CB){
	cos_vec_node *temp, *temp2;
	temp = CB->head;
	while(temp){
		temp2 = temp->next;
		free(temp);
		temp = temp2;
	}
}

int cosList_size(cos_vec_list *CB){
	return (CB->N);
}

void cosList_clear2(cos_vec_list *CB, cos_vec_node *cur){
	cos_vec_node *temp, *temp2;
	temp = cur;
	while(temp->next){
		if((temp->next->item_index) == (cur->item_index)){
			temp2 = temp->next;
			temp->next = temp->next->next;
			free(temp2);
			CB->N--;
		}
		else{
			temp = temp->next;
		}
	}
}

void clearCosList(cos_vec_list *CB){
	cos_vec_node *temp;
	temp = CB->head;
	while(temp){
		cosList_clear2(CB,temp);
		temp = temp->next;
	}
}

void cosRange(cos_vec_list *CB, double *query, int d, double R, cos_vec_list *results){
	cos_vec_node *temp;
	/////////////
	temp = CB->head;
	while(temp){
		if(cosine_dist(query,temp->vec,d) < R)
			cosList_add(results,temp->vec,temp->item_index);
		temp = temp->next;
	}
}

///////////////////

void cos_vec_table_init(int k, cos_vec_table *CT){
	int i;
	CT->size = 1 << k;//table size is 2^k
	CT->buckets = malloc((CT->size)*sizeof(cos_vec_list));
	for(i=0;i<(CT->size);i++)
		cosList_init(&((CT->buckets)[i]));
}

void cos_vec_table_destroy(cos_vec_table *CT){
	int i;
	for(i=0;i<(CT->size);i++)
		cosList_destroy(&((CT->buckets)[i]));
	free(CT->buckets);
}

void cos_vec_table_insert(cos_vec_table *CT, cos_hashFun *G, double *item, int item_id){
	int bucket_num;
	bucket_num = hash_cos_vec(G,item);
	cosList_add(&((CT->buckets)[bucket_num]),item,item_id);
}

///////////////////////////////////////
void cos_vector(int L, cos_vec_table *tables, cos_hashFun *Gs, double R, double *query,
		int d, cos_vec_list *results){
	//////////////////
	int i, hash_val;
	////
	for(i=0;i<L;i++){
	    hash_val = hash_cos_vec(&(Gs[i]),query);
	    cosRange(&(((tables[i]).buckets)[hash_val]),query,d,R,results);
	}
}

//////////////////////

void InitCosine(int L, int k, cos_vec_table **Tables, cos_hashFun **Gs,
		const double **ratings, int numOfUsers, int numOfItems){
	///////////
	int i, j;
	//////////////
	*Gs = malloc(L*sizeof(cos_hashFun));
	*Tables = malloc(L*sizeof(cos_vec_table));
	for(i=0;i<L;i++){
		init_cos_vecHash(&((*Gs)[i]),numOfItems,k);
		cos_vec_table_init(k,&((*Tables)[i]));
		for(j=0;j<numOfUsers;j++)
			cos_vec_table_insert(&((*Tables)[i]),&((*Gs)[i]),ratings[j],j);
	}
}



void FreeCosine(int L, cos_vec_table *Tables, cos_hashFun *Gs){
	int i;
	for(i=0;i<L;i++){
		destroy_cos_vectorHF(&(Gs[i]));
		cos_vec_table_destroy(&(Tables[i]));
	}
	free(Gs);
	free(Tables);
}













#include "hamming.h"
#include "distsNmore.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rands.h"
#include <unistd.h>


/******************* HAMMING ********************/
void init_hamHash(ham_hash *HH, int L_val, int k_val, int d){
    int i;
    /////////////
    HH->L = L_val;
    HH->k = k_val;
    HH->g = malloc(k_val * sizeof(int)); /*allocate space for hash function*/
    for(i=0;i<k_val;i++)
        (HH->g)[i] = uniform_randInt(1,d); /*pick random "h" functions*/
}
/////////////////

int hash_string(ham_hash *HH, const char *input){
    char *hash_output;
    int i, dec_val;
    ////////////////////
    hash_output = malloc((HH->k) * sizeof(char)); /*allocate space for the k output bits*/
    ///////////////
    for(i=0;i<(HH->k);i++) /*calculate hash ouput*/
        hash_output[i] = input[(HH->g)[i] - 1];
    dec_val = bin_to_dec(hash_output,HH->k); /*calculate decimal value of binary string*/
    /////////////////////////////
    free(hash_output); //free temp string
    return dec_val;
}
//////////////

void destroy_ham(ham_hash *HH){
    free(HH->g);
}

int isDuplicateHam(ham_hash *H, int n, int k, ham_hash *G){
    int i, j, duplicate;
    for(i=0;i<n;i++){
        duplicate = 1; //init as duplicate
        for(j=0;j<k;j++){
            if((G->g[j])!=((H[i].g[j]))){
                duplicate = 0;
                break;
            }
        }
        if(duplicate)
            return 1;
    }
    return 0; //there is not duplicate
}

/////////////
/***********************************************/

void hamList_init(ham_list *HB){
    HB->head = NULL;
    HB->N = 0;
}
////////////

void hamList_add(ham_list *HB, char *s, int index){
    ham_node *temp;
    /*add new entry at head of bucket*/
    temp = malloc(sizeof(ham_node));
    temp->item = s;
    temp->item_index = index;
    temp->next = HB->head;
    HB->head = temp;
    HB->N++; /*increase number of entries in bucket*/
}
/////////

void hamList_destroy(ham_list *HB){
    ham_node *temp, *temp2;
    temp = HB->head;
    while(temp){
        temp2 = temp->next;
        free(temp);
        temp = temp2;
    }
}

int hamList_size(ham_list *HB){
    return (HB->N);
}

void hamList_clear(ham_node *cur){
    ham_node *temp, *temp2;
    temp = cur;
    while(temp->next){
        if((temp->next->item_index) == (cur->item_index)){
            temp2 = temp->next;
            temp->next = temp->next->next;
            free(temp2);
        }
        else{
            temp = temp->next;
        }
    }
}

////////

int minHamDist(ham_list *HB, char *query, ham_node **min_node,int L){
    ham_node *temp;
    int min_dist, t, count;
    min_dist = -1; //init with -1 to encode "infinity"
    (*min_node) = NULL;
    temp = HB->head;
    count = 0;//num of nodes we have checked, if >3L, stop
    while(temp){
        if(count > 3*L)//if we checked too many nodes, return best thus far
            return min_dist;
        if(min_dist == -1){
            min_dist = hamming_distance(query,temp->item,strlen(query));
            (*min_node) = temp;
            count++;
        }
        else{
            t = hamming_distance(query,temp->item,strlen(query));
            if(t<min_dist){
                min_dist = t;
                (*min_node) = temp;
            }
            count++;
        }
        temp = temp->next;
    }
    return min_dist;
}

void hamRange(ham_list *HB, char *query, double R, ham_list *results){
    ham_node *temp;
    ///////////////
    temp = HB->head;
    while(temp){
        if((double)(hamming_distance(query,temp->item,strlen(query))) < R)
            hamList_add(results,temp->item,temp->item_index);
        temp = temp->next;
    }
}


////////

void ham_table_init(int k, ham_table *HT){
    int i;
    HT->size = 1 << k; /*table size is 2^k*/
    HT->buckets = malloc((HT->size) * sizeof(ham_list));
    for(i=0;i<(HT->size);i++)
        hamList_init(&((HT->buckets)[i]));
}

void ham_table_destroy(ham_table *HT){
    int i;
    for(i=0;i<(HT->size);i++)
        hamList_destroy(&((HT->buckets)[i]));
    free(HT->buckets);
}

void ham_table_insert(ham_table *HT, ham_hash *G, const char *item, int item_id){
    int bucket_num;
    bucket_num = hash_string(G,item);
    hamList_add(&((HT->buckets)[bucket_num]),item,item_id);
}

/*********************************************/

void Hamming_Search(int L, int k, const char **bitstrings, int colNum, int rowNum, double R, char **queries, int Qnum, const char **itemNames, char **QNames, FILE *out){
    int i, j, aNN, appDist, appTime, exhDist, exhTime;
    ham_hash *Gs;
    ham_table *Tables;
    ham_list range_results;
    ////////
    ham_node  *temp;
    ///////
    /********* BUILD HASH-TABLES ***********/
    Gs = malloc(L * sizeof(ham_hash));
    Tables = malloc(L * sizeof(ham_table));
    for(i=0;i<L;i++){ //initialize stuff
        init_hamHash(&(Gs[i]),L,k,colNum);
        //check for same hash functions and if so, regenerate
        if(i>0){
            while(isDuplicateHam(Gs,i,k,&(Gs[i]))){
                destroy_ham(&(Gs[i]));
                init_hamHash(&(Gs[i]),L,k,colNum);
            }
        }
        //////////////////////////////////////////////
        ham_table_init(k,&(Tables[i])); /* create table of size 2^k */
        /** store input points **/
        for(j=0;j<rowNum;j++)
            ham_table_insert(&(Tables[i]),&(Gs[i]),bitstrings[j],j+1);//j+1 item ID
    }
    /**************************************/
    for(i=0;i<Qnum;i++){ //for each query
        fprintf(out,"Query: %s\n",QNames[i]);
        if(R){ //case of range search (R > 0)
            hamList_init(&range_results); //we have already defined a list for results to use in case of Range
            range_ham(L,Tables,Gs,R,queries[i],&range_results);
            ////////////////////////////////////////////////////////////
            /*writing results to file*/
            temp = range_results.head;
            fprintf(out,"R-near neighbors:\n");
            while(temp){
                fprintf(out,"%s\n",itemNames[(temp->item_index)-1]);
                hamList_clear(temp); //remove same node from rest of list
                temp = temp->next;
            }
            //////////////////////////////////////////////////////////////
            hamList_destroy(&range_results); //destroy results list
        }
        else
            fprintf(out,"R-near neighbors: -\n");
        appTime = GetClockTimeInMilliSec(); //TIME_START
        aNN = approx_NN_ham(L,Tables,Gs,queries[i],&appDist);
        appTime = GetClockTimeInMilliSec() - appTime; //TIME_END
        if(aNN != -1)
            fprintf(out,"Nearest neighbor: %s\ndistanceLSH: %d\n",itemNames[aNN-1],appDist);
        else
            fprintf(out,"Nearest neighbor: not found\ndistanceLSH: +oo\n");
        //}
        exhaustiveSearchHammingSpace(bitstrings,rowNum,colNum,(const char*)queries[i],&exhDist,&exhTime);
        fprintf(out,"distanceTrue: %d\ntLSH: ",exhDist);
        PrintTime(out,appTime);
        fprintf(out,"\ntTrue: ");
        PrintTime(out,exhTime);fprintf(out,"\n\n");
        ///////////////////////////////////
    }
    /*****************************************/
    /************ FREE RESOURCES *************/
    for(i=0;i<L;i++){
        destroy_ham(&(Gs[i]));
        ham_table_destroy(&(Tables[i]));
    }
    free(Gs); //free hash functions
    free(Tables);
}




int approx_NN_ham(int L, ham_table *tables, ham_hash *G, char *query, int *aDist){
    int min_dist, temp_dist, i, sizeee;
    ham_node *min_node, *temp_min_node; /*will store pointer to node with min distance*/

    min_dist = -1;
    for(i=0;i<L;i++){ /*search corresponding table*/
        if(min_dist == -1)
            min_dist = minHamDist(&(((tables[i]).buckets)[hash_string(&(G[i]),query)]),query,&min_node,L);
        else{
            temp_dist = minHamDist(&(((tables[i]).buckets)[hash_string(&(G[i]),query)]),query,&temp_min_node,L);
            if((temp_dist != -1) && (temp_dist < min_dist)){
                min_dist = temp_dist;
                min_node = temp_min_node;
            }
        }
    }
    (*aDist) = min_dist;
    if(min_dist == -1) //if no neighbor found, return -1
        return -1;
    return (min_node->item_index); //else return index of nearest neighbor's ID
}




void range_ham(int L, ham_table *tables, ham_hash *G, double R, char *query, ham_list *results){
    int i;
    ////////////
    for(i=0;i<L;i++)
        hamRange(&(((tables[i]).buckets)[hash_string(&(G[i]),query)]),query,R,results);
}


/*********************************************/
/************************************************/
/*********** clustering functions ***************/
void InitHam(int L, int k, ham_table **Tables, ham_hash **Gs, const char **bitstrings, int colNum, int rowNum){
    int i, j;
    /////////////////////////////////////////////
    /********* BUILD HASH-TABLES ***********/
    *Gs = malloc(L * sizeof(ham_hash));
    *Tables = malloc(L * sizeof(ham_table));
    for(i=0;i<L;i++){ //initialize stuff
        init_hamHash(&((*Gs)[i]),L,k,colNum);
        //check for same hash functions and if so, regenerate
        if(i>0){
            while(isDuplicateHam((*Gs),i,k,&((*Gs)[i]))){
                destroy_ham(&((*Gs)[i]));
                init_hamHash(&((*Gs)[i]),L,k,colNum);
            }
        }
        //////////////////////////////////////////////
        ham_table_init(k,&((*Tables)[i])); /* create table of size 2^k */
        /** store input points **/
        for(j=0;j<rowNum;j++)
            ham_table_insert(&((*Tables)[i]),&((*Gs)[i]),bitstrings[j],j+1);//j+1 item ID
    }
}

void FreeHam(int L, ham_hash *Gs, ham_table *Tables){
    int i;
    /************ FREE RESOURCES *************/
    for(i=0;i<L;i++){
        destroy_ham(&(Gs[i]));
        ham_table_destroy(&(Tables[i]));
    }
    free(Gs); //free hash functions
    free(Tables);
}

void hamRange2(ham_list *HB, char *query, double R, ham_list *results, int *state_array, int iter){
    ham_node *temp;
    int state;
    ///////////////
    temp = HB->head;
    while(temp){
        state = state_array[temp->item_index];
        if((state == -1) || ((state >= 0) && (state < iter))){//if centroid or already assigned in previous range search, skip
            temp = temp->next;
            continue;
        }
        if((double)(hamming_distance(query,temp->item,strlen(query))) < R)
            hamList_add(results,temp->item,temp->item_index);
        temp = temp->next;
    }
}

void range_ham2(int L, ham_table *tables, ham_hash *G, double R, char *query, ham_list *results, int *state_array, int iter){
    int i;
    ////////////
    for(i=0;i<L;i++)
        hamRange2(&(((tables[i]).buckets)[hash_string(&(G[i]),query)]),query,R,results,state_array,iter);
}

/************************************************/

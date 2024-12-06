#include "vector.h"
#include "distsNmore.h"
#include <stdlib.h>
#include <stdio.h>
#include "rands.h"
#include <math.h>


void init_vecHash(vector_hashFun *VH, int d, int k, int w, int *Ri){
    int i;
    ////////
    VH->d = d;
    VH->k = k;
    VH->Ri = Ri;//Ris allocated and initialized in main vector function
    VH->Hi = malloc(k*sizeof(h_func)); //allocate space for 'h' functions
    for(i=0;i<k;i++){
        ((VH->Hi)[i]).w = w;
        ((VH->Hi)[i]).t = uniform_randR(w);//t is uniform in [0,w)
        ((VH->Hi)[i]).vec = malloc(d*sizeof(double));
        gauss_rand(((VH->Hi)[i]).vec,d); //randomize vector
    }
}

int hash_vec(vector_hashFun *VH, double *input, long long *objID, int table_size){
    int i;
    unsigned long M = (1ULL << 32) - 5;
    (*objID) = 0;
    for(i=0;i<(VH->k);i++)
        (*objID) += ((long long)((VH->Ri)[i]) * (long long)floorl((inner_product(input,((VH->Hi)[i]).vec,VH->d) + (((VH->Hi)[i]).t))/4.0)) % M;
    return (int)((*objID) % ((long long)table_size));
}

void destroy_vectorHF(vector_hashFun *VH){
    int i;
    for(i=0;i<(VH->k);i++)
        free(((VH->Hi)[i]).vec);
    free(VH->Hi);
}

int isDuplicateVec(vector_hashFun *H, int n, vector_hashFun *G){
    int i, j, duplicate;
    for(i=0;i<n;i++){
        duplicate = 1;
        for(j=0;j<(G->k);j++){
            if( (G->Hi[j].t != H[i].Hi[j].t) || (difVectors(G->Hi[j].vec,H[i].Hi[j].vec,G->d)) ){
                duplicate = 0;
                break;
            }
        }
        if(duplicate)
            return 1;
    }
    return 0;
}

int difVectors(double *vA, double *vB, int d){
    int i;
    for(i=0;i<d;i++){
        if(vA[i] != vB[i])
            return 1;
    }
    return 0;
}
//////////////////////////////////
/**************************************************/

void vecList_init(vec_list *VB){
    VB->head = NULL;
    VB->N = 0;
}

void vecList_add(vec_list *VB, double *v, int index, long long obj){
    vec_node *temp;
    temp = malloc(sizeof(vec_node));
    temp->item_index = index;
    temp->objID = obj;
    temp->vec = v;
    temp->next = VB->head;
    VB->head = temp;
    VB->N++;
}

void vecList_destroy(vec_list *VB){
    vec_node *temp, *temp2;
    temp = VB->head;
    while(temp){
        temp2 = temp->next;
        free(temp);
        temp = temp2;
    }
}

int vecList_size(vec_list *VB){
    return (VB->N);
}

void vecList_clear(vec_node *cur){
    vec_node *temp, *temp2;
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

void vecList_clear2(vec_list *VB, vec_node *cur){
	vec_node *temp, *temp2;
	temp = cur;
	while(temp->next){
	    if((temp->next->item_index) == (cur->item_index)){
	        temp2 = temp->next;
	        temp->next = temp->next->next;
	        free(temp2);
	        VB->N--;//reduce size of list
	    }
	    else{
	        temp = temp->next;
	    }
	}
}

void clearVecList(vec_list *VB){
	vec_node *temp;
	temp = VB->head;
	while(temp){
		vecList_clear2(VB,temp);
		temp = temp->next;
	}
}

double minVecDist(vec_list *VB, double *query, long long QobjID, int d, vec_node **min_node, int L){
    int count;
    vec_node *temp;
    double min_dist, temp_dist;
    min_dist = -1.0; //init with -1 to encode "infinity"
    (*min_node) = NULL;
    temp = VB->head;
    count = 0;
    while(temp){
        if(count > 3*L)
            return min_dist;
        if(temp->objID != QobjID){ //if they have different object IDs, don't compare
            temp = temp->next;
            continue;
        }
        if(min_dist == -1.0){
            min_dist = euclidean_distance(query,temp->vec,d);
            (*min_node) = temp;
            count++;
        }
        else{
            temp_dist = euclidean_distance(query,temp->vec,d);
            if(temp_dist < min_dist){
                min_dist = temp_dist;
                (*min_node) = temp;
            }
            count++;
        }
        temp = temp->next;
    }
    return min_dist;
}

void vecRange(vec_list *VB, double *query, long long QobjID, int d, double R, vec_list *results){
	vecRangeYolo(VB, query, QobjID, d, R, results, 0);
}

void vecRangeYolo(vec_list *VB, double *query, long long QobjID, int d, double R, vec_list *results, int yolo){
    vec_node *temp;
    //////////////
    temp = VB->head;
    while(temp){
        if (yolo && temp->objID != QobjID){
            temp = temp->next;
            continue;
        }
        if(euclidean_distance(query,temp->vec,d) < R)
            vecList_add(results,temp->vec,temp->item_index,temp->objID);
        temp = temp->next;
    }
}

//////////////
void vec_table_init(int n, vec_table *VT){
    int i;
    VT->size = n; //'n' cells ('n' is not necessarily size of input(but w/e is passed as argument...))
    VT->buckets = malloc(n*sizeof(vec_list));
    for(i=0;i<n;i++)
        vecList_init(&((VT->buckets)[i]));
}

void vec_table_destroy(vec_table *VT){
    int i;
    for(i=0;i<(VT->size);i++)
        vecList_destroy(&((VT->buckets)[i]));
    free(VT->buckets);
}

void vec_table_insert(vec_table *VT, vector_hashFun *G, double *item, int item_id){
    int bucket_num;
    long long objID;
    bucket_num = hash_vec(G,item,&objID,VT->size);
    vecList_add(&((VT->buckets)[bucket_num]),item,item_id,objID);
}

/**************************************************/

void Vector_Search(int L, int k, const double **vectors, int colNum, int rowNum, double R, double **queries, int Qnum, const char **itemNames, char **QNames, FILE *out){
    int i, j, aNN, *Ri, appTime, exhTime;
    double appDist, exhDist;
    vector_hashFun *Gs;
    vec_table *Tables;
    vec_list range_results;
    /////////////////////////////
    vec_node *temp;
    /////////////////////////////
    /******** BUILD HASH-TABLES ********/
    Gs = malloc(L*sizeof(vector_hashFun));
    Tables = malloc(L*sizeof(vec_table));
    ////////////////
    Ri = malloc(k*sizeof(int)); //generate Ris, same for all hash functions
    for(i=0;i<k;i++)
        Ri[i] = rand();
    ///////////////
    for(i=0;i<L;i++){
        init_vecHash(&(Gs[i]),colNum,k,4,Ri); //w = 4
        //check for same vecHashes, if so regenerate
        if(i>0){
            while(isDuplicateVec(Gs,i,&(Gs[i]))){
                destroy_vectorHF(&(Gs[i]));
                init_vecHash(&(Gs[i]),colNum,k,4,Ri);
            }
        }
        //////////////////////////
        vec_table_init(rowNum/2,&(Tables[i]));//create table of size n/2
        /** store input points **/
        for(j=0;j<rowNum;j++)
            vec_table_insert(&(Tables[i]),&(Gs[i]),vectors[j],j+1);
    }
    /****************************************/
    for(i=0;i<Qnum;i++){
        fprintf(out,"Query: %s\n",QNames[i]);
        if(R){
            vecList_init(&range_results);//init list used to store results
            range_vector(L,Tables,Gs,R,queries[i],colNum,&range_results);
            /////////////////
            /*writing results to file*/
            temp = range_results.head;
            fprintf(out,"R-near neighbors:\n");
            while(temp){ //print range results
                fprintf(out,"%s\n",itemNames[(temp->item_index)-1]);
                vecList_clear(temp);//clear rest of list from duplicates
                temp = temp->next;
            }
            /////////////////
            vecList_destroy(&range_results); //destroy results list
        }
        else
            fprintf(out,"R-near neighbors: -\n");
        appTime = GetClockTimeInMilliSec(); //time_start
        aNN = approx_NN_vector(L,Tables,Gs,queries[i],colNum,&appDist);
        appTime = GetClockTimeInMilliSec() - appTime; //time_end
        //////////////////////
        if(aNN != -1.0)
            fprintf(out,"Nearest neighbor: %s\ndistanceLSH: %f\n",itemNames[aNN-1],sqrt(appDist));
        else
            fprintf(out,"Nearest neighbor: not found\ndistanceLSH: +oo\n");
        exhaustiveSearchVectorSpace(vectors,rowNum,colNum,(const double*)queries[i],&exhDist,&exhTime);
        fprintf(out,"distanceTrue: %f\ntLSH: ",exhDist);
        PrintTime(out,appTime);
        fprintf(out,"\ntTrue: ");
        PrintTime(out,exhTime);fprintf(out,"\n\n");
    }
    /****************************************/
    /************ FREE RESOURCES ************/
    for(i=0;i<L;i++){
        destroy_vectorHF(&(Gs[i]));
        vec_table_destroy(&(Tables[i]));
    }
    free(Gs); //free hash functions
    free(Tables); //free hash tables
    free(Ri); //free Ris
}


int approx_NN_vector(int L, vec_table *tables, vector_hashFun *Gs, double *query, int d, double *aDist){
    int i, hash_val;
    long long obj;
    double min_dist, temp_dist;
    vec_node *min_node, *temp_min_node;

    min_dist = -1.0;
    for(i=0;i<L;i++){
        hash_val = hash_vec(&(Gs[i]),query,&obj,tables[i].size); //calculate hash_value and object ID
        if(min_dist == -1)
            min_dist = minVecDist(&(((tables[i]).buckets)[hash_val]),query,obj,d,&min_node,L);
        else{
            temp_dist = minVecDist(&(((tables[i]).buckets)[hash_val]),query,obj,d,&temp_min_node,L);
            if((temp_dist != -1.0) && (temp_dist < min_dist)){
                min_dist = temp_dist;
                min_node = temp_min_node;
            }
        }
    }
    (*aDist) = min_dist;
    if(min_dist == -1.0)
        return -1;
    else
        return (min_node->item_index);
}



void range_vector(int L, vec_table *tables, vector_hashFun *Gs, double R, double *query, int d, vec_list *results){
	range_vector_yolo(L, tables, Gs, R, query, d, results, 0);
}

void range_vector_yolo(int L, vec_table *tables, vector_hashFun *Gs, double R, double *query, int d, vec_list *results, int yolo) {
    int i, hash_val;
    long long obj;
    for(i=0;i<L;i++){
        hash_val = hash_vec(&(Gs[i]),query,&obj,tables[i].size);
        vecRangeYolo(&(((tables[i]).buckets)[hash_val]),query,obj,d,R,results, yolo);
    }
}


/************************************************/
/*********** clustering functions ***************/
void InitVec(int L, int k, vec_table **Tables, vector_hashFun **Gs, const double **vectors, int colNum, int rowNum){
    int i, j, *Ri;
    ////////////////////////////
    /******** BUILD HASH-TABLES ********/
    *Gs = malloc(L*sizeof(vector_hashFun));
    *Tables = malloc(L*sizeof(vec_table));
    ////////////////
    Ri = malloc(k*sizeof(int)); //generate Ris, same for all hash functions
    for(i=0;i<k;i++)
        Ri[i] = rand();
    ///////////////
    for(i=0;i<L;i++){
    	init_vecHash(&((*Gs)[i]),colNum,k,4,Ri); //w = 4
        //init_vecHash(&(Gs[i]),colNum,k,4,Ri); //w = 4
        //check for same vecHashes, if so regenerate
        if(i>0){
            while(isDuplicateVec(*Gs,i,&((*Gs)[i]))){
                destroy_vectorHF(&((*Gs)[i]));
                init_vecHash(&((*Gs)[i]),colNum,k,4,Ri);
            }
        }
        //////////////////////////
        vec_table_init(rowNum/8,&((*Tables)[i]));//create table of size n/8
        /** store input points **/
        for(j=0;j<rowNum;j++)
            vec_table_insert(&((*Tables)[i]),&((*Gs)[i]),vectors[j],j);
    }
}

void FreeVec(int L, vec_table *Tables, vector_hashFun *Gs){
    int i;
	/************ FREE RESOURCES ************/
    free((Gs[0]).Ri);//PENDING CHECK
    for(i=0;i<L;i++){
        destroy_vectorHF(&(Gs[i]));
        vec_table_destroy(&(Tables[i]));
    }
    free(Gs); //free hash functions
    free(Tables); //free hash tables
}
void vecRange2(vec_list *VB, double *query, int QobjID, int d, double R, vec_list *results, int *state_array, int iter){
    vec_node *temp;
    int state;
    //////////////
    temp = VB->head;
    while(temp){
        state = state_array[temp->item_index];
        if(temp->objID != QobjID){
            temp = temp->next;
            continue;
        }
        if((state == -1) || ((state >= 0) && (state < iter))){//if centroid or already assigned in previous range search, skip
            temp = temp->next;
            continue;
        }
        if(euclidean_distance(query,temp->vec,d) < R)
            vecList_add(results,temp->vec,temp->item_index,temp->objID);
        temp = temp->next;
    }
}

void range_vector2(int L, vec_table *tables, vector_hashFun *Gs, double R, double *query, int d, vec_list *results, int *state_array, int iter){
    int i, hash_val;
    long long obj;
    for(i=0;i<L;i++){
        hash_val = hash_vec(&(Gs[i]),query,&obj,tables[i].size);
        vecRange2(&(((tables[i]).buckets)[hash_val]),query,obj,d,R,results,state_array,iter);
    }
}
/**************************************************/



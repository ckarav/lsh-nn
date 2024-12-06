#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include "rands.h"
#include "distsNmore.h"

/***************** MATRIX *******************/

void init_matHash(matrix_hashF *MH, int k, int n, const double **m_array){
    int i;
    MH->k = k;
    MH->Hi = malloc(k*sizeof(h_x1x2t1t2));//allocate space for 'k' h functions
    for(i=0;i<k;i++){
        MH->Hi[i].x1 = uniform_randInt(1,n);//x1 uniformly distributed in [1,n]
        MH->Hi[i].x2 = uniform_randInt(1,n);//x2 like x1...
        while((MH->Hi[i].x1)==(MH->Hi[i].x2))//make sure not same x1 and x2
            MH->Hi[i].x2 = uniform_randInt(1,n);
        /////////////////
        MH->Hi[i].t1 = calculate_median(MH->Hi[i].x1,MH->Hi[i].x2,n,m_array);//calculate t1 value
        MH->Hi[i].t2 = -1; //encode "infinity"...
        MH->Hi[i].d_x1x2 = m_array[(MH->Hi[i].x1)-1][(MH->Hi[i].x2)-1];
    }
}

int hash_item(matrix_hashF *MH, double *item_dist){
    char *hash_output;
    int i, dec_val, x1, x2;
    double hx1x2, dx1x2;
    ////////////
    hash_output = malloc((MH->k)*sizeof(char));//allocate space for the 'k' output bits
    ///////////
    for(i=0;i<(MH->k);i++){
        x1 = MH->Hi[i].x1;
        x2 = MH->Hi[i].x2;
        dx1x2 = MH->Hi[i].d_x1x2;
        hx1x2 = ((item_dist[x1-1]*item_dist[x1-1]) + (item_dist[x2-1]*item_dist[x2-1]) - (dx1x2*dx1x2))/(2*dx1x2);
        if(hx1x2 >= ((double)(MH->Hi[i].t1)))
            hash_output[i] = '1';
        else
            hash_output[i] = '0';
    }
    dec_val = bin_to_dec(hash_output,MH->k);
    //////////
    free(hash_output);
    return dec_val;
}


void destroy_matrixHF(matrix_hashF *MH){
    free(MH->Hi);
}

int isDuplicateMat(matrix_hashF *H, int n, matrix_hashF *G){
    int i, j, duplicate;
    for(i=0;i<n;i++){
        duplicate = 1;
        for(j=0;j<(G->k);j++){
            if((G->Hi[j].x1 != H[i].Hi[j].x1) || (G->Hi[j].x2 != H[i].Hi[j].x2) || (G->Hi[j].t1 != H[i].Hi[j].t1)){
                duplicate = 0;
                break;
            }
        }
        if(duplicate)
            return 1;
    }
    return 0;
}

double calculate_median(int x1, int x2, int n, const double **d){
    int i;
    double median;
    double *hx1x2; //array that holds set of number we want the median of
    hx1x2 = malloc(n*sizeof(double));
    for(i=0;i<n;i++)
        hx1x2[i] = ((d[i][x1-1]*d[i][x1-1])+(d[i][x2-1]*d[i][x2-1])-(d[x1-1][x2-1]*d[x1-1][x2-1]))/(2*d[x1-1][x2-1]);
    quicksort(n,hx1x2);/*sort array*/
    if(n%2) //if odd number of points
        median = hx1x2[n/2]; //median is middle element
    else //if even number of points
        median = (hx1x2[n/2 - 1] + hx1x2[n/2])/2.0; //median is the mean of the 2 medians
    /////////////
    free(hx1x2);
    return median;
}
//////////////////
/*********************************************/
void matList_init(matrix_list *MB){
    MB->head = NULL;
    MB->N = 0;
}

void matList_add(matrix_list *MB, int index){
    matrix_node *temp;
    temp = malloc(sizeof(matrix_node));
    temp->itemID = index;
    temp->next = MB->head;
    MB->head = temp;
    MB->N++;
}

void matList_destroy(matrix_list *MB){
    matrix_node *temp, *temp2;
    temp = MB->head;
    while(temp){
        temp2 = temp->next;
        free(temp);
        temp = temp2;
    }
}

int matList_size(matrix_list *MB){
    return (MB->N);
}

void matList_clear(matrix_node *cur){
    matrix_node *temp, *temp2;
    temp = cur;
    while(temp->next){
        if((temp->next->itemID) == (cur->itemID)){
            temp2 = temp->next;
            temp->next = temp->next->next;
            free(temp2);
        }
        else{
            temp = temp->next;
        }
    }
}

///////

double minMatDist(matrix_list *MB, double *item_dist, matrix_node **min_node, int L){
    matrix_node *temp;
    double min_dist, t;
    int count;
    min_dist = -1;
    (*min_node) = NULL;
    temp = MB->head;
    count = 0;
    while(temp){
        if(count > 3*L) //if we checked too many nodes, return best thus far
            return min_dist;
        if(min_dist == -1){
            min_dist = item_dist[(temp->itemID)-1];
            (*min_node) = temp;
            count++;
        }
        else{
            t = item_dist[(temp->itemID)-1];
            if(t<min_dist){
                min_dist = t;
                (*min_node) = temp;
            }
            count++;
        }
        temp = temp->next;
    }
    //printf("min_dist is %d\n",min_dist);
    return min_dist;
}
////////

void matrixRange(matrix_list *MB, double *item_dist, double R, matrix_list *results){
    matrix_node *temp;
    temp = MB->head;
    while(temp){
        if(item_dist[(temp->itemID)-1] < R)
            matList_add(results,temp->itemID);
        temp = temp->next;
    }
}
///////////////

void mat_table_init(int k, matrix_table *MT){
    int i;
    MT->size = 1 << k; //table size is 2^k
    MT->buckets = malloc((MT->size)*sizeof(matrix_list));//allocate space for buckets
    for(i=0;i<(MT->size);i++)
        matList_init(&((MT->buckets)[i]));
}

void mat_table_destroy(matrix_table *MT){
    int i;
    for(i=0;i<(MT->size);i++)
        matList_destroy(&((MT->buckets)[i]));
    free(MT->buckets);
}

void mat_table_insert(matrix_table *MT, matrix_hashF *G, int item_index, double *item_dist){
    int bucket_num;
    bucket_num = hash_item(G,item_dist);
    matList_add(&((MT->buckets)[bucket_num]),item_index);
}

/***********************************************/

void Matrix_Search(int L, int k, const double **dists, int n, double R, double **queries, int Qnum, const char **itemNames, char **QNames, FILE *out){
    int i, j, aNN, appTime, exhTime, exhDist;
    double appDist;
    matrix_hashF *Gs;
    matrix_table *Tables;
    matrix_list range_results;
    matrix_node *temp;
    /////////////////////////////
    /******** BUILD HASH-TABLES ********/
    Gs = malloc(L*sizeof(matrix_hashF));
    Tables = malloc(L*sizeof(matrix_table));
    ////////////////////////////
    for(i=0;i<L;i++){
        init_matHash(&(Gs[i]),k,n,dists);
        //check for same hash functions and if so, regenerate
        if(i>0){
            while(isDuplicateMat(Gs,i,&(Gs[i]))){
                destroy_matrixHF(&(Gs[i]));
                init_matHash(&(Gs[i]),k,n,dists);
            }
        }
        ///////////////////////////////
        mat_table_init(k,&(Tables[i]));//create table of size 2^k
        for(j=0;j<n;j++)
            mat_table_insert(&(Tables[i]),&(Gs[i]),j+1,dists[j]);
    }
    /***********************************/
    for(i=0;i<Qnum;i++){//for each query
        if(R){
            matList_init(&range_results);
            range_matrix(L,Tables,Gs,R,queries[i],&range_results);
            /////////////////////////
            /*writing results to file*/
            temp = range_results.head;
            fprintf(out,"R-near neighbors:\n");
            while(temp){
                fprintf(out,"%s\n",itemNames[(temp->itemID)-1]);
                matList_clear(temp);//clear rest of list from duplicates
                temp = temp->next;
            }
            ///////////////////////
            matList_destroy(&range_results);
        }
        else
            fprintf(out,"R-near neighbors: -\n");
        appTime = GetClockTimeInMilliSec(); //time_start
        aNN = approx_NN_matrix(L,Tables,Gs,queries[i],&appDist);
        appTime = GetClockTimeInMilliSec() - appTime; //time_end
        if(aNN != -1)
            fprintf(out,"Nearest neighbor: %s\ndistanceLSH: %f\n",itemNames[aNN-1],appDist);
        else
            fprintf(out,"Nearest neighbor: not found\ndistanceLSH: +oo\n");
        exhaustiveSearchMatrixSpace(n,(const int*)queries[i],&exhDist,&exhTime);
        fprintf(out,"distanceTrue: %d\ntLSH: ",exhDist);
        PrintTime(out,appTime);
        fprintf(out,"\ntTrue: ");
        PrintTime(out,exhTime);fprintf(out,"\n\n");
    }
    /**********************************/
    /********* FREE RESOURCES *********/
    for(i=0;i<L;i++){
        destroy_matrixHF(&(Gs[i]));
        mat_table_destroy(&(Tables[i]));
    }
    free(Gs);
    free(Tables);
}


int approx_NN_matrix(int L, matrix_table *tables, matrix_hashF *G, double *query, double *aDist){
    double min_dist, temp_dist;
    int i, hash_val;
    matrix_node *min_node, *temp_min_node;

    min_dist = -1;
    for(i=0;i<L;i++){
        hash_val = hash_item(&(G[i]),query);//calculate hash_value (cell index of table[i])
        if(min_dist == -1)
            min_dist = minMatDist(&(((tables[i]).buckets)[hash_val]),query,&min_node,L);
        else{
            temp_dist = minMatDist(&(((tables[i]).buckets)[hash_val]),query,&temp_min_node,L);
            if((temp_dist != -1) && (temp_dist < min_dist)){
                min_dist = temp_dist;
                min_node = temp_min_node;
            }
        }
    }
    (*aDist) = min_dist;
    if(min_dist == -1)
        return -1;
    return (min_node->itemID);
}

void range_matrix(int L, matrix_table *tables, matrix_hashF *G, double R, double *query, matrix_list *results){
    int i, hash_val;
    for(i=0;i<L;i++){
        hash_val = hash_item(&(G[i]),query);
        matrixRange(&(((tables[i]).buckets)[hash_val]),query,R,results);
    }
}

/************************************************/
/*********** clustering functions ***************/
void InitMatrix(int L, int k, matrix_table **Tables, matrix_hashF **Gs, const double **dists, int n){
    int i, j;
    ///////////////
    /******** BUILD HASH-TABLES ********/
    *Gs = malloc(L*sizeof(matrix_hashF));
    *Tables = malloc(L*sizeof(matrix_table));
    ////////////////////////////
    for(i=0;i<L;i++){
        init_matHash(&((*Gs)[i]),k,n,dists);
        //check for same hash functions and if so, regenerate
        if(i>0){
            while(isDuplicateMat((*Gs),i,&((*Gs)[i]))){
                destroy_matrixHF(&((*Gs)[i]));
                init_matHash(&((*Gs)[i]),k,n,dists);
            }
        }
        ///////////////////////////////
        mat_table_init(k,&((*Tables)[i]));//create table of size 2^k
        for(j=0;j<n;j++)
            mat_table_insert(&((*Tables)[i]),&((*Gs)[i]),j+1,dists[j]);
    }
}

void FreeMatrix(int L, matrix_table *Tables, matrix_hashF *Gs){
    int i;
    /********* FREE RESOURCES *********/
    for(i=0;i<L;i++){
        destroy_matrixHF(&(Gs[i]));
        mat_table_destroy(&(Tables[i]));
    }
    free(Gs);
    free(Tables);
}

void matrixRange2(matrix_list *MB, double *item_dist, double R, matrix_list *results, int *state_array, int iter){
    matrix_node *temp;
    int state;
    ///////////////////
    temp = MB->head;
    while(temp){
        state = state_array[temp->itemID];
        if((state == -1) || ((state >= 0) && (state < iter))){//if centroid or already assigned in previous range search, skip
            temp = temp->next;
            continue;
        }
        if(item_dist[(temp->itemID)-1] < R)
            matList_add(results,temp->itemID);
        temp = temp->next;
    }
}


void range_matrix2(int L, matrix_table *tables, matrix_hashF *G, double R, double *query, matrix_list *results, int *state_array, int iter){
    int i, hash_val;
    for(i=0;i<L;i++){
        hash_val = hash_item(&(G[i]),query);
        matrixRange2(&(((tables[i]).buckets)[hash_val]),query,R,results,state_array,iter);
    }
}

/*********************************************/

#include "distsNmore.h"
#include <sys/time.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <pthread.h>
#include <unistd.h>
#include "rands.h"


double cosine_dist(double *vectorA, double *vectorB, int dim){
	return (1.0-(inner_product(vectorA,vectorB,dim) / (length(vectorA, dim)*length(vectorB, dim))));

}

int hamming_distance(char *a, char *b, int len){
    int i, dist;
    dist = 0;
    for(i=0;i<len;i++){
        if(a[i] != b[i])
            dist++;
    }
    return dist;
}


double euclidean_distance(double *vectorA, double *vectorB, int dim){
    int i;
    double sum = 0.0;
    for(i=0;i<dim;i++)
        sum += (vectorA[i] - vectorB[i]) * (vectorA[i] - vectorB[i]);
    //no need to calculate sqrt... (dist1 > dist2 <==> dist1^2 > dist2^2)
    return sum;
}
void transform_conformations(double **confA, int N){
	int i;
	double xc0, xc1, xc2;
	xc0=xc1=xc2 = 0;
	for(i=0;i<N;i++){ //calculate center of conformation
		xc0 += confA[i][0]/(double)N;
		xc1 += confA[i][1]/(double)N;
		xc2 += confA[i][2]/(double)N;
	}
	for(i=0;i<N;i++){ //update conformation
		confA[i][0] -= xc0;
		confA[i][1] -= xc1;
		confA[i][2] -= xc2;
	}
}

double cRMSD(gsl_matrix *X, gsl_matrix *Y, int N){
	int i, signum;
	double dist = -1.0, detQ, trace;
	gsl_matrix *transpX_mult_Y, *V, *Q, *Q_LUdecomp, *X_mult_Q;
	gsl_vector *S, *work;
	gsl_permutation *p;
	///////
	/** allocate resources for temporary structures **/
	p = gsl_permutation_alloc(3); //PENDING (size 3?)
	S = gsl_vector_alloc(3);
	work = gsl_vector_alloc(3);
	transpX_mult_Y = gsl_matrix_alloc(3,3);
	X_mult_Q = gsl_matrix_alloc(N,3);
	V = gsl_matrix_alloc(3,3);
	Q = gsl_matrix_alloc(3,3);
	Q_LUdecomp = gsl_matrix_alloc(3,3);
	///////////////////////////////////////////////////
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,(const gsl_matrix*)X,(const gsl_matrix*)Y,0.0,transpX_mult_Y);
	gsl_linalg_SV_decomp(transpX_mult_Y,V,S,work);//calculate SVD
	//at this point, transpX_mult_Y array holds 'U'
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,(const gsl_matrix*)transpX_mult_Y,(const gsl_matrix*)V,0.0,Q);
	gsl_matrix_memcpy(Q_LUdecomp,(const gsl_matrix*)Q);//copy matrix Q to Q_LUdecomp
	gsl_linalg_LU_decomp(Q_LUdecomp,p,&signum);//calculate LU decomposition of Q
	detQ = gsl_linalg_LU_det(Q_LUdecomp,signum);//calculate determinant of Q
	if(detQ < 0){
		//update Q
		double temp;
		for(i=0;i<3;i++){ //update last column of 'U' matrix(a.k.a. transpX_mult_Y)
			temp = gsl_matrix_get((const gsl_matrix*)transpX_mult_Y,i,2);
			gsl_matrix_set(transpX_mult_Y,i,2,temp*(-1.0));
		}
		gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,(const gsl_matrix*)transpX_mult_Y,(const gsl_matrix*)V,0.0,Q);
	}
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,(const gsl_matrix*)X,(const gsl_matrix*)Q,0.0,X_mult_Q);
	gsl_matrix_sub(X_mult_Q,(const gsl_matrix*)Y);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,(const gsl_matrix*)X_mult_Q,(const gsl_matrix*)X_mult_Q,0.0,Q);//use Q to hold a temporary array(actual Q no longer needed)
	trace = gsl_matrix_get(Q,0,0) + gsl_matrix_get(Q,1,1) + gsl_matrix_get(Q,2,2);//calculate trace of matrix
	trace = sqrt(trace);//get square root of trace (which is the Frobenius norm we wanted)
	dist = trace / ((double)N);
	////////////////////////////////
	/********* free resources *********/
	gsl_permutation_free(p);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_matrix_free(transpX_mult_Y);
	gsl_matrix_free(X_mult_Q);
	gsl_matrix_free(V);
	gsl_matrix_free(Q);
	gsl_matrix_free(Q_LUdecomp);
	////////////////////////////////////
	return dist;
}

multiDist general_distance(METRIC_SPACES metricSpace, int objA, int objB, const multiMatrix *data, int colNum){
    multiDist dist;
    if(metricSpace == HAMMING_SPACE)
        dist.i = hamming_distance(data->c[objA],data->c[objB],colNum);
    else if(metricSpace == VECTOR_SPACE)
        dist.d = euclidean_distance(data->d[objA],data->d[objB],colNum);
    else if(metricSpace == MATRIX_SPACE)
    	dist.d = (data->d)[objA][objB];
    else if(metricSpace == COSINE)
    	dist.d = cosine_dist(data->d[objA],data->d[objB],colNum);
    else //case of Proteins (conformations dataset)
    	dist.d = cRMSD(data->b[objA],data->b[objB],colNum);
    return dist;
}

int lessThan(METRIC_SPACES metricSpace, multiDist mD1, multiDist mD2){
    if((metricSpace == MATRIX_SPACE) || (metricSpace == VECTOR_SPACE) || (metricSpace == BIO_SPACE) || (metricSpace == COSINE)){
        if(mD1.d < mD2.d)
            return 1;
        return 0;
    }
    else{ //hamming
        if(mD1.i < mD2.i)
            return 1;
        return 0;
    }
}

int lessThan2(METRIC_SPACES metricSpace, multiLongDist mD1, multiLongDist mD2){
    if((metricSpace == MATRIX_SPACE) || (metricSpace == VECTOR_SPACE) || (metricSpace == BIO_SPACE) || (metricSpace == COSINE)){
        if(mD1.d < mD2.d)
            return 1;
        return 0;
    }
    else{ //hamming
        if(mD1.i < mD2.i)
            return 1;
        return 0;
    }
}

multiLongDist addDist(METRIC_SPACES metricSpace, multiDist mD1, multiLongDist mD2){
    multiLongDist newDist;
    if((metricSpace == MATRIX_SPACE) || (metricSpace == VECTOR_SPACE) || (metricSpace == BIO_SPACE) || (metricSpace == COSINE))
        newDist.d = (mD1.d)/10000.0 + mD2.d;
    else //hamming
        newDist.i = mD1.i + mD2.i;
    return newDist;
}

void zeroInit(METRIC_SPACES metricSpace, multiLongDist *mlD){
    if((metricSpace == MATRIX_SPACE) || (metricSpace == VECTOR_SPACE) || (metricSpace == BIO_SPACE) || (metricSpace == COSINE))
        mlD->d = 0.0;
    else //hamming
        mlD->i = 0;
}

////////////////////////////////////

double length(double *vec, int d){
	int i;
	double sum;
	sum = 0.0;
	for(i=0;i<d;i++)
		sum += vec[i] * vec[i];
	return sqrt(sum);
}

int bin_to_dec(char *s, int length){
    int i, sum;
    sum = 0;
    for(i=1;i<=length;i++)
        sum += (1 << (length-i)) * (s[i-1] - '0');
    return sum;
}



double inner_product(double *vA, double *vB, int d){
    int i;
    double sum = 0.0;
    for(i=0;i<d;i++)
        sum += vA[i] * vB[i];
    return sum;
}

//////////////////////////////////////////////
/**
Code for Quicksort from Itroduction to Programming slides
**/

void swapd(double *a, double *b)       /* Just exchange two doubles */
{ int temp;
  temp = *a;
  *a = *b;
  *b = temp; }


void quicksort_body(double *x, int up, int down)
{ int start, end;
  start = up;              /* Save start position of small elements */
  end = down;                /* Save end position of large elements */
  while (up < down) {            /* Pivot element is at up position */
    while (x[down] >= x[up] && up < down)      /* Let down elements */
      down--;              /* larger than pivot stay where they are */
    if (up != down) {                    /* If pivot is not reached */
      swapd(&x[up], &x[down]);   /* echange it with smaller element */
      up++;     /* Pivot is at down position, move up a bit further */
    }
    while (x[up] <= x[down] && up < down)        /* Let up elements */
    up++;                 /* smaller than pivot stay where they are */
    if (up != down) {                    /* If pivot is not reached */
      swapd(&x[up], &x[down]);   /* exchange it with larger element */
      down--;   /* Pivot is at up position, move down a bit further */
    } }       /* Now up = down is the position of the pivot element */
  if (start < up-1) /* Is there at least one element left of pivot? */
    quicksort_body(x, start, up-1); /* Quick(sort) smaller elements */
  if (end > down+1)/* Is there at least one element right of pivot? */
    quicksort_body(x, down+1, end);  /* Quick(sort) larger elements */
}

void quicksort(int n, double *x)
{ quicksort_body(x, 0, n-1);    /* Call recursive quicksort to sort */
}          /* elements of the array from position 0 to position n-1 */


void* safeMalloc(size_t sizeOfElement) {
	void* space = NULL;
	space = malloc(sizeOfElement);
	if (space == NULL) {
		printf("Out of memory error!\n");
		exit(-1);
	}
	return space;
}

///////////////////////////////
int GetClockTimeInMilliSec()
{
	struct timeval t2; gettimeofday(&t2, NULL);
	return t2.tv_sec * 1000 + t2.tv_usec / 1000;
}

void PrintTime(FILE* output, int milli_sec)
{
	int v = milli_sec;
	int hours = v / (1000 * 60 * 60); v %= (1000 * 60 * 60);
	int minutes = v / (1000 * 60); v %= (1000 * 60);
	int seconds = v / 1000; v %= 1000;
	int milli_seconds = v;
	int first = 1;
	if (hours) { if (!first) fprintf(output, ":"); fprintf(output, "%dh", hours); first = 0; }
	if (minutes) { if (!first) fprintf(output, ":"); fprintf(output, "%dm", minutes); first = 0; }
	if (seconds) { if (!first) fprintf(output, ":"); fprintf(output, "%ds", seconds); first = 0; }
	if (milli_seconds) { if (!first) fprintf(output, ":"); fprintf(output, "%dms", milli_seconds); first = 0; }
	else fprintf(output,"%dms",milli_seconds);
}

int exhaustiveSearchVectorSpace(const double** vectors, int rowNum, int colNum, const double* p, double* distanceTrue, int* tTrue) {
	int i, j;
	double dist;
	double minDist = 0;
	int minIndex = 0;
	int timestamp = GetClockTimeInMilliSec();

	for (i = 0; i < colNum; i++) minDist += pow(vectors[0][i] - p[i], 2.0);

	for (i = 1; i < rowNum; i++) {
		dist = 0;
		for (j = 0; j < colNum; j++) dist += pow(vectors[i][j] - p[j], 2.0);

		if (dist < minDist) {
			minDist = dist;
			minIndex = i;
		}
	}

	*tTrue = GetClockTimeInMilliSec() - timestamp;
	if (minDist != 0)
		*distanceTrue = sqrt(minDist);
	else
		*distanceTrue = 0;

	return minIndex;
}

int exhaustiveSearchHammingSpace(const char** bitstrings, int rowNum, int colNum, const char* p, int* distanceTrue, int* tTrue) {
	int i, j;
	double dist;
	double minDist = 0;
	int minIndex = 0;
	int timestamp = GetClockTimeInMilliSec();

	for (i = 0; i < colNum; i++) {
		if (bitstrings[0][i] != p[i]) minDist++;
	}

	for (i = 1; i < rowNum; i++) {
		dist = 0;
		for (j = 0; j < colNum; j++) {
			if (bitstrings[i][j] != p[j]) dist++;
		}

		if (dist < minDist) {
			minDist = dist;
			minIndex = i;
		}
	}

	*tTrue = GetClockTimeInMilliSec() - timestamp;
	*distanceTrue = minDist;
	return minIndex;
}

int exhaustiveSearchMatrixSpace(int colNum, const int* p, int* distanceTrue, int* tTrue) {
	int i;
	int minDist = p[0];
	int minIndex = 0;
	int timestamp = GetClockTimeInMilliSec();

	for (i = 1; i < colNum; i++) {
		if (p[i] < minDist) {
			minDist = p[i];
			minIndex = i;
		}
	}

	*tTrue = GetClockTimeInMilliSec() - timestamp;
	*distanceTrue = minDist;
	return minIndex;
}

// Shared variables between threads.
pthread_mutex_t counterMtx = PTHREAD_MUTEX_INITIALIZER;
int counter;
double** dists;
double* avgRatingsRef;

typedef struct DataStruct {
	double*** conformations;
	multiMatrix* datasetMatrix;
	int numConform;
	IntPair* indexPairs;
	int r;

	double** recommendations;
	int userNum;
	int itemNum;
} DataStruct;
typedef DataStruct* DataPtr;

void* dRMSD_job1(void* dataPtr) {
	int i = 0, j, index;
	DataPtr ptr = (DataPtr)dataPtr;

	// For each conformation, calculate the distance between the points of the (index)-th pair.
	while (counter < ptr->numConform) {
		pthread_mutex_lock(&counterMtx); // Lock the mutex before trying to access the counter.
		index = counter++; // Ensure that this thread's index is unique by getting the next available index.
		pthread_mutex_unlock(&counterMtx);
		if (index >= ptr->numConform) break;

		for (i = 0; i < ptr->r; i++) {
			dists[index][i] = 0;
			for (j = 0; j < 3; j++) {
				double tmp = ptr->conformations[index][ptr->indexPairs[index].l][j] - ptr->conformations[index][ptr->indexPairs[index].r][j];
				dists[index][i] += tmp * tmp;
			}
			dists[index][i] = sqrt(dists[index][i]);
		}
	}
	return NULL;
}

void* dRMSD_job2(void* dataPtr) {
	int col = 0, j, row;
	DataPtr ptr = (DataPtr)dataPtr;

	// For each conformation, calculate the distance between the points of the (col)-th pair.
	while (counter < ptr->numConform) {
		pthread_mutex_lock(&counterMtx); // Lock the mutex before trying to access the counter.
		row = counter++; // Ensure that this thread's row is unique by getting the next available row.
		pthread_mutex_unlock(&counterMtx);
		if (row >= ptr->numConform) break;

		for (col = 0; col < ptr->numConform; col++) {
			(ptr->datasetMatrix->d)[row][col] = 0;
			if (row == col) continue;

			for (j = 0; j < ptr->r; j++) {
				double tmp = dists[row][j] - dists[col][j];
				(ptr->datasetMatrix->d)[row][col] += tmp * tmp;
			}
			(ptr->datasetMatrix->d)[row][col] = (ptr->datasetMatrix->d)[row][col] / (double)ptr->r;
		}
	}
	return NULL;
}

void dRMSD(double*** conformations, multiMatrix* datasetMatrix, int numConform, int N, int r) {
	int i = 0, j, rnd1, rnd2, cores = sysconf(_SC_NPROCESSORS_ONLN);
	IntPair* indexPairs = safeMalloc(numConform * sizeof(IntPair));
	pthread_t* thread_id = safeMalloc(cores * sizeof(pthread_t));
	DataStruct data;
	dists = safeMalloc(numConform * sizeof(double*));
	for (j = 0; j < numConform; j++) dists[j] = safeMalloc(r * sizeof(double));

	// Pre-calculate the indexes of each pair of points.
	while (i < numConform) {
		do {
			rnd1 = uniform_randInt(0, N - 1);
			rnd2 = uniform_randInt(0, N - 1);
		} while (rnd1 == rnd2);

		// Ensure that the same pair of points does not already exist in the indexPairs array.
		for (j = 0; j < i; j++) {
			if (((indexPairs[j].l == rnd1) && (indexPairs[j].r == rnd2)) || ((indexPairs[j].l == rnd2) && (indexPairs[j].r == rnd1))) continue;
		}
		indexPairs[i].l = rnd1;
		indexPairs[i].r = rnd2;
		i++;
	}

	// Prepare the data for the threads.
	data.conformations = conformations;
	data.datasetMatrix = datasetMatrix;
	data.numConform = numConform;
	data.indexPairs = indexPairs;
	data.r = r;

	// Start the threads and wait for all of them to finish the distance calculation for each conformation's pair of points stored in indexPairs.
	counter = 0;
	for (i = 0; i < cores; i++) pthread_create(&thread_id[i], NULL, dRMSD_job1, &data);
	for (i = 0; i < cores; i++) pthread_join(thread_id[i], NULL);

	// Start the threads and wait for all of them to finish the distance calculation between each conformation pair.
	counter = 0; // Reset the counter.
	for (i = 0; i < cores; i++) pthread_create(&thread_id[i], NULL, dRMSD_job2, &data);
	for (i = 0; i < cores; i++) pthread_join(thread_id[i], NULL);

	// Deallocate the memory.
	for (i = 0; i < numConform; i++) free(dists[i]);
	free(dists);
	free(indexPairs);
}

void* avgR_job(void* dataPtr) {
	int userID, itemID, ratedItemsNum;
	DataPtr ptr = (DataPtr)dataPtr;

	// For each user, calculate the average rating of the items he rated.
	while (counter < ptr->userNum) {
		pthread_mutex_lock(&counterMtx); // Lock the mutex before trying to access the counter.
		userID = counter++; // Ensure that this thread's row is unique by getting the next available row.
		pthread_mutex_unlock(&counterMtx);
		if (userID >= ptr->userNum) break;

		ratedItemsNum = 0;
		avgRatingsRef[userID] = 0;
		for (itemID = 0; itemID < ptr->itemNum; itemID++) {
			if (ptr->recommendations[userID][itemID] != 0) {
				avgRatingsRef[userID] += ptr->recommendations[userID][itemID];
				ratedItemsNum++;
			}
		}
		avgRatingsRef[userID] = avgRatingsRef[userID] / (double)ratedItemsNum;
	}
	return NULL;
}

void calculateAvgRating(double** recommendations, double *avgR, int userNum, int itemNum) {
	DataStruct data;
	int i, cores = sysconf(_SC_NPROCESSORS_ONLN);
	pthread_t* thread_id = safeMalloc(cores * sizeof(pthread_t));

	// Prepare the data for the threads.
	data.recommendations = recommendations;
	data.userNum = userNum;
	data.itemNum = itemNum;

	// Start the threads and wait for all of them to finish the distance calculation for each conformation's pair of points stored in indexPairs.
	counter = 0;
	avgRatingsRef = avgR;
	for (i = 0; i < cores; i++) pthread_create(&thread_id[i], NULL, avgR_job, &data);
	for (i = 0; i < cores; i++) pthread_join(thread_id[i], NULL);
}

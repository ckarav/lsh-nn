#ifndef _DISTS_
#define _DISTS_

#include <stdio.h>
#include <stdlib.h> // For size_t.
#include "multi.h"
#include "enums.h"
#include <gsl/gsl_matrix.h>

typedef struct IntPair {
	int l; // left
	int r; // right
} IntPair;


double cosine_dist(double *vectorA, double *vectorB, int dim);

int hamming_distance(char *a, char *b, int len);

double euclidean_distance(double *vectorA, double *vectorB, int dim);

void transform_conformations(double **confA, int N);
double cRMSD(gsl_matrix *X, gsl_matrix *Y, int N);

void dRMSD(double*** conformations, multiMatrix* datasetMatrix, int numConform, int N, int r);

void calculateAvgRating(double** recommendations, double *avgR, int userNum, int itemNum);

multiDist general_distance(METRIC_SPACES matricSpace, int objA, int objB, const multiMatrix *data, int colNum);

int lessThan(METRIC_SPACES metricSpace, multiDist mD1, multiDist mD2); //returns 1 if mD1 < mD2
int lessThan2(METRIC_SPACES metricSpace, multiLongDist mlD1, multiLongDist mlD2);
multiLongDist addDist(METRIC_SPACES metricSpace, multiDist mD1, multiLongDist mD2); //returns mD1+mD2
void zeroInit(METRIC_SPACES metricSpace, multiLongDist *mlD); //sets the right mlD var to 0

////////////////////////////////////////////////
/*HELPING FUNCTIONS*/

double length(double *vec, int d);//returns of vector in d-dimensional space
int bin_to_dec(char *s, int length); /*takes string of 0's and 1's and returns value in decimal*/
double inner_product(double *vA, double *vB, int d); /*returns the i_p of 2 d-dimensional vectors*/


void quicksort(int n, double *x);
void quicksort_body(double *x, int up, int down);
void swapd(double *a, double *b);

// Safely allocates (sizeOfElement) bytes of space and terminates the program if the space was not allocated.
void* safeMalloc(size_t sizeOfElement);

int GetClockTimeInMilliSec();

// Print the time from milliseconds to the HH:MM:SS:MS format.
void PrintTime(FILE* output, int milli_sec);

/*
Searches every vector in the vectors array to find the one closest to p.
The parameter p should be a vector of colNum-dimension.
The parameter distanceTrue will be set to the value of the minimum distance after the exhaustive search is completed.
The parameter tTrue will be set to the amount of time it took to complete the exhaustive search.
*/
int exhaustiveSearchVectorSpace(const double** vectors, int rowNum, int colNum, const double* p, double* distanceTrue, int* tTrue);

/*
Searches every bitstring in the bitstrings array to find the one closest to p.
The parameter p should be a bitstring of (colNum) length.
The parameter distanceTrue will be set to the value of the minimum distance after the exhaustive search is completed.
The parameter tTrue will be set to the amount of time it took to complete the exhaustive search.
*/
int exhaustiveSearchHammingSpace(const char** bitstrings, int rowNum, int colNum, const char* p, int* distanceTrue, int* tTrue);

/*
Searches every distance in the matrix to find the item closest to p.
The parameter distanceTrue will be set to the value of the minimum distance after the exhaustive search is completed.
The parameter tTrue will be set to the amount of time it took to complete the exhaustive search.
*/
int exhaustiveSearchMatrixSpace(int colNum, const int* p, int* distanceTrue, int* tTrue);

#endif // _DISTS_

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "init.h"
#include "rands.h"
#include "distsNmore.h"

void findMinVectorSpace(const double** vectors, int* centroidIndexes, double* D, int rowNum, int colNum, int clusterNum) {
	int i, j, k;
	double dist;

	for (i = 0; i < rowNum; i++) {
		dist = 0;
		D[i] = 0;
		for (j = 0; j < colNum; j++) D[i] += (vectors[i][j] - vectors[centroidIndexes[0]][j]) * (vectors[i][j] - vectors[centroidIndexes[0]][j]);

		for (k = 1; k < clusterNum; k++) {
			for (j = 0; j < colNum; j++)
				dist += (vectors[i][j] - vectors[centroidIndexes[k]][j]) * (vectors[i][j] - vectors[centroidIndexes[k]][j]);

			if (dist < D[i]) D[i] = dist;
		}
		if (D[i] != 0) D[i] = sqrt(D[i]);
	}
}

void findMinBioSpace(const multiMatrix* gslMatrix, int* centroidIndexes, double* D, int rowNum, int colNum, int clusterNum) {
	int i, k;
	double dist;

	for (i = 0; i < rowNum; i++) {
		dist = 0;
		D[i] = general_distance(BIO_SPACE, i, centroidIndexes[0], gslMatrix, colNum).d;

		for (k = 1; k < clusterNum; k++) {
			dist += general_distance(BIO_SPACE, i, centroidIndexes[k], gslMatrix, colNum).d;
			if (dist < D[i]) D[i] = dist;
		}
		if (D[i] != 0) D[i] = sqrt(D[i]);
	}
}

void findMinHammingSpace(const char** bitstrings, int* centroidIndexes, double* D, int rowNum, int colNum, int clusterNum) {
	int i, j, k;
	double dist;

	for (i = 0; i < rowNum; i++) {
		dist = 0;
		D[i] = 0;
		for (j = 0; j < colNum; j++) {
			if (bitstrings[i][j] != bitstrings[centroidIndexes[0]][j]) D[i]++;
		}

		for (k = 1; k < clusterNum; k++) {
			for (j = 0; j < colNum; j++) {
				if (bitstrings[i][j] != bitstrings[centroidIndexes[k]][j]) dist++;
			}

			if (dist < D[i]) D[i] = dist;
		}
	}
}

void findMinMatrixSpace(const double** matrix, int* centroidIndexes, double* D, int colNum, int clusterNum) {
	int i, k;

	for (i = 0; i < colNum; i++) {
		D[i] = matrix[centroidIndexes[0]][i];

		for (k = 1; k < clusterNum; k++) {
			if (matrix[centroidIndexes[k]][i] < D[i]) {
				D[i] = matrix[centroidIndexes[k]][i];
			}
		}
	}
}

int binSearch(double* array, int arrayLength, double target) {
	int low = 0;
	int mid = -1;
	int high = arrayLength - 1;

	while (low <= high) {
		mid = (low + high) / 2;

		if (array[mid] < target)
			low = mid + 1;
		else if (array[mid] == target)
			return mid;
		else
			high = mid - 1;
	}

	while (target > array[mid] && mid < arrayLength - 1) mid++; // Find the closest element to the target.
	return mid;
}

void calculateVectorDistanceMatrix(const double** datasetMatrix, int rowNum, int colNum, double** distMatrix) {
	int i, j, k;

	for (i = 0; i < rowNum; i++) {
		for (j = i; j < rowNum; j++) {
			distMatrix[i][j - i] = 0;

			// If the element is on the main diagonal, its distance (from itself) is zero.
			if (i == j) continue;

			// The distance between the vectors i and j is calculated.
			for (k = 0; k < colNum; k++) {
				distMatrix[i][j - i] += (datasetMatrix[i][k] - datasetMatrix[j][k]) * (datasetMatrix[i][k] - datasetMatrix[j][k]);
			}
			if (distMatrix[i][j - i] != 0) distMatrix[i][j - i] = sqrt(distMatrix[i][j - i]); // The distance is the square root of the calculated value.
		}
	}
}

void calculateHammingDistanceMatrix(const char** datasetMatrix, int rowNum, int colNum, int** distMatrix) {
	int i, j, k;

	for (i = 0; i < rowNum; i++) {
		for (j = i; j < colNum; j++) {
			// If the element is on the main diagonal, its distance (from itself) is zero.
			if (i == j) {
				distMatrix[i][j - i] = 0;
				continue;
			}

			// The distance between the bitstrings i and j is calculated.
			distMatrix[i][j - i] = 0;
			for (k = 0; k < colNum; k++) {
				if (datasetMatrix[i][k] != datasetMatrix[j][k]) distMatrix[i][j - i]++;
			}
		}
	}
}

void calculateConcentrationArray(METRIC_SPACES metricSpace, const multiDistMatrix* distMatrix, int rowNum, double* v) {
	int i, j, t;
	double innerDistSum;

	for (i = 0; i < rowNum; i++) {
		v[i] = 0;
		for (j = 0; j < rowNum; j++) {
			innerDistSum = 0;
			for (t = 0; t < rowNum; t++) {
				if (metricSpace != HAMMING_SPACE) {
					if (j <= t)
						innerDistSum += distMatrix->d[j][t];
					else
						innerDistSum += distMatrix->d[t][j];
				} else {
					if (j <= t)
						innerDistSum += distMatrix->i[j][t];
					else
						innerDistSum += distMatrix->i[t][j];
				}
			}

			if (innerDistSum != 0) {
				if (metricSpace != HAMMING_SPACE) {
					if (i <= j)
						v[i] += distMatrix->d[i][j] / innerDistSum;
					else
						v[i] += distMatrix->d[j][i] / innerDistSum;
				} else {
					if (i <= j)
						v[i] += distMatrix->i[i][j] / innerDistSum;
					else
						v[i] += distMatrix->i[j][i] / innerDistSum;
				}
			}
		}
	}
}

void k_medoids_init(METRIC_SPACES metricSpace, const multiMatrix* datasetMatrix, int rowNum, int colNum, int clusterNum, int* centroidIndexes) {
	double* D;
	double* partialSums;
	double max;
	int i, k;

	D = safeMalloc(rowNum * sizeof(double));
	partialSums = safeMalloc(rowNum * sizeof(double));

	centroidIndexes[0] = uniform_randInt(0, rowNum);
	for (k = 1; k < clusterNum; k++) {
		memset(partialSums, 0, rowNum); // Reset the partial sums.

		if (metricSpace == VECTOR_SPACE) {
			findMinVectorSpace((const double**)datasetMatrix->d, centroidIndexes, D, rowNum, colNum, k);
		} else if (metricSpace == HAMMING_SPACE) {
			findMinHammingSpace((const char**)datasetMatrix->c, centroidIndexes, D, rowNum, colNum, k);
		} else if (metricSpace == MATRIX_SPACE) {
			findMinMatrixSpace((const double**)datasetMatrix->d, centroidIndexes, D, colNum, k);
		} else if (metricSpace == BIO_SPACE) {
			findMinBioSpace(datasetMatrix, centroidIndexes, D, rowNum, colNum, k);
		}

		// Find the max element.
		max = D[0];
		for (i = 1; i < rowNum; i++) {
			if (D[i] > max) max = D[i];
		}

		// Calculate the normalized partial sums (partial sums divided by the maximum element of the D array.
		partialSums[0] += (D[0] * D[0]) / (max * max);
		for (i = 1; i < rowNum; i++) partialSums[i] += ((D[i] * D[i]) / (max * max)) + partialSums[i - 1];
		centroidIndexes[k] = binSearch(partialSums, rowNum, uniformRandDouble(partialSums[rowNum - 1]));
	}

	free(D);
	free(partialSums);
}

void concentrate_init(METRIC_SPACES metricSpace, const multiMatrix* datasetMatrix, int rowNum, int colNum, int clusterNum, int* centroidIndexes) {
	multiDistMatrix distMatrix;
	multiDistMatrix distMatrixReference; // This is used to reference the datasetMatrix in case the metric space is the MATRIX_SPACE.
	double* v; // The concentration values are stored here.
	int i, j;
	double curMin;

	// Allocate space for the triangular matrix and the concentration array.
	v = safeMalloc(rowNum * sizeof(double));
	if (metricSpace == VECTOR_SPACE) {
		distMatrix.d = safeMalloc(rowNum * sizeof(double*));
		for (i = 0; i < rowNum; i++) distMatrix.d[i] = safeMalloc((rowNum - i) * sizeof(double));
	} else if (metricSpace == HAMMING_SPACE) {
		distMatrix.i = safeMalloc(rowNum * sizeof(int*));
		for (i = 0; i < rowNum; i++) distMatrix.i[i] = safeMalloc((rowNum - i) * sizeof(int));
	}

	/*
	 * Create the distance matrix (only if the metric space is not the MATRIX_SPACE, because then it is already prepared).
	 * After the distance matrix is created the concentration array values are calculated.
	 */
	if (metricSpace == VECTOR_SPACE) {
		calculateVectorDistanceMatrix((const double**)datasetMatrix->d, rowNum, colNum, distMatrix.d);
		calculateConcentrationArray(metricSpace, &distMatrix, rowNum, v);
	} else if (metricSpace == HAMMING_SPACE) {
		calculateHammingDistanceMatrix((const char**)datasetMatrix->c, rowNum, colNum, distMatrix.i);
		calculateConcentrationArray(metricSpace, &distMatrix, rowNum, v);
	} else {
		distMatrixReference.d = datasetMatrix->d; // Reference the datasetMatrix.
		calculateConcentrationArray(metricSpace, &distMatrixReference, rowNum, v);
	}

	// Find the (#clusterNum) smallest elements.
	for (i = 0; i < clusterNum; i++) {
		// Find the first next element that does not already exist in the centroidIndexes array.
		for (j = 0; j < rowNum; j++) {
			if (v[j] != -1) {
				curMin = v[j];
				break;
			}
		}
		// Find the next smallest element that does not already exist in the centroidIndexes array.
		for (j = 0; j < rowNum; j++) {
			if ((v[j] <= curMin) && (v[j] != -1)) {
				curMin = v[j]; // The curMin is updated with the newest smallest element.
				centroidIndexes[i] = j; // The index of the smallest element is kept.
			}
		}
		v[centroidIndexes[i]] = -1; // Mark the current centroid to avoid it for the next the loop.
	}

	// Deallocate the space for the triangular matrix and the concentration array.
	if (metricSpace == VECTOR_SPACE) {
		for (i = 0; i < rowNum; i++) free((distMatrix.d)[i]);
		free(distMatrix.d);
	} else if (metricSpace == HAMMING_SPACE) {
		for (i = 0; i < rowNum; i++) free((distMatrix.i)[i]);
		free(distMatrix.i);
	}
	free(v);
}

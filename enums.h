#ifndef __ENUMS__
#define __ENUMS__

typedef enum METRIC_SPACES { VECTOR_SPACE, HAMMING_SPACE, MATRIX_SPACE, FUNCTION_SPACE, BIO_SPACE, COSINE } METRIC_SPACES;
typedef enum METRICS { EUCLIDEAN, MANHATTAN, COSINEz } METRICS;

typedef enum INITIALIZATION_METHODS {
	K_MEDOIDS = 1, // K-medoids++
	CONCENTRATE = 2 // Concentrate (Park-Jun)
} INITIALIZATION_METHODS;

typedef enum ASSIGNMENT_METHODS {
	PAM = 1, // PAM assignment
	LSH_DBH = 2 // Assignment by LSH/DBH: Reverse Approach
} ASSIGNMENT_METHODS;

typedef enum UPDATE_METHODS {
	LLOYD = 1, // Update a la Lloyd's (improved PAM)
	CLARANS = 2 // CLARANS Update
} UPDATE_METHODS;

#endif

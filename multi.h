#ifndef __MULTI__
#define __MULTI__

#include <gsl/gsl_matrix.h>

/****************************/
/*** forward declarations ***/
struct ham_hash;
struct ham_table;

struct vector_hashFun;
struct vec_table;

struct matrix_hashF;
struct matrix_table;

struct hi_func;
struct cos_vec_table;
/****************************/


// Union that contains all the different data types for each type of dataset.
typedef union multiMatrix {
	double** d; // double**
	gsl_matrix** b;// double*** to hold all 3D conformations of proteins
	char** c; // char**
} multiMatrix;

// Union that contains all the different data types used to calculate the distance between elements and place them into a 2D matrix.
typedef union multiDistMatrix {
	double** d; // double**
	int** i; // int**
} multiDistMatrix;

// Union that contains all the different data types used to calculate the distance between elements.
typedef union multiDist {
	double d; // double
	int i; // int
} multiDist;

typedef union multiLongDist {
	double d; // long double
	long long i; // long long
} multiLongDist;

typedef union multiHashes{
    struct ham_hash *Gh;
    struct vector_hashFun *Gv;
    struct matrix_hashF *Gm;
    struct hi_func *Gc;
} multiHashes;

typedef union multiTables{
    struct ham_table *HTables;
    struct vec_table *VTables;
    struct matrix_table *MTables;
    struct cos_vec_table* CTables;
} multiTables;

#endif

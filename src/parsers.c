#include <string.h>
#include <linux/limits.h>
#include <dlfcn.h>
#include <math.h>

#include "testing.h"
#if (__TEST__)
#include <CUnit/CUnit.h>
#endif

#include "parsers.h"
#include "distsNmore.h"
#include "hamming.h"
#include "vector.h"
#include "matrix.h"

#include "init.h"
#include "assignment.h"
#include "update.h"
#include "evaluation.h"
#include "tree.h"
#include "recommend.h"

// Those values may be changed by the user.
#define _r 67785 // (N * (N - 1) / 2) - 1
#define _hashFuncNum 4
#define _L 5
#define _k 5

// If it equal to 1 the NN LSH recommendation method will run. In any other case, the Clustering will run.
#define RECOMMENDATION_RUN_NN_LSH 1
// If it is equal to 1 and RECOMMENDATION_RUN_NN_LSH is equal to 0, then the Vector metric will be used. In any other case, the Cosine metric will be used.
#define RECOMMENDATION_RUN_CLUSTERING_VECTOR 1

#if (__TEST__)
static int timeStamp = -1;

void TEST_initTiming() {
	int timeDifference;

	if (timeStamp < 0) { // No starting time.
		CU_FAIL("No starting time!\n");
		return;
	}
	timeDifference = GetClockTimeInMilliSec() - timeStamp;
	if (timeDifference == 0) {
		CU_PASS("The time difference between the 2 timestamps is equal to zero!\n");
	} else {
		PrintTime(stdout, timeDifference);
		printf("\n");
	}
}
#endif

typedef struct CentroidPointsTree {
	int index;
	TreePtr tree;
} CentroidPointsTree;

void parse(FILE* input, FILE* output, int validation) {
	FILE* conformFile = NULL, *experimentsFile = NULL, *tmpFile = NULL;
	int isValid, update_flag;
	char path[PATH_MAX];
	char buffer[100];
	char answer;
	char* temp;
	int i, j, m, startTime, endTime;

	int P = -1;
	int userNum, itemNum;
	TreePtr usersTree, itemsTree;
	int userID, itemID, rating;

	int numConform = -1; // Amount of conformations
	int N = -1; // Amount of points in a conformation
	int r = _r, hashFuncNum = _hashFuncNum, L = _L, k = _k; // The definitions are copied to some local variables.
	int* centroidIndexes = NULL, *pointsClust = NULL;
	CentroidPointsTree* trees;
	long double avgS; // Average Silhouette value.
	int T = 1; // Amount of experiments

	multiLongDist oldCost;
	multiMatrix datasetMatrix, gslMatrix;
	multiTables mTs;
	multiHashes mHs;
	fpos_t pos;

	double*** conformations = NULL, ***conformationsBak = NULL; // [conformation_id, point_id, point_coordinate_id] (numConform x N x 3)
	gslMatrix.b = NULL;
	double** recommendations = NULL;
	double* avgR = NULL;
	int** fiveBest = NULL;

	if (input == NULL) {
		do {
			printf("* Please give the path to the input file:\n> ");
			scanf("%s", path);

			if ((input = fopen(path, "r")) == NULL) {
				printf("* The path is invalid!\n");
				isValid = 0;
			} else {
				isValid = 1;
			}
		} while (isValid == 0);
	}

	fgetpos(input, &pos);
	fscanf(input, "%[^\n]\n", buffer);
	if (strchr(buffer, '\t') == NULL) {
		if (sscanf(buffer, "%d", &i) != 1) return;
		numConform = P = i;

		fgetpos(input, &pos);
		fscanf(input, "%[^\n]\n", buffer);
		if (strchr(buffer, '\t') == NULL) {
			if (sscanf(buffer, "%d", &N) != 1) return;

			// This format corresponds to a conformations file.

			// Allocate space for all the required structs.
			conformations = safeMalloc(numConform * sizeof(double**));
			conformationsBak = safeMalloc(numConform * sizeof(double**));
			datasetMatrix.d = safeMalloc(numConform * sizeof(double*));
			gslMatrix.b = safeMalloc(numConform * sizeof(gsl_matrix*));
			pointsClust = safeMalloc(numConform * sizeof(int));
			for (i = 0; i < numConform; i++) {
				(datasetMatrix.d)[i] = safeMalloc(numConform * sizeof(double));
				(gslMatrix.b)[i] = gsl_matrix_alloc(N,3);
				conformations[i] = safeMalloc(N * sizeof(double*));
				conformationsBak[i] = safeMalloc(N * sizeof(double*));
				for (j = 0; j < N; j++) {
					conformations[i][j] = safeMalloc(3 * sizeof(double));
					conformationsBak[i][j] = safeMalloc(3 * sizeof(double));
				}
			}
		} else { // This format corresponds to a recommendations file.
			initTree(&usersTree);
			initTree(&itemsTree);
		}
	} else { // This format corresponds to a recommendations file.
		P = 20;
		initTree(&usersTree);
		initTree(&itemsTree);
	}

	if ((usersTree == NULL) && (conformations == NULL)) return;

	if ((usersTree != NULL) && (output == NULL)) {
		while ((tmpFile = fopen(path, "w")) == NULL) {
			printf("* Please give the path to the output file:\n> ");
			if ((temp = fgets(path, PATH_MAX - 1, stdin)) == NULL) return;
			path[strlen(buffer) - 1] = '\0'; // Replace '\n' with '\0'.
			fseek(stdin, 0, SEEK_END);
		}
	}

	if (usersTree != NULL) { // Parse a file that follows the recommendations format.
		fsetpos(input, &pos);
		while (fscanf(input, "%[^\n]\n", buffer) == 1) {
			if (sscanf(buffer, "%d\t%d\t%*d", &i, &j) != 2) continue;
			treeInsert(usersTree, i);
			treeInsert(itemsTree, j);
		}
		userNum = getTreeSize(usersTree);
		itemNum = getTreeSize(itemsTree);

		// Allocate space for the structs.
		recommendations = safeMalloc(userNum * sizeof(double*));
		pointsClust = safeMalloc(userNum * sizeof(int));
		avgR = safeMalloc(userNum * sizeof(double));
		fiveBest = safeMalloc(userNum * sizeof(int*));
		for (i = 0; i < userNum; i++) {
			recommendations[i] = safeMalloc(itemNum * sizeof(double));
			fiveBest[i] = safeMalloc(5 * sizeof(int));
			for (j = 0; j < itemNum; j++) recommendations[i][j] = 0;
			for (j = 0; j < 5; j++) fiveBest[i][j] = -1;
		}

		fsetpos(input, &pos);
		while (fscanf(input, "%[^\n]\n", buffer) == 1) {
			if (sscanf(buffer, "%d\t%d\t%d", &userID, &itemID, &rating) != 3) continue;

			if ((userID < 1) || (userID > userNum)) {
				printf("Non pre-processed input. All users must have an identifier between 1 and %d", userNum);
				return;
			}
			if ((itemID < 1) || (itemID > itemNum)) {
				printf("Non pre-processed input. All items must have an identifier between 1 and %d", itemNum);
				return;
			}
			recommendations[userID - 1][itemID - 1] = (double)rating;
		}
		calculateAvgRating(recommendations, avgR, userNum, itemNum); // Calculate the average rating for each user's ratings.

		while (hashFuncNum <= 0) {
			printf("* Please give the parameter hashFuncNum (default = 4):\n> ");
			scanf("%d%*[^\n]%*c", &hashFuncNum);
			fflush(stdin);
		}
		while (L <= 0) {
			printf("* Please give the parameter L (default = 5):\n> ");
			scanf("%d%*[^\n]%*c", &L);
			fflush(stdin);
		}

		datasetMatrix.d = recommendations;

		#if (RECOMMENDATION_RUN_NN_LSH == 1)
			initStructs(L, hashFuncNum, VECTOR_SPACE, &mHs, &mTs, &datasetMatrix, itemNum, userNum);
			startTime = GetClockTimeInMilliSec();
			NN_lsh_recommend(VECTOR_SPACE, L, &mHs, &mTs, userNum, itemNum, recommendations, P, fiveBest, avgR);
			startTime = GetClockTimeInMilliSec() - startTime;
			destroyStructs(L, VECTOR_SPACE, &mHs, &mTs);

			fprintf(output, "NN LSH (Euclidean)\n\n");
			for (i = 0; i < userNum; i++) {
				fprintf(output, "%d", i + 1);
				for (j = 0; j < 5; j++) {
					fprintf(output, "\t%d", (fiveBest[i][j] + 1));
				}
				fprintf(output, "\n");
			}
			fprintf(output, "Execution time: "); PrintTime(output, startTime); fprintf(output, "\n\n");

			initStructs(L, hashFuncNum, COSINE, &mHs, &mTs, &datasetMatrix, itemNum, userNum);
			startTime = GetClockTimeInMilliSec();
			NN_lsh_recommend(COSINE, L, &mHs, &mTs, userNum, itemNum, recommendations, P, fiveBest, avgR);
			startTime = GetClockTimeInMilliSec() - startTime;
			destroyStructs(L, COSINE, &mHs, &mTs);

			fprintf(output, "NN LSH (Cosine)\n\n");
			for (i = 0; i < userNum; i++) {
				fprintf(output, "%d", i + 1);
				for (j = 0; j < 5; j++) {
					fprintf(output, "\t%d", (fiveBest[i][j] + 1));
				}
				fprintf(output, "\n");
			}
			fprintf(output, "Execution time: "); PrintTime(output, startTime); fprintf(output, "\n");
		#else
			centroidIndexes = safeMalloc(k * sizeof(int)); // Allocate the space for storing the centroid indexes after initialization.
			#if (RECOMMENDATION_RUN_CLUSTERING_VECTOR == 1)
				while ((k <= 0) || (k >= itemNum)) {
					printf("* Please give the parameter k:\n> ");
					scanf("%d%*[^\n]%*c", &k);
				}

				startTime = GetClockTimeInMilliSec();
				k_medoids_init(VECTOR_SPACE, &datasetMatrix, userNum, itemNum, k, centroidIndexes);
				startTime = GetClockTimeInMilliSec() - startTime;

				// Reset the pointsClust array by setting 0 to every position except the ones that correspond to the centroid indexes and will be set as -1.
				for (i = 0; i < numConform; i++) pointsClust[i] = 0;
				for (i = 0; i < k; i++) pointsClust[centroidIndexes[i]] = -1;

				update_flag = 1;
				endTime = GetClockTimeInMilliSec();
				while (update_flag == 1) {
					oldCost = PAM_assignment(VECTOR_SPACE, pointsClust, &datasetMatrix, userNum, itemNum, k, centroidIndexes);
					update_flag = Lloyd_update(VECTOR_SPACE, pointsClust, &datasetMatrix, userNum, itemNum, k, centroidIndexes, oldCost);
				}
				endTime = GetClockTimeInMilliSec() - endTime;
				avgS = silhouette2(VECTOR_SPACE, &datasetMatrix, userNum, itemNum, k, pointsClust, centroidIndexes);
				printf("k = %d\navgS = %3.10Lf\n", k, avgS);

				cluster_recommend(COSINE, userNum, itemNum, datasetMatrix.d, fiveBest, avgR, pointsClust);
				fprintf(output, "Clustering (Cosine)\n\n");
				for (i = 0; i < userNum; i++) {
					fprintf(output, "%d", i + 1);
					for (j = 0; j < 5; j++) {
						fprintf(output, "\t%d", (fiveBest[i][j] + 1));
					}
					fprintf(output, "\n");
				}
				fprintf(output, "Execution time: "); PrintTime(output, startTime + endTime); fprintf(output, "\n");
			#else
				while ((k <= 0) || (k >= itemNum)) {
					printf("* Please give the parameter k:\n> ");
					scanf("%d%*[^\n]%*c", &k);
				}

				startTime = GetClockTimeInMilliSec();
				k_medoids_init(COSINE, &datasetMatrix, userNum, itemNum, k, centroidIndexes);
				startTime = GetClockTimeInMilliSec() - startTime;

				// Reset the pointsClust array by setting 0 to every position except the ones that correspond to the centroid indexes and will be set as -1.
				for (i = 0; i < numConform; i++) pointsClust[i] = 0;
				for (i = 0; i < k; i++) pointsClust[centroidIndexes[i]] = -1;

				update_flag = 1;
				endTime = GetClockTimeInMilliSec();
				while (update_flag == 1) {
					oldCost = PAM_assignment(COSINE, pointsClust, &datasetMatrix, userNum, itemNum, k, centroidIndexes);
					update_flag = Lloyd_update(COSINE, pointsClust, &datasetMatrix, userNum, itemNum, k, centroidIndexes, oldCost);
				}
				endTime = GetClockTimeInMilliSec() - endTime;
				avgS = silhouette2(COSINE, &datasetMatrix, userNum, itemNum, k, pointsClust, centroidIndexes);
				printf("k = %d\navgS = %3.10Lf\n", k, avgS);

				cluster_recommend(COSINE, userNum, itemNum, datasetMatrix.d, fiveBest, avgR, pointsClust);
				fprintf(output, "Clustering (Cosine)\n\n");
				for (i = 0; i < userNum; i++) {
					fprintf(output, "%d", i + 1);
					for (j = 0; j < 5; j++) {
						fprintf(output, "\t%d", (fiveBest[i][j] + 1));
					}
					fprintf(output, "\n");
				}
				fprintf(output, "Execution time: "); PrintTime(output, startTime + endTime); fprintf(output, "\n");
			#endif
			free(centroidIndexes);
		#endif
	} else { // Parse a file that follows the conformations format.
		// Create the output files.
		strncpy(path, "conform.dat", 12);
		while ((conformFile = fopen(path, "w")) == NULL) {
			printf("* Please give the path to the conform output file:\n> ");
			if ((temp = fgets(path, PATH_MAX - 1, stdin)) == NULL) return;
			path[strlen(path) - 1] = '\0'; // Replace '\n' with '\0'.
			fseek(stdin, 0, SEEK_END);
		}
		strncpy(path, "experim.dat", 12);
		while ((experimentsFile = fopen(path, "w")) == NULL) {
			printf("* Please give the path to the experiments output file:\n> ");
			if ((temp = fgets(path, PATH_MAX - 1, stdin)) == NULL) return;
			path[strlen(path) - 1] = '\0'; // Replace '\n' with '\0'.
			fseek(stdin, 0, SEEK_END);
		}
		strncpy(path, "tmp_file", 9);
		while ((tmpFile = fopen(path, "w+")) == NULL) { // Open the file for reading and writing.
			printf("* Please give the path for placing a temporary file:\n> ");
			if ((temp = fgets(path, PATH_MAX - 1, stdin)) == NULL) return;
			path[strlen(path) - 1] = '\0'; // Replace '\n' with '\0'.
			fseek(stdin, 0, SEEK_END);
		}

		for (i = 0; i < numConform; i++) {
			for (j = 0; j < N; j++) {
				// Parse the input file.
				fscanf(input, "%[^\n]\n", buffer);
				if (sscanf(buffer, "%lf\t%lf\t%lf", &(conformationsBak[i][j][0]), &(conformationsBak[i][j][1]), &(conformationsBak[i][j][2])) != 3) return;
				for (m = 0; m < 3; m++) conformations[i][j][m] = conformationsBak[i][j][m]; // Copy the conformations to the original matrix.
			}
			transform_conformations(conformations[i], N); // Transform the conformation.
			for (j = 0; j < N; j++) {
				for (m = 0; m < 3; m++) gsl_matrix_set((gslMatrix.b)[i], j, m, conformations[i][j][m]); // Copy the transformation to the gslMatrix as well.
			}
		}

		// Re-copy the conformations to the original matrix.
		for (i = 0; i < numConform; i++) {
			for (j = 0; j < N; j++) {
				for (m = 0; m < 3; m++) conformations[i][j][k] = conformationsBak[i][j][k];
				free(conformationsBak[i][j]);
			}
			free(conformationsBak[i]);
		}
		free(conformationsBak);

		while (1) {
			///////////////////////////////////
			// Begin clustering using cRMSD. //
			///////////////////////////////////
			while (1) {
				while ((k <= 0) || (k >= numConform)) {
					printf("* Please give the parameter k E (0, %d):\n> ", numConform);
					scanf("%d%*[^\n]%*c", &k);
				}
				centroidIndexes = safeMalloc(k * sizeof(int)); // Allocate the space for storing the centroid indexes after initialization.
				trees = safeMalloc(k * sizeof(CentroidPointsTree));

				k_medoids_init(BIO_SPACE, &gslMatrix, numConform, N, k, centroidIndexes);

				// Reset the pointsClust array by setting 0 to every position except the ones that correspond to the centroid indexes and will be set as -1.
				for (i = 0; i < numConform; i++) pointsClust[i] = 0;
				for (i = 0; i < k; i++) pointsClust[centroidIndexes[i]] = -1;

				update_flag = 1;
				while (update_flag == 1) {
					oldCost = PAM_assignment(BIO_SPACE, pointsClust, &gslMatrix, numConform, N, k, centroidIndexes);
					update_flag = Lloyd_update(BIO_SPACE, pointsClust, &gslMatrix, numConform, N, k, centroidIndexes, oldCost);
				}
				avgS = silhouette2(BIO_SPACE, &gslMatrix, numConform, N, k, pointsClust, centroidIndexes);
				fprintf(conformFile, "%d\n%3.10Lf\n", k, avgS);

				// Initialize the array of binary trees.
				for (i = 0; i < k; i++) {
					trees[i].index = -2;
					initTree(&(trees[i].tree));
				}
				for (i = 0; i < numConform; i++) {
					isValid = 0;
					if (pointsClust[i] == -1) {
						for (j = 0; j < k; j++) {
							if (trees[j].index != i) continue;
							treeInsert(trees[j].tree, i);
							isValid = 1;
							break;
						}
						if (isValid) continue;
						for (j = 0; j < k; j++) {
							if (trees[j].index != -2) continue;
							trees[j].index = i;
							treeInsert(trees[j].tree, i);
							break;
						}
					}
					for (j = 0; j < k; j++) {
						if (trees[j].index == pointsClust[i]) {
							treeInsert(trees[j].tree, i);
							break;
						}
						if (trees[j].index == -2) {
							trees[j].index = pointsClust[i];
							treeInsert(trees[j].tree, i);
							break;
						}
					}
				}
				for (i = 0; i < k; i++) {
					treePrint(trees[i].tree, conformFile); // Print the binary tree's nodes in ascending order.
					fprintf(conformFile, "\n");
					destroyTree(&(trees[i].tree)); // Destroy the binary tree.
				}
				free(trees); // Deallocate the space for the array of binary trees.

				printf("* Do you want to run another test for a different k (amount of clusters)? (y/n)\n> ");
				scanf("%c%*[^\n]%*c", &answer);
				if ((answer == 'y') || (answer == 'Y')) {
					fprintf(conformFile, "\n\n*****************************\n\n");
					k = -1;
					free(centroidIndexes);
				} else {
					break;
				}
			}

			///////////////////////////////////
			// Begin clustering using dRMSD. //
			///////////////////////////////////
			while (1) {
				while ((r <= 0) || (r >= N * (N - 1) / 2)) {
					printf("* Please give the parameter r E (0, %d):\n> ", N * (N - 1) / 2);
					scanf("%d%*[^\n]%*c", &r);
				}
				dRMSD(conformations, &datasetMatrix, numConform, N, r);

				startTime = GetClockTimeInMilliSec();
				k_medoids_init(MATRIX_SPACE, &datasetMatrix, numConform, numConform, k, centroidIndexes);
				startTime = GetClockTimeInMilliSec() - startTime;

				// Reset the pointsClust array by setting 0 to every position except the ones that correspond to the centroid indexes and will be set as -1.
				for (i = 0; i < numConform; i++) pointsClust[i] = 0;
				for (i = 0; i < k; i++) pointsClust[centroidIndexes[i]] = -1;

				update_flag = 1;
				endTime = GetClockTimeInMilliSec();
				while (update_flag == 1) {
					oldCost = PAM_assignment(MATRIX_SPACE, pointsClust, &datasetMatrix, numConform, numConform, k, centroidIndexes);
					update_flag = Lloyd_update(MATRIX_SPACE, pointsClust, &datasetMatrix, numConform, numConform, k, centroidIndexes, oldCost);
				}
				endTime = GetClockTimeInMilliSec() - endTime;
				avgS = silhouette2(MATRIX_SPACE, &datasetMatrix, numConform, numConform, k, pointsClust, centroidIndexes);
				fprintf(tmpFile, "%d\t%3.10Lf\t", r, avgS); PrintTime(tmpFile, startTime + endTime); fprintf(tmpFile, "\n");

				printf("* Do you want to run another test for a different r? (y/n)\n> ");
				scanf("%c%*[^\n]%*c", &answer);
				if ((answer == 'y') || (answer == 'Y')) {
					r = -1;
					T++;
				} else {
					fprintf(experimentsFile, "%d\n%d\n", k, T);
					fseek(tmpFile, 0, SEEK_SET); // Reset the reading/writing position of the temporary file.
					while (fgets(buffer, 100, tmpFile) != NULL)
						fputs(buffer, experimentsFile); // Copy the tmp file's contents to the experiments output file.
					break;
				}
			}

			// End of testing runs.
			printf("* Do you want to run the program again? (y/n)\n> ");
			scanf("%c%*[^\n]%*c", &answer);
			if ((answer == 'y') || (answer == 'Y')) {
				fprintf(conformFile, "\n\n*****************************\n\n");
				fprintf(experimentsFile, "\n\n*****************************\n\n");
				if (freopen(path, "w+", tmpFile) == NULL) {
					printf("* ERROR! The temporary file could not be reopened.");
					return;
				}
				r = hashFuncNum = L = k = -1;
				T = 1;
				free(centroidIndexes);
			} else {
				fclose(tmpFile); // Close the temporary file stream.
				if (remove(path) != 0) printf("* ERROR! The temporary file could not be cleaned."); // Delete the temporary file.
				break;
			}
		}
	}

	// Deallocate memory and close FILE streams.
	fclose(input);
	if (recommendations != NULL) {
		for (i = 0; i < userNum; i++) free(recommendations[i]); // FIX : CRASH
		free(recommendations);
		free(pointsClust);
		free(avgR);
		fclose(output);
	}

	if (conformations != NULL) {
		for (i = 0; i < numConform; i++) {
			for (j = 0; j < N; j++) free(conformations[i][j]);
			free(conformations[i]);
			free((datasetMatrix.d)[i]);
		}
		free(conformations);
		free(datasetMatrix.d);

		free(pointsClust);
		fclose(conformFile);
		fclose(experimentsFile);
	}
}

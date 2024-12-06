#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "testing.h"
#if (__TEST__)
	#include <CUnit/Basic.h>
#endif
#include "parsers.h"

#if (__TEST__)
	void TEST_initTiming();
#endif

int main(int argc, char* argv[]) {
	FILE* dataset = NULL; // Input file stream containing the dataset.
	FILE* output = NULL; // Output file stream for results.
	int validate = 0;
	int i = 1;

#if (__TEST__)
	CU_pSuite pSuite = NULL;

	if (CU_initialize_registry() != CUE_SUCCESS) return CU_get_error(); // Initialize the CUnit test registry.
	if ((pSuite = CU_add_suite("Test suite 1", NULL, NULL)) == NULL) {
		CU_cleanup_registry();
		return CU_get_error();
	}
	if ((CU_add_test(pSuite, "timing test", TEST_initTiming)) == NULL) {
		CU_cleanup_registry();
		return CU_get_error();
	}

	CU_basic_set_mode(CU_BRM_VERBOSE);
	CU_basic_run_tests();
	CU_cleanup_registry();
#endif

	srand(time(NULL));
	// Parse the given parameters.
	while (i < argc) {
		if (strncmp(argv[i], "-d", 2) == 0) {
			if (i + 1 > argc) break;
			dataset = fopen(argv[i + 1], "r");
			i++;
		} else if (strncmp(argv[i], "-o", 2) == 0) {
			if (i + 1 > argc) break;
			output = fopen(argv[i + 1], "w");
			i++;
		} else if (strncmp(argv[i], "-validate", 9) == 0) {
			validate = 1;
		}
		i++;
	}

	parse(dataset, output, validate); // Run the whole program logic.

	return 0;
}

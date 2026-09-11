//////////////////////////////////////////////////////////////////////////////
// Project: Sum numbers
// Author: Josiah Lawrence
// Date: 8 September 2026
// Description: This project sums all integers between 0 and a given number.
// It was developed at Chico State to practice parallel programming using 
// OpenMP and reduction variables.
//////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

int main(int argc, char *argv[]) {
	long long start_num = 0;
	long long end_num = 0;
	long long sum = 0;
	int thread_count = 1;
	int num_iterations = 1;

	/* Parse args */
	if (argc < 2) {
		printf("Usage: sumnums start_num end_num --threads thread_count --iterations num_iterations\n");
		exit(-1);
	} else if(argv[1] < 0 || argv[2] < 0 || argv[1] > argv[2]) {
		// printf("%d %d", argv[1], argv[2]);
		printf("Error: Number too small\n");
		exit(-1);
	} else {
		start_num = strtoll(argv[1], NULL, 10);
		end_num = strtoll(argv[2], NULL, 10);
	}

	for(long long i = 0; i < argc; i++) {
		if(!strcmp(argv[i], "--threads") && i < argc-1) {
			thread_count = atoi(argv[i+1]);
		} 
		if(!strcmp(argv[i], "--iterations") && i < argc-1) {
			num_iterations = atoi(argv[i+1]);
		} 
	}

	/* OpenMP directive */
#	pragma omp parallel for num_threads(thread_count) reduction(+: sum) collapse(2)
	/* Primary loop */
	for(int i = 0; i < num_iterations; i++) {
		for(long long j = start_num; j <= end_num; j++) {
			sum += j;
		}
	}

	/* Print results */
	printf("Sum = %lld\n", sum);

	return 0;
}


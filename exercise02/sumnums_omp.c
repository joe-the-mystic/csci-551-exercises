//////////////////////////////////////////////////////////////////////////////
// Project: Sum numbers
// Author: Josiah Lawrence
// Date: 8 September 2026
// Description: This project sums all integers between two given numbers.
// It was developed at Chico State to practice parallel programming using
// OpenMP and reduction variables. This version divides the count into
// subranges and assigns each thread to a single subrange.
//////////////////////////////////////////////////////////////////////////////

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int
main(int argc, char* argv[])
{
  long long lowest = 0;
  long long highest = 0;
  long long global_sum = 0;
  int thread_count = 1;
  int num_iterations = 1;
  long long num_terms = 0;
  struct timespec now, start;
  double fnow, fstart;

  clock_gettime(CLOCK_MONOTONIC, &start);
  fstart = (double)start.tv_sec + (double)start.tv_nsec / 1000000000.0;



  /* Parse args and define vars */
  if (argc < 3) {
    printf("Usage: sumnums lowest highest --threads thread_count "
           "--iterations num_iterations\n");
    exit(-1);
  } else {
    lowest = strtoll(argv[1], NULL, 10);
    highest = strtoll(argv[2], NULL, 10);
  }

  if (lowest < 0 || highest < 0 || lowest > highest) {
    // printf("%d %d", argv[1], argv[2]);
    printf("Error: Number too small\n");
    exit(-1);
  }

  for (long long i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "--threads") && i < argc - 1) {
      thread_count = atoi(argv[i + 1]);
    }
    if (!strcmp(argv[i], "--iterations") && i < argc - 1) {
      num_iterations = atoi(argv[i + 1]);
    }
  }

  num_terms = highest - lowest + 1;

  /* Counting loop */
  for (int iter = 0; iter < num_iterations; iter++) {
    global_sum = 0;
    printf("\nSumming %lld to %lld: iteration %d\n", lowest, highest, iter);
#pragma omp parallel for num_threads(thread_count) reduction(+ : global_sum)
    for (int th_idx = 0; th_idx < thread_count; th_idx++) {
      long long start = lowest + (th_idx * num_terms) / thread_count;
      long long end = lowest + ((th_idx + 1) * num_terms / thread_count - 1);
      long long local_sum = 0;
      for (long long i = start; i <= end; i++) {
        global_sum += i;
        local_sum += i;
      }
      printf("Thread %d summing range %lld to %lld: sum = %lld\n",
             th_idx,
             start,
             end,
             local_sum);
    }
  }

  printf("\nTEST COMPLETE:\n");
  printf("  Sum Total = %lld\n", global_sum);

  // Sum of arithmetic progression:
  // Verify that global_sum is (n*(2a+(n-1))/2, where a is the
  // lowest number, and n is the number of terms.
  long long verified_sum = (num_terms) * (2 * lowest + (num_terms - 1)) / 2;
  printf("  Verify: (n*(2a+(n-1))/2 = %lld\n", verified_sum);

  clock_gettime(CLOCK_MONOTONIC, &now);
  fnow = (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
  printf("Time elapsed: %f", fnow - fstart);
  
  return 0;
}

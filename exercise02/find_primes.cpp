//////////////////////////////////////////////////////////////////////////////
// Project: Find Primes
// Author: Josiah Lawrence
// Date: 12 September 2026
// Description: This project finds all prime numbers between 1 and a given
// number. It outputs the largest three primes, and the number of primes on the
// interval. In order to preserve time- and space-efficiency as well as
// thread safety, it first runs the sieve on a bytemap, then packs the
// data into a more compact bitmap. It was developed at Chico State to
// practice parallel programming using OpenMP.
//////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#define COUNT_PRIMES true
#define COMPUTE_PRIME_FACTORS true
#define WORD_LENGTH 8
#define THREAD_COUNT 24

/* A function to compact a binary bytemap into a bitmap.*/
void
compact(const std::vector<uint8_t>& bytemap,
        std::vector<uint8_t>& bitmap,
        unsigned long long highest)
{

# pragma omp parallel for num_threads(THREAD_COUNT)
  for (unsigned long long byte_idx = 0; byte_idx < bitmap.size(); byte_idx++) {
    uint8_t block;
    unsigned long long i;
    block = 0; // Prepare a block of 8 0's
    for (unsigned int bit_idx = 0; bit_idx < WORD_LENGTH; bit_idx++) {
      i = byte_idx * WORD_LENGTH + bit_idx;
      if (i <= highest && bytemap[i] == 1) {
        // Turn on the bit corresponding to the number i
        block |= (uint8_t(1) << bit_idx);
      }
    }
    bitmap[byte_idx] = block;
  }
}

/* Performs a prime sieve on the numbers up to the given number.
 *  Returns a bitmap where each bit shows whether its index is prime.*/
void
sieve(std::vector<uint8_t>& is_prime_bits, unsigned long long highest)
{
  std::vector<uint8_t> is_prime_bytes(highest + 1, 1);

  // By definition:
  is_prime_bytes[0] = 0;
  is_prime_bytes[1] = 0;

  // Find the upper limit for the sieve (sqrt(highest))
  unsigned long long limit = (unsigned long long)(sqrt(double(highest)));

  // Crossout loop
  // I used dynamic scheduling here because earlier threads do significantly
  // more work.
#pragma omp parallel for num_threads(THREAD_COUNT) schedule(dynamic, 1)
  for (unsigned long long i = 2; i <= limit; i++) {
    if (is_prime_bytes[i] == 1) {
      // Invalidate all multiples of current lowest prime
      for (unsigned long long j = i * i; j <= highest; j += i) {
        // Easily parallelizable: no race condition occurs unless two threads
        // write to the same index, in which case they both write the same value
        is_prime_bytes[j] = 0;
      }
    }
  }

  // Compact bytemap into bitmap
  compact(is_prime_bytes, is_prime_bits, highest);
}

/* Function to read a bit from an 8-bit binary bitmap */
bool
read_bit(const std::vector<uint8_t>& bitmap, unsigned long long read_idx)
{
  unsigned long long byte_idx = read_idx / WORD_LENGTH;
  unsigned long long bit_idx = read_idx % WORD_LENGTH;

  // Check for a 1 bit in the specified position.
  return (bitmap[byte_idx] & (uint8_t(1) << bit_idx)) > 0;
}

int
main(int argc, char* argv[])
{
  struct timespec now, start;
  double fnow, fstart;

  clock_gettime(CLOCK_MONOTONIC, &start);
  fstart = (double)start.tv_sec + (double)start.tv_nsec / 1000000000.0;

  /* Parse args and define vars */
  unsigned long long highest = 0;
  unsigned long long prime_count = 0;
  std::vector<unsigned long long> prime_list;

  if (argc < 2) {
    std::cout << "Usage: findprimes <highest>" << std::endl;
    exit(-1);
  } else {
    highest = strtoull(argv[1], NULL, 10);
  }

  std::vector<uint8_t> is_prime_bits(highest / 8 + 1);

  // Stage one: perform sieve on bytemap and compact into bitmap
  std::cout << "Looking for primes on interval 0 through " << highest << "."
            << std::endl;
  sieve(is_prime_bits, highest);
  clock_gettime(CLOCK_MONOTONIC, &now);
  fnow = (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
  std::cout << "Sieve completed. " << std::endl
            << "Elapsed time: " << (fnow - fstart) << std::endl;

#if COUNT_PRIMES
  // Stage two: count all primes and list the top three.
#pragma omp parallel for num_threads(THREAD_COUNT) reduction(+ : prime_count)
  for (unsigned long long i = 2; i <= highest; i++) {
    if (read_bit(is_prime_bits, i)) {
      prime_count++;
    }
  }
  std::cout << "Number of primes found: " << prime_count << std::endl;

  // Store the top 3 primes in prime_list and print
  prime_count = 0;
  prime_list.clear();
  unsigned long long i = highest;
  while (prime_count < 3 && i > 0) {
    if (read_bit(is_prime_bits, i)) {
      prime_list.push_back(i);
      prime_count++;
    }
    i--;
  }
  sort(prime_list.begin(), prime_list.end());
  std::cout << "Highest 3 primes: ";
  for (unsigned int i = 0; i < prime_list.size(); i++) {
    std::cout << prime_list[i] << " ";
  }
  std::cout << std::endl;

  // Find the product of the 2 highest primes
  std::cout << "Product of 2 highest primes: "
            << (prime_list[1] * prime_list[2]) << std::endl;

  clock_gettime(CLOCK_MONOTONIC, &now);
  fnow = (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
  std::cout << "Elapsed time: " << (fnow - fstart) << std::endl;
#endif

#if COMPUTE_PRIME_FACTORS
  // Stage three: find a set of prime factors of the given number(highest),
  // such that p1 * p2 = highest. Useful for computing the prime factors of a
  // large semiprime.
  std::vector<unsigned long long> prime_factors;
  for (unsigned long long i = 2; i * i <= highest; i++) {
    if (highest % i == 0 && read_bit(is_prime_bits, i) &&
        read_bit(is_prime_bits, highest / i)) {
      prime_factors.push_back(i);
      prime_factors.push_back(highest / i);
    }
  }
  // Example:
  // highest = 16
  // 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
  // x

  // Print prime factors
  std::cout << "Prime factors: ";
  for (unsigned long long i = 0; i < prime_factors.size(); i++) {
    std::cout << prime_factors[i] << " ";
  }
  std::cout << std::endl;

  clock_gettime(CLOCK_MONOTONIC, &now);
  fnow = (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
  std::cout << "Elapsed time: " << (fnow - fstart) << std::endl;

#endif
}

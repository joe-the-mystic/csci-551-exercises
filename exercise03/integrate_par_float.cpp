//////////////////////////////////////////////////////////////////////////////
// Project: Integrate
// Author: Josiah Lawrence
// Date: 28 September 2026
// Description: This project numerically integrates a function over a given
// interval, using MPI to divide the work in parallel. It was developed at
// CSU, Chico to practive distributed memory programming. This version uses
// a constant step size, and has been modified to work with MPI.
//////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
using namespace std;

// num_steps * step_size = interval
// num_steps = interval / step_size

#define DEBUG false
#define OUTPUT_PRECISION 7 // Scale of 6 to 15 recommended
#define NUM_STEPS_DEFAULT 180000

// Train Problem
#define INTERVAL_FROM_ZERO_DEFAULT (1800)
#define TSCALE (1800 / (2 * M_PI))
#define ASCALE 0.236589076381454

// Basic interval from 0 to pi
//#define INTERVAL_FROM_ZERO_DEFAULT (M_PI)
//#define TSCALE 1
//#define ASCALE 1

// Interval from 0 to 10
// #define INTERVAL_FROM_ZERO_DEFAULT 10
// #define TSCALE 1
// #define ASCALE 1


// Function Declarations
float accel(float x);
float vel(float x);
float pos(float x);
float left_riemann_sum(
  float (*integrand)(float), 
  const float &start, 
  const float &end, 
  const unsigned int &nstep, 
  const float &step_size  
);


// MAIN ///////////////////////////////////////////////////////////////////////
int
main(int argc, char* argv[])
{
  
  int my_rank, comm_sz, global_n, local_n, base, n_rem, first_step;
  float start, end, step_size, local_start, local_end;
  float local_int_area = 0.0; 
  float total_int_area = 0.0;
  float (*integrand)(float) = accel; // Function pointer to the integrand
  float (*antiderivative)(float) = vel; // FP to the antiderivative
  
  // Set float precision on cout
  cout << std::fixed << std::setprecision(OUTPUT_PRECISION);

  // Initialize MPI and get communicator info
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

  if (my_rank == 0) {
    if (argc >= 4) {
      start = atof(argv[1]);
      end = atof(argv[2]);
      global_n = atoi(argv[3]);
    } else {
      start = 0.0;
      end = INTERVAL_FROM_ZERO_DEFAULT;
      global_n = NUM_STEPS_DEFAULT;
      cout << "Usage: ./integrate_par <start> <end> <num_steps>" << endl
      << "Using defaults instead: " << start << " " << end << " " << global_n << endl;
    }
    for (int dest = 1; dest < comm_sz; dest++) {
      MPI_Send(&start, 1, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
      MPI_Send(&end, 1, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
      MPI_Send(&global_n, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
    }
  } else { /* my_rank != 0 */
    MPI_Recv(&start, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&end, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&global_n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  if (DEBUG) {
    cout << "Hello from process " << my_rank 
    << " of " << comm_sz << "!" << endl;
  }

  if (my_rank == 0) {
    cout << "Integrating from " << start << " to " << end << endl
    << "Number of steps: " << global_n << endl;
  } 

  /* Set up local sections to be computed in multiprocess. 
   * If global_n % comm_sz != 0, then leftover steps are assigned to ranks 0..n_rem. */
  base = global_n / comm_sz;
  n_rem = global_n % comm_sz; // remainder 
  local_n = base + (my_rank < n_rem ? 1 : 0); // Give an extra step to the first n_rem procs
  first_step = my_rank * base + min(my_rank, n_rem);
  step_size = (end - start) / global_n;
  local_start = start + first_step * step_size;
  local_end = local_start + local_n * step_size;

  /* Use a Riemann sum to estimate the definite integral using the area under
   * the curve.
   * First we sum the sub-intervals inside each section to get each process's local sum. */
  local_int_area = left_riemann_sum(integrand, local_start, local_end, local_n, step_size);
  if (DEBUG) {
    cout << "Process " << my_rank << " summed from " << local_start << " to " << local_end << ": " << local_int_area << endl;
  }

  /* Next we combine all local sums using MPI reduce. */
  MPI_Reduce(&local_int_area, &total_int_area, 1, MPI_FLOAT, MPI_SUM, 0,
         MPI_COMM_WORLD);
  
  
  if (my_rank == 0) {
    /* Print Riemann sum and check the numerical integration using the 2nd
     * fundamental theorem of Calculus. */
    cout << "Riemann sum: " << total_int_area << endl
         << "Verify by antiderivative: " << antiderivative(end) - antiderivative(start) << endl;
  }

  /* Shut down MPI */
  MPI_Finalize();
}

/* Acceleration function */
float
accel(float t)
{
  // return 10.0;
  return (sin(t / TSCALE) * ASCALE); // acceleration
}

/* Velocity function */
float
vel(float t)
{
  // return 10 * t;
  return (1.0 - cos(t / TSCALE)) * ASCALE * TSCALE; // velocity
}

/* Position function */
float
pos(float t)
{
  // return 5 * pow(t, 2);
  return (-TSCALE * (sin(t / TSCALE) + t)) * ASCALE * TSCALE; // position
}

/* Computes the middle riemann sum over the given interval */
float
left_riemann_sum(
  float (*integrand)(float), 
  const float &start, 
  const float &end, 
  const unsigned int &nstep, 
  const float &step_size  ) 
{
  
  float sum = 0.0;

  if (DEBUG) {
    cout << "Step size = " << step_size << " for " << nstep << " steps over interval "
         << start << " to " << end << "." << endl;
  }

  for (unsigned int idx = 0; idx < nstep; idx++) {
    float x = 0.0, fx = 0.0;
    x = start + (idx * step_size);
    fx = integrand(x);

    // Add the rectangle at the middle of each subinterval.
    sum += fx * step_size;
  }

  return sum;
}

// /* Computes the trapezoidal sum over the given interval. */
// double
// trapezoidal_sum(double start, double end)
// {
//   int nstep = lround((end - start) / STEP_SIZE);
//   double sum = 0.0;
//   double dx = (end - start) / nstep;
//
//   cout << "Step size = " << dx << " for " << nstep
//   << " steps over interval " << start << " to " << end << "." << endl;
//
//   for (int idx = 0; idx < nstep; idx++) {
//     double x = 0.0;
//     x = start + (idx * dx);
//
//     // Add the area of the trapezoid on the current subinterval.
//     sum += ( ((integrand(x) + integrand(x+dx)) / 2) * dx );
//   }
//
//   return sum;
// }

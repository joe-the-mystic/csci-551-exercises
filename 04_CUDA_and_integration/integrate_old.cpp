//////////////////////////////////////////////////////////////////////////////
// Project: Integrate
// Author: Josiah Lawrence
// Date: 14 October 2026
// Description: This project numerically integrates a function defined by data
// points in a CSV, using MPI to divide the work in parallel. It was developed
// at CSU, Chico to practice distributed memory programming. This version
// performs a double integration, handling a loop-carried dependency.
//////////////////////////////////////////////////////////////////////////////

#include "default_profile.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#define OUTPUT_PRECISION 9
#define STEPS_PER_UNIT 10000 // Number of steps per unit on the X-axis
#define DEBUG false
#define INITIAL_V 0
#define INITIAL_X 0
#define CLOCK true

// FUNCTIONS //////////////////////////////////////////////////////////////////

/* Integrates a discrete function defined by a vector. */
double
integrate_vec(
  const vector<double> &lookup_table, 
  const size_t start_idx,
  const size_t end_idx, 
  const unsigned int degree = 1, // 1st or 2nd degree integration 
  const double offset = 0) // Used only for degree 2
{
  double sum = 0.0;
  double step_size = 1.0 / STEPS_PER_UNIT;
  size_t nstep = end_idx - start_idx;

  if (degree == 0 || degree > 2) {
    cerr << "Error: invalid degree for integration: " << degree << endl;
    return 0;
  }

  if (end_idx >= lookup_table.size()) {
    cerr << "Error: end index " << end_idx << " out of bounds." << endl;
    return 0;
  } else if (start_idx >= lookup_table.size()) {
    cerr << "Error: start index " << start_idx << " out of bounds." << endl;
    return 0;
  } else if (start_idx >= end_idx) {
    return 0;
  }

  if (DEBUG) {
    cout << "Step size = " << step_size << " for " << nstep
         << " steps over interval " << start_idx / STEPS_PER_UNIT << " to "
         << end_idx / STEPS_PER_UNIT << "." << endl;
  }

  switch (degree) {
    case 1: {
      sum = 0.0;
      for (size_t i = start_idx; i < end_idx; i++) {
        double ffx = (((lookup_table[i] + lookup_table[i+1]) / 2) * step_size);
        sum += ffx;
      }
      break;
    }
    case 2: {
      double ffx = offset;
      sum = 0.0;
      for (size_t i = start_idx; i < end_idx; i++) {
        double dffx = ((lookup_table[i] + lookup_table[i+1]) / 2) * step_size;
        double ffx_next = ffx + dffx;
        sum += ((ffx + ffx_next) / 2) * step_size;
        ffx = ffx_next;
      }
      break;
    }
  }

  return sum;
} // INTEGRATE_VEC //

size_t
interpolate(vector<double> &in_vec, vector<double> &out_vec, unsigned int steps_per_unit) 
{
  if (in_vec.size() < 1) {
    cerr << "Vector too small: " << endl;
    return 1; 
  }
  
  size_t num_vals;
  double val, next_val, slope;
  
  num_vals = 0;
  for (unsigned int i = 0; i < in_vec.size() - 1; i++) {
    val = in_vec[i];
    next_val = in_vec[i + 1];
    slope = (next_val - val); // Slope is rise over run with run=1
    
    for (unsigned int j = 0; j < STEPS_PER_UNIT; j++) {
      out_vec.push_back(val + (j * (slope / STEPS_PER_UNIT)));
      num_vals++;
    }
  }
  out_vec.push_back(in_vec[in_vec.size()-1]); // get last element
  num_vals++;

  return num_vals;
}

void
interpolate_slice(
  const vector<double> &in_vec, 
  vector<double> &out_vec,
  size_t start_idx, 
  size_t end_idx)
{
  out_vec.resize(end_idx - start_idx + 1);
  for (size_t i = start_idx; i <= end_idx; i++) {
    // i is the global index, seg+j is the slice index
    size_t seg = i / STEPS_PER_UNIT;
    size_t j   = i % STEPS_PER_UNIT;
    if (seg + 1 >= in_vec.size()) {
      out_vec[i - start_idx] = in_vec.back();
    } else {
      double slope = in_vec[seg + 1] - in_vec[seg];
      out_vec[i - start_idx] = in_vec[seg] + j * (slope / STEPS_PER_UNIT);
    }
  }
}


double
find_max_v(
  const vector<double>& lookup_table, 
  const size_t start_idx, 
  const size_t end_idx,
  const double offset = 0)
{
  double sum;
  double max;
  double step_size = 1.0 / STEPS_PER_UNIT;

  if (end_idx >= lookup_table.size()) {
    cerr << "Error: end index " << end_idx << " out of bounds." << endl;
    return 0;
  } else if (start_idx >= lookup_table.size()) {
    cerr << "Error: start index " << start_idx << " out of bounds." << endl;
    return 0;
  } else if (start_idx >= end_idx) {
    return 0;
  }

  sum = offset;
  max = offset;
  for (size_t i = start_idx; i < end_idx; i++) {
    double ffx = (((lookup_table[i] + lookup_table[i + 1]) / 2) * step_size);
    sum += ffx;
    if (sum > max) max = sum;
  }

  return max;
}

/* Takes in data from a CSV and stores it in a vector. 
 * Returns size of vector. */
unsigned int
read_csv_to_vec(const string &fname_in, vector<double> &vec_out) 
{
  ifstream fin;
  if (!vec_out.empty()) vec_out.clear();

  fin.open(fname_in);
  if (!fin.is_open()) {
    cerr << "Couldn't read file: " << fname_in << endl;
    return 0; 
  }

  // Just takes the first column of every row
  for (string line; getline(fin, line);) {
    stringstream ss(line);
    string buf;
    getline(ss, buf, ',');
    if (!buf.empty()){
      vec_out.push_back(stod(buf));
    }
  }

  fin.close();
  return vec_out.size();
} // READ_CSV_TO_VEC //



// MAIN ///////////////////////////////////////////////////////////////////////

int 
main(int argc, char* argv[]) {
  cout << fixed << setprecision(OUTPUT_PRECISION);

  // Integration vars
  double local_initial_v = 0, local_delta_v, local_delta_x, local_max_v;
  vector<double> global_table, acc_table;
  size_t table_size = 0;
  size_t global_num_steps, local_num_steps, chunk_size, chunk_residual;
  size_t local_start_idx, local_end_idx;
  double global_final_v = 0, global_final_x = 0, global_max_v = 0;

  // MPI Vars
  int my_rank, comm_sz;

  // Initialize clock
#if CLOCK
  struct timespec now, start;
  double fnow, fstart, finterp, finit, finteg;
  clock_gettime(CLOCK_MONOTONIC, &start);
  fstart = (double)start.tv_sec + (double)start.tv_nsec / 1000000000.0;
#endif

// Initialize MPI and get communicator info
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
if (DEBUG && my_rank == 0) cout << "Comm size: " << comm_sz << endl;


// Check clock time
#if CLOCK
  if (my_rank == 0) {
    clock_gettime(CLOCK_MONOTONIC, &now);
    fnow = (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
    finit = fnow;
    cout << "Initialization time: " << (finit - fstart) << endl;
  }
#endif

  if (argc > 2) {
    if (my_rank == 0) 
      cout << "Usage: ./integrate in_file.csv" << endl;
    MPI_Finalize();
    return 1;
  }

  if (argc == 2) {
    // Get table from file
    if (my_rank == 0) {
      table_size = read_csv_to_vec(argv[1], global_table);
    }
    MPI_Bcast(&table_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    global_table.resize(table_size);
    MPI_Bcast(global_table.data(), table_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  } else {
    global_table = DefaultProfile;
  }

  // Divide up chunks for each process. 
  global_num_steps = (global_table.size() - 1) * STEPS_PER_UNIT;
  if (comm_sz > static_cast<int>(global_num_steps)) {
    chunk_size = 1;
    chunk_residual = 0;
    // Procs without work idle
    if (my_rank >= static_cast<int>(global_num_steps)) {
      chunk_size = 0;
      chunk_residual = 0;
      MPI_Finalize();
      return 0;
    }
  } else {
    chunk_size = global_num_steps / comm_sz;
    chunk_residual = global_num_steps % comm_sz;
  }
  local_start_idx = my_rank * chunk_size;
  local_end_idx = (my_rank + 1) * chunk_size;
  // The last process gets the residual.
  if (my_rank == comm_sz - 1) local_end_idx += chunk_residual;
  local_num_steps = local_end_idx - local_start_idx;


  interpolate_slice(global_table, acc_table, local_start_idx, local_end_idx);
  

  // Check clock time
#if CLOCK
  if (my_rank == 0) {
    clock_gettime(CLOCK_MONOTONIC, &now);
    fnow = (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
    finterp = fnow;
    cout << "Interpolation time: " << (finterp - finit) << endl;
  }
#endif

  
  // Compute local change in velocity over the chunk
  local_delta_v = integrate_vec(acc_table, 0, local_num_steps);

  // Get initial velocity for each rank by using Exscan
  // For each rank r, Exscan gets the sum of the values from ranks 0..(r-1).
  local_initial_v = 0.0;
  MPI_Exscan(&local_delta_v, &local_initial_v, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (my_rank == 0) local_initial_v = 0.0;
  local_initial_v += INITIAL_V;
  
  // Compute local change in position over the chunk.
  local_delta_x = integrate_vec(acc_table, 0, local_num_steps, 2, local_initial_v);

  // Compute local maximum velocity
  local_max_v = find_max_v(acc_table, 0, local_num_steps, local_initial_v);

  // Compute final velocity
  MPI_Reduce(&local_delta_v, &global_final_v, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (my_rank == 0) global_final_v += INITIAL_V;
  
  // Compute final position
  MPI_Reduce(&local_delta_x, &global_final_x, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (my_rank == 0) global_final_x += INITIAL_X;

  // Compute global max velocity
  MPI_Reduce(&local_max_v, &global_max_v, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  // Check clock time
#if CLOCK
  if (my_rank == 0) {
    clock_gettime(CLOCK_MONOTONIC, &now);
    fnow = (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
    finteg = fnow;
    cout << "Integration time: " << (finteg - finterp) << endl;
  }
#endif
  
  // Print results
  if (my_rank == 0) {
    cout << "Results:\n" << "\tFinal Velocity: " << global_final_v << "\n"
         << "\tFinal Position: " << global_final_x << "\n"
         << "\tMaximum Velocity: " << global_max_v << endl;
  }

  // Check clock time
#if CLOCK
  if (my_rank == 0) {
    clock_gettime(CLOCK_MONOTONIC, &now);
    fnow = (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
    cout << "Total time elapsed: " << (fnow - fstart) << endl;
  }
#endif
  
  MPI_Finalize();
  return 0;
} // MAIN //

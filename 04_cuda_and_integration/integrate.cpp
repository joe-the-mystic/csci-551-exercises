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
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#define OUTPUT_PRECISION 9
#define STEPS_PER_UNIT 1000 // Number of steps per unit on the X-axis
#define DEBUG false
#define INITIAL_V 0
#define INITIAL_X 0
#define CLOCK true

// FUNCTIONS //////////////////////////////////////////////////////////////////

void
interpolate_slice(
  const vector<long double> &in_vec, 
  vector<long double> &out_vec,
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
      long double slope = in_vec[seg + 1] - in_vec[seg];
      out_vec[i - start_idx] = in_vec[seg] + j * (slope / STEPS_PER_UNIT);
    }
  }
}

/* Takes in data from a CSV and stores it in a vector. 
 * Returns size of vector. */
unsigned int
read_csv_to_vec(const string &fname_in, vector<long double> &vec_out) 
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
  long double local_initial_v = 0, local_delta_v, local_delta_x, local_max_v;
  vector<long double> acc_table;
  size_t table_size = 0;
  size_t global_num_steps, local_num_steps, chunk_size, chunk_residual;
  size_t local_start_idx, local_end_idx;
  long double global_final_v = 0, global_final_x = 0, global_max_v = 0;
  long double step_size;

  // MPI Vars
  int my_rank, comm_sz;

  // Initialize clock
#if CLOCK
    struct timespec now, start;
    double fnow, fstart, finit, finteg;
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
  MPI_Barrier(MPI_COMM_WORLD);
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
    // If a file is specified on the command line, get table from file.
    // Otherwise use the global_table defined in default_profile.h.
    if (my_rank == 0) {
      table_size = read_csv_to_vec(argv[1], global_table);
    }
    MPI_Bcast(&table_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    global_table.resize(table_size);
    MPI_Bcast(global_table.data(), table_size, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
  } 

  // Divide up chunks for each process. 
  global_num_steps = (global_table.size() - 1) * STEPS_PER_UNIT;
  if (comm_sz > static_cast<int>(global_num_steps)) {
    if (my_rank == 0) cerr << "Table size " << global_num_steps << " smaller than number of processes. Aborting..." << endl;
    MPI_Finalize();
    return 1;
  } else {
    chunk_size = global_num_steps / comm_sz;
    chunk_residual = global_num_steps % comm_sz;
  }
  local_start_idx = my_rank * chunk_size;
  local_end_idx = (my_rank + 1) * chunk_size;
  // The last process gets the residual.
  if (my_rank == comm_sz - 1) local_end_idx += chunk_residual;
  local_num_steps = local_end_idx - local_start_idx;
  step_size = 1.0 / STEPS_PER_UNIT;

  interpolate_slice(global_table, acc_table, local_start_idx, local_end_idx);

  // Compute chunk results
  // Relative to local initial v:
  local_delta_v = 0.0;
  local_delta_x = 0.0;
  local_max_v = 0.0;

  if (DEBUG) {
    cout << "Step size = " << step_size << " for " << local_num_steps
         << " steps over interval " << local_start_idx / STEPS_PER_UNIT << " to "
         << local_end_idx / STEPS_PER_UNIT << "." << endl;
  }

  for (size_t i = 0; i < local_num_steps; i++) {
    long double v_next = local_delta_v + 0.5 * (acc_table[i] + acc_table[i+1]) * step_size;
    local_delta_x += 0.5 * (local_delta_v + v_next) * step_size;
    local_delta_v = v_next;
    if (local_delta_v > local_max_v) local_max_v = local_delta_v;
  }
  
  // Get initial velocity for each rank by using Exscan
  // For each rank r, Exscan gets the sum of the values from ranks 0..(r-1).
  MPI_Exscan(&local_delta_v, &local_initial_v, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (my_rank == 0) local_initial_v = 0.0;
  local_initial_v += INITIAL_V;

  // Add offsets
  local_delta_x += local_initial_v * local_num_steps * step_size;
  local_max_v += local_initial_v;

  // Compute final velocity
  MPI_Reduce(&local_delta_v, &global_final_v, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (my_rank == 0) global_final_v += INITIAL_V;
  
  // Compute final position
  MPI_Reduce(&local_delta_x, &global_final_x, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (my_rank == 0) global_final_x += INITIAL_X;

  // Compute global max velocity
  MPI_Reduce(&local_max_v, &global_max_v, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  // Check clock time
#if CLOCK
  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank == 0) {
    clock_gettime(CLOCK_MONOTONIC, &now);
    fnow = (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
    finteg = fnow;
    cout << "Integration time: " << (finteg - finit) << endl;
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
  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank == 0) {
    clock_gettime(CLOCK_MONOTONIC, &now);
    fnow = (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
    cout << "Total time elapsed: " << (fnow - fstart) << endl;
  }
#endif
  
  MPI_Finalize();
  return 0;
} // MAIN //

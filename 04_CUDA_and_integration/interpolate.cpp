//////////////////////////////////////////////////////////////////////////////
// Project: Interpolate
// Author: Josiah Lawrence
// Date: 12 October 2026
// Description: This program reads values from a CSV to get a function profile,
// then uses linear interpolation to turn 1 Hz data into 10 Hz data.
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#define READ_COLUMN 1
#define OUTPUT_PRECISION 9
#define STEPS_PER_UNIT 1000 // Number of steps to be interpolated between each point.

int 
main(int argc, char* argv[]) {
  vector<double> in_vec; 
  vector<double> out_vec;
  if (argc < 4) {
    cerr << "Usage: ./interpolate in_file.csv out_file.csv outfile.h" << endl;
    return 1;
  }
  string fname_in = argv[1];
  string fname_out = argv[2];
  string fname_h_out = argv[3];
  ifstream fin;
  ofstream fout;
  ofstream fhout;

  /* Open and read the given column of the CSV. Store as Vector. */
  fin.open(fname_in);
  if (!fin.is_open()) {
    cerr << "Couldn't read file: " << fname_in << endl;
    return 1; 
  }

  for (string line; getline(fin, line);) {
    stringstream ss(line);
    string dummy, buf;
    for (unsigned int i = 1; i < READ_COLUMN; i++) {
      getline(ss, dummy, ',');
    }
    getline(ss, buf, ',');
    if (!buf.empty()){
      in_vec.push_back(stod(buf));
    }
  }

  /* Create a new vector with interpolated data. */
  if (in_vec.size() < 1) {
    cerr << "Vector too small: " << endl;
    return 1; 
  }

  for (unsigned int i = 0; i < in_vec.size() - 1; i++) {
    double val = in_vec[i];
    double next_val = in_vec[i + 1];
    double slope = (next_val - val); // Slope is rise over run
    
    for (unsigned int j = 0; j < STEPS_PER_UNIT; j++) {
      out_vec.push_back(val + (j * (slope / STEPS_PER_UNIT)));
    }


  }
  out_vec.push_back(in_vec[in_vec.size()-1]); // get last element

  fin.close();

  /* Write output vector to CSV */
  fout.open(fname_out);
  fout << fixed << setprecision(OUTPUT_PRECISION);
  if (!fout.is_open()) {
    cerr << "Couldn't create file: " << fname_out << endl;
    return 1; 
  }
  for (unsigned int i = 0; i < out_vec.size(); i++) {
    fout << out_vec[i] << '\n';
  }
  fout.close();

  /* Write output vector to header file */
  fhout.open(fname_h_out);
  fhout << fixed << setprecision(OUTPUT_PRECISION);
  if (!fhout.is_open()) {
    cerr << "Couldn't create file: " << fname_h_out << endl;
    return 1;
  }
  fhout << "double DefaultProfile[] = {\n\t";
  for (unsigned int i = 0; i < out_vec.size(); i++) {
    fhout << out_vec[i];
    if (i != out_vec.size() - 1) {
      fhout << (i%3==0 ? ",\n\t" : ", ");
    } else {
      fhout << "\n};" << endl;
    }
  }
  fhout.close();
  

}
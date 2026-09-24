//////////////////////////////////////////////////////////////////////////////
// Project: Copy PGM
// Author: Josiah Lawrence
// Date: 15 September 2026
// Description: This project parses and reads all PGM files in a given
// directory, then rotates them and writes them to a new directory. It was 
// developed at Chico State to practice vector and image processing in parallel.
//////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using std::vector;

#define THREAD_COUNT 1

struct Image
{
  std::string name;
  unsigned int width = 0;
  unsigned int height = 0;
  int maxval = 0;
  vector<vector<unsigned char>> data; // index by [row][col]
};

// Skips whitespace and comment lines in the PGM header
static void
skipWhitespaceAndComments(std::istream& in_fstream)
{
  char c;
  while (in_fstream.peek()) {
    c = in_fstream.peek();
    if (c == '#') {
      // skip to end of line
      std::string buf;
      std::getline(in_fstream, buf);
    } else if (std::isspace(static_cast<unsigned char>(c))) {
      in_fstream.get(c);
    } else {
      break;
    }
  }
}

int
readPGM(Image& img, const std::string& filename) 
{
  std::string magic;

  std::ifstream file(filename, std::ios::binary);
  if (!file) {
    std::cerr << "Error: could not open file " << filename << std::endl;
    exit(-1);
  }

  // Extract magic number
  // Should read "P5" for PGM file
  file >> magic;
  if (magic != "P5") {
    std::cerr << "Error: unsupported file format (PGM only)." << std::endl;
    exit(-1);
  }

  // Get image properties
  skipWhitespaceAndComments(file);
  file >> img.width;
  skipWhitespaceAndComments(file);
  file >> img.height;
  skipWhitespaceAndComments(file);
  file >> img.maxval;
  file.get();

  if (img.maxval != 255) {
    std::cerr << "Error: incorrect maxval: " << img.maxval << std::endl;
    exit(-1);
  }

  // Read in binary data
  img.data.resize(img.height, vector<unsigned char>(img.width, 0));
//#pragma omp parallel for num_threads(THREAD_COUNT)
  for (unsigned int i = 0; i < img.height; i++) {
    unsigned char buf;
    for (unsigned int j = 0; j < img.width; j++) {
      file.read(reinterpret_cast<char*>(&buf), sizeof(unsigned char));
      img.data[i][j] = buf; 
    }
  }

  file.close();
  return 0;
}

int
writePGM(const std::string& filename, const Image& img)
{
  std::string magic = "P5";
  std::ofstream file(filename, std::ios::binary);

  if (!file) {
    std::cerr << "Error: could not open file." << std::endl;
    exit(-1);
  }

  // Construct Header
  file << magic << std::endl
  << "# Created by Joe Lawrence" << std::endl
  << img.width << " " << img.height << std::endl
  << img.maxval << std::endl;

  // Write binary data
//#pragma omp parallel for num_threads(THREAD_COUNT)
  for (unsigned int j = 0; j < img.height; j++) {
    unsigned char buf;
    for (unsigned int i = 0; i < img.width; i++) {
      buf = img.data[j][i];
      file.write(reinterpret_cast<const char *>(&buf), sizeof(unsigned char));
    }
  }

  file.close();
  return(0);
}

// Transposes the input image by swapping the columns and rows.
int 
transposePGM (Image& img) 
{
  // Create a new vector with height and width swapped
  vector<vector<unsigned char>> new_data(img.width, vector<unsigned char>(img.height, 0));

  // Assign values to the new vector so V2[i][j] = V1[j][i]
  for (unsigned int j = 0; j < img.height; j++) {
    for (unsigned int i = 0; i < img.width; i++) {
      new_data[i][j] = img.data[j][i];
    }
  }
  img.data = new_data;
  int tmp = img.width;
  img.width = img.height;
  img.height = tmp;
  return 0;
}

// Pivots the input image by swapping the values about the central
// vertical axis
int
verticalPivotPGM (Image& img) 
{
  // Create a new identical empty vector
  vector<vector<unsigned char>> new_data(img.height, vector<unsigned char>(img.width, 0));
  // Assign values to the new vector so the columns are assigned in reverse order
  for (unsigned int j = 0; j < img.height; j++) {
    for (unsigned int i = 0; i < img.width; i++) {
      new_data[j][(img.width-1)-i] = img.data[j][i];
    }
  }
  img.data = new_data;
  return 0;
}

// Pivots the input image by swapping the values about the central
// horizontal axis
int
horizontalPivotPGM (Image& img) 
{
  // Create a new identical empty vector
  vector<vector<unsigned char>> new_data(img.height, vector<unsigned char>(img.width, 0));
  // Assign values to the new vector so the rows are assigned in reverse order
  for (unsigned int j = 0; j < img.height; j++) {
    for (unsigned int i = 0; i < img.width; i++) {
      new_data[(img.height-1)-j][i] = img.data[j][i];
    }
  }
  img.data = new_data;
  return 0;
}

// Rotate a PGM right 90 degrees
int
rightRotatePGM (Image& img) {
  transposePGM(img);
  verticalPivotPGM(img);
  return 0;
}

// Rotate a PGM left 90 degrees
int
leftRotatePGM (Image& img) {
  transposePGM(img);
  horizontalPivotPGM(img);
  return 0;
}

int
main(int argc, char* argv[]) {
  std::string in_filename;
  std::string out_filename;
  Image img;
  bool verbose = false;

  // Parse Args
  if (argc < 3) {
    std::cout << "Usage: ./rotate_pgm <in_filename> <out_filename> --{options}" << std::endl
    << "Options: --verbose|-v" << std::endl;
    exit(-1);
  } else {
    in_filename = argv[1];
    out_filename = argv[2];
  }
  for (int i = 3; i < argc; i++) {
    if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "-v") {
      verbose = true;
    }
  }

  // Read file as image
  readPGM(img, in_filename);

  if (verbose) {
    // Print file header info
    std::cout << "File type: PGM" << std::endl
      << "Width: " << img.width << " Height: " << img.height << std::endl
      << "Max color value: " << img.maxval << std::endl;
  }

  // Transpose PGM
  rightRotatePGM(img);

  // Write out image to file
  writePGM(out_filename, img);
  std::cout << "File copied from " << in_filename << 
    " to " << out_filename << "." << std::endl;

}



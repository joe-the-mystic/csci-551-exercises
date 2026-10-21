//////////////////////////////////////////////////////////////////////////////
// Project: Copy PGM
// Author: Josiah Lawrence
// Date: 15 September 2026
// Description: This project parses and reads a .pgm file, then writes out a
// copy of it. It was developed at Chico State to practice vector and image
// processing in parallel.
//////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using std::vector;

#define THREAD_COUNT 1

struct Image
{
  int width = 0;
  int height = 0;
  int maxval = 0;
  vector<vector<unsigned char>> data; // index by [row][col]
};

// Skips whitespace and comment lines in the PGM header
static void
skipWhitespaceAndComments(std::istream& input_fs)
{
  char c;
  while (input_fs.peek()) {
    c = input_fs.peek();
    if (c == '#') {
      // skip to end of line
      std::string buf;
      std::getline(input_fs, buf);
    } else if (std::isspace(static_cast<unsigned char>(c))) {
      input_fs.get(c);
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
    std::cerr << "Error: could not open file." << std::endl;
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

  std::cout << img.width << " "
    << img.height << " "
    << img.maxval << std::endl;

  if (img.maxval != 255) {
    std::cerr << "Error: incorrect maxval: " << img.maxval << std::endl;
    exit(-1);
  }

  // Read in binary data
  img.data.resize(img.height, vector<unsigned char>(img.width, 0));
//#pragma omp parallel for num_threads(THREAD_COUNT)
  for (int i = 0; i < img.height; i++) {
    unsigned char buf;
    for (int j = 0; j < img.width; j++) {
      file.read(reinterpret_cast<char*>(&buf), sizeof(unsigned char));
      img.data[i][j] = buf; 
    }
  }

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
  for (int i = 0; i < img.height; i++) {
    unsigned char buf;
    for (int j = 0; j < img.width; j++) {
      buf = img.data[i][j];
      file.write(reinterpret_cast<const char *>(&buf), sizeof(unsigned char));
    }
  }

  file.close();
  return(0);
}

int
main(int argc, char* argv[]) {
  std::string in_filename;
  std::string out_filename;
  Image img;
  bool verbose = false;

  // Parse Args
  if (argc < 3) {
    std::cout << "Usage: ./copy_pgm <in_filename> <out_filename> --{options}" << std::endl
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

  // Write out image to file
  writePGM(out_filename, img);
  std::cout << "File copied from " << in_filename << 
    " to " << out_filename << "." << std::endl;

}


// Matrix flip functions? Transpose and Pivot.
#include "historical_pruned_landmark_labeling.h"
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
using namespace std;

int main(int argc, char **argv) {
  if (argc != 3) {
    cerr << "usage: query GRAPH(input) INDEX(output)" << endl;
    exit(EXIT_FAILURE);
  }

  historical_pruned_landmark_labeling hpll;
  hpll.construct_index(argv[1]);

  ofstream ofs(argv[2]);
  HPLL_CHECK(ofs);
  boost::archive::binary_oarchive oa(ofs);
  oa << hpll;
  ofs.close();

  exit(EXIT_SUCCESS);
}

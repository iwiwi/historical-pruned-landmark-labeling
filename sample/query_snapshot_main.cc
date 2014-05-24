#include "historical_pruned_landmark_labeling.h"
#include <iostream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
using namespace std;

int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "usage: query_snapshot INDEX" << endl;
    exit(EXIT_FAILURE);
  }

  historical_pruned_landmark_labeling hpll;

  ifstream ifs(argv[1]);
  HPLL_CHECK(ifs);
  boost::archive::binary_iarchive ia(ifs);
  ia >> hpll;

  for (int u, v, t; cin >> u >> v >> t; ) {
    int d = hpll.query_snapshot(u, v, t);
    cout << d << endl;
  }
}

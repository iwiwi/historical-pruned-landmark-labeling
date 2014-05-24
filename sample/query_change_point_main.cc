#include "historical_pruned_landmark_labeling.h"
#include <iostream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
using namespace std;

int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "usage: query_change_point INDEX" << endl;
    exit(EXIT_FAILURE);
  }

  historical_pruned_landmark_labeling hpll;

  ifstream ifs(argv[1]);
  HPLL_CHECK(ifs);
  boost::archive::binary_iarchive ia(ifs);
  ia >> hpll;

  for (int u, v; cin >> u >> v; ) {
    vector<pair<int, int>> cp;
    hpll.query_change_points(u, v, cp);
    for (auto p : cp) {
      cout << p.first << ":" << p.second << "\t";
    }
    cout << endl;
  }
}

#include "historical_pruned_landmark_labeling.h"
#include <iostream>
using namespace std;

int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "usage: query GRAPH" << endl;
    exit(EXIT_FAILURE);
  }

  historical_pruned_landmark_labeling hpll;
  hpll.construct_index(argv[1]);
  for (int u, v, t; cin >> u >> v >> t; ) {
    int d = hpll.query_snapshot(u, v, t);
    cout << d << endl;
  }
}

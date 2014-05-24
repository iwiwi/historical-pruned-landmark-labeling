#include "historical_pruned_landmark_labeling.h"
#include "../lib/gtest/gtest.h"
#include <climits>
#include <cstdlib>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
using namespace std;

TEST(hpll, triangle) {
  vector<tuple<int, int, int>> es ={
    make_tuple(1, 0, 1),
    make_tuple(2, 1, 2),
    make_tuple(3, 0, 2),
  };

  historical_pruned_landmark_labeling hpll;
  hpll.construct_index(es);

  ASSERT_EQ(hpll.query_snapshot(0, 1, 0), -1);
  ASSERT_EQ(hpll.query_snapshot(0, 2, 0), -1);
  ASSERT_EQ(hpll.query_snapshot(1, 2, 0), -1);

  ASSERT_EQ(hpll.query_snapshot(0, 1, 1), 1);
  ASSERT_EQ(hpll.query_snapshot(0, 2, 1), -1);
  ASSERT_EQ(hpll.query_snapshot(1, 2, 1), -1);

  ASSERT_EQ(hpll.query_snapshot(0, 1, 2), 1);
  ASSERT_EQ(hpll.query_snapshot(0, 2, 2), 2);
  ASSERT_EQ(hpll.query_snapshot(1, 2, 2), 1);

  ASSERT_EQ(hpll.query_snapshot(0, 1, 3), 1);
  ASSERT_EQ(hpll.query_snapshot(0, 2, 3), 1);
  ASSERT_EQ(hpll.query_snapshot(1, 2, 3), 1);

  ASSERT_EQ(hpll.query_snapshot(0, 1, 123123123), 1);
  ASSERT_EQ(hpll.query_snapshot(0, 2, 123123123), 1);
  ASSERT_EQ(hpll.query_snapshot(1, 2, 123123123), 1);

  ASSERT_EQ(hpll.query_snapshot(0, 123123123, 0), -1);
  ASSERT_EQ(hpll.query_snapshot(-1, 0, 0), -1);

  {
    vector<pair<int, int>> cp, ans =
      {make_pair(0, -1), make_pair(2, 2), make_pair(3, 1)};
    hpll.query_change_points(0, 2, cp);
    ASSERT_EQ(cp, ans);
  }
}

TEST(hpll, random) {
  const int kNumVertices = 50, kNumInitialEdges = 10, kNumAdditionalEdges = 50;
  const int kNumTrials = 10, kNumQueries = 100;
  const int kNumEdges = kNumInitialEdges + kNumAdditionalEdges;

  for (int trial = 0; trial < kNumTrials; ++trial) {
    vector<tuple<int, int, int>> es(kNumEdges);
    for (int i = 0; i < kNumEdges; ++i) {
      es[i] = make_tuple(rand(), rand() % kNumVertices, rand() % kNumVertices);
    }
    sort(es.begin(), es.end());

    historical_pruned_landmark_labeling a1;
    a1.construct_index(vector<tuple<int, int, int>>(es.begin(), es.begin() + kNumInitialEdges));

    ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << a1;

    historical_pruned_landmark_labeling a2;
    istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> a2;

    // Serialization test
    for (int j = 0; j < kNumQueries; ++j) {
      int v = rand() % kNumVertices;
      int w = rand() % kNumVertices;
      int t = rand();
      ASSERT_EQ(a1.query_snapshot(v, w, t), a2.query_snapshot(v, w, t));
    }

    // Online incremental update test
    for (int i = kNumInitialEdges; i < kNumEdges; ++i) {
      const auto &e = es[i];
      a1.insert_edge(get<1>(e), get<2>(e), get<0>(e));
      a2.insert_edge(get<1>(e), get<2>(e), get<0>(e));

      historical_pruned_landmark_labeling a3;
      a3.construct_index(vector<tuple<int, int, int>>(es.begin(), es.begin() + i + 1));

      for (int j = 0; j < kNumQueries; ++j) {
        int v = rand() % kNumVertices;
        int w = rand() % kNumVertices;
        int t = rand();
        ASSERT_EQ(a3.query_snapshot(v, w, t), a1.query_snapshot(v, w, t));
        ASSERT_EQ(a3.query_snapshot(v, w, t), a2.query_snapshot(v, w, t));
      }
    }
  }
}

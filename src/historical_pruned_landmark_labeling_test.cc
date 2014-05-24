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
  const int kNumTrials = 10, kNumVertices = 50, kNumEdges = 200, kNumQueries = 100000;

  for (int trial = 0; trial < kNumTrials; ++trial) {
    vector<tuple<int, int, int>> es(kNumEdges);
    for (int i = 0; i < kNumEdges; ++i) {
      es[i] = make_tuple(rand(), rand() % kNumVertices, rand() % kNumVertices);
    }

    historical_pruned_landmark_labeling a1;
    a1.construct_index(es);

    ostringstream oss;
    boost::archive::binary_oarchive oa(oss);
    oa << a1;

    historical_pruned_landmark_labeling a2;
    istringstream iss(oss.str());
    boost::archive::binary_iarchive ia(iss);
    ia >> a2;

    for (int i = 0; i < kNumQueries; ++i) {
      int v = rand() % kNumVertices;
      int w = rand() % kNumVertices;
      int t = rand();
      ASSERT_EQ(a1.query_snapshot(v, w, t), a2.query_snapshot(v, w, t));
    }
  }
}

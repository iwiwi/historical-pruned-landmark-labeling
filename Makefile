CXX = g++
CXXFLAGS = -g -Wall -Wextra -Ilib -pthread -std=c++0x -fopenmp -O3

all: bin bin/query_snapshot bin/query_change_point

bin:
	mkdir -p bin

bin/query_snapshot: sample/query_snapshot_main.cc src/historical_pruned_landmark_labeling.cc
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^
	
bin/query_change_point: sample/query_change_point_main.cc src/historical_pruned_landmark_labeling.cc
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

bin/test: src/historical_pruned_landmark_labeling.cc src/historical_pruned_landmark_labeling_test.cc lib/gtest/gtest-all.cc lib/gtest/gtest_main.cc
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: test clean

test: bin/test
	bin/test

clean:
	rm -rf bin
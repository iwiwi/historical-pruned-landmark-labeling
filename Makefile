CXX = g++
CXXFLAGS = -g -Wall -Wextra -Ilib -pthread -std=c++0x -fopenmp

all: bin bin/construct_index bin/query_snapshot bin/query_change_point

bin:
	mkdir -p bin

bin/construct_index: sample/construct_index_main.cc src/historical_pruned_landmark_labeling.cc
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^ -lboost_serialization
	
bin/query_snapshot: sample/query_snapshot_main.cc src/historical_pruned_landmark_labeling.cc
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^ -lboost_serialization
	
bin/query_change_point: sample/query_change_point_main.cc src/historical_pruned_landmark_labeling.cc
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^ -lboost_serialization

bin/test: src/historical_pruned_landmark_labeling.cc src/historical_pruned_landmark_labeling_test.cc lib/gtest/gtest-all.cc lib/gtest/gtest_main.cc
	$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_serialization

.PHONY: test clean

test: bin bin/test
	bin/test

clean:
	rm -rf bin
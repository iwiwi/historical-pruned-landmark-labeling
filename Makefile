CXX = g++
CXXFLAGS = -g -Wall -Wextra -Ilib -pthread -fopenmp -std=c++0x

all: bin bin/query

bin:
	mkdir -p bin

bin/query: sample/query_main.cc src/historical_pruned_landmark_labeling.cc
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

bin/test: src/historical_pruned_landmark_labeling.cc src/historical_pruned_landmark_labeling_test.cc lib/gtest/gtest-all.cc lib/gtest/gtest_main.cc
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: test clean

test: bin/test
	bin/test

clean:
	rm -rf bin
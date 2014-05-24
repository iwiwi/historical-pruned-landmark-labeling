Historical Pruned Landmark Labeling
===================================

When analyzing historical networks, for which timestamps
of edges are also available, in addition to the
latest snapshot, the shortest paths and distances on previous
snapshots or transition of them by time are also of interest.
By using our *historical pruned landmark labeling*,
we can answer two kinds of queries:
a snapshot query asks the distance on a specified previous snapshot,
and a change-point query asks all the moments when the distance between two vertices has changed.

Please see [our paper at WWW'14](http://www-imai.is.s.u-tokyo.ac.jp/~takiba/papers/dpll_www14.pdf)
for the detailed problem definition, algorithms and experimental results.



## Usage
Given a graph with edge time stamps, it first constructs an index. Then, using the index, it can quickly answer distance between two vertices.

### From CUI Interfaces

    $ make
    $ bin/construct_index sample/example.tsv tmp.dat
    $ bin/query_change_point tmp.dat <<< "0 2"
    0:-1	200:2	300:1

* Execute `make` to build programs.
* Execute `bin/construct_index` to construct an index.
* Execute `bin/query_snapshot` to process snapshot queries.
* Execute `bin/query_change_point` to process change-point queries.

### From Your Program

    historical_pruned_landmark_labeling hpll;
    hpll.construct_index(edge_list);
    cout << hpll.query_snapshot(0, 2, 200) << endl;

* Call `construct_index` to construct an index from a graph (an edge list or a file).
* Call `query_snapshot` to process a snapshot query.
* Call `query_change_point` to process a change-point query.

For further information, please see `historical_pruned_landmark_labeling.h`, samples and tests.

### Details

* In a graph file, each line should contain three integers describing an edge (see `samples/example.tsv`).
* Vertices should be described by integers starting from zero.
* CUI interfaces read and process queries until EOF.
* Execute `make test` to run tests (*google-gtest* is required).
* Only serialization depends on boost (boost::serialization). Note that, if you do not need to use serialization, you can use this library on systems without boost.

## Reference

* Takuya Akiba, Yoichi Iwata, and Yuichi Yoshida, **[Dynamic and Historical Shortest-Path Distance Queries on Large Networks by Pruned Landmark Labeling](http://www-imai.is.s.u-tokyo.ac.jp/~takiba/papers/dpll_www14.pdf)**.
In *WWW 2014 (Research track full paper)*.
#include "historical_pruned_landmark_labeling.h"
#include <omp.h>
#include <climits>
#include <xmmintrin.h>
#include <algorithm>
#include <cstdint>
using namespace std;

#define rep(i, n) for (int i = 0; i < (int)(n); i++)

namespace {
// Return the number of threads that would be executed in parallel regions
int get_max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

// Set the number of threads that would be executed in parallel regions
void set_num_threads(int num_threads) {
#ifdef _OPENMP
  omp_set_num_threads(num_threads);
#else
  if (num_threads != 1) {
    HPLL_CHECK(!"compile with -fopenmp");
  }
#endif
}

// Return my thread ID
int get_thread_id() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

template<typename T>
struct parallel_vector {
  parallel_vector(size_t size_limit)
      : v(get_max_threads(), vector<T>(size_limit)),
        n(get_max_threads(), 0) {}

  void push_back(const T &x) {
    int id = get_thread_id();
    v[id][n[id]++] = x;
  }

  void clear() {
    rep (i, get_max_threads()) n[i] = 0;
  }

  vector<vector<T> > v;
  vector<size_t> n;
};
}  // namespace

void historical_pruned_landmark_labeling::construct_index(const char *filename) {
  std::ifstream ifs(filename);
  HPLL_CHECK(ifs);
  construct_index(ifs);
}

void historical_pruned_landmark_labeling::construct_index(istream &ifs) {
  std::vector<std::tuple<int, int, int> > es;
  for (int t, v, w; ifs >> t >> v >> w; ) es.emplace_back(t, v, w);
  HPLL_CHECK(!ifs.bad());
  construct_index(es);
}

void historical_pruned_landmark_labeling::construct_index(const vector<tuple<int, int, int>> &es) {
  // Setup the graph
  V = 0;
  for (const auto &e : es) {
    V = max({V, get<1>(e) + 1, get<2>(e) + 1});
  }
  adj.assign(V, vector<edge_t>());
  for (const auto &e : es) {
    HPLL_CHECK(get<0>(e) >= 0);
    adj[get<1>(e)].push_back((edge_t){get<2>(e), get<0>(e)});
    adj[get<2>(e)].push_back((edge_t){get<1>(e), get<0>(e)});
  }

  // Prepare
  labels.clear();
  labels.resize(V);
  rep (v, V) labels[v].push_back(((label_entry_t){INT_MAX, 0, 0}));

  // crr_time[v] = t  <=>  can reach |v| with distance |d|   on or after time |t|
  // nxt_time[v] = t  <=>  can reach |v| with distance |d+1| on or after time |t|
  vector<int> crr_time(V, INT_MAX), nxt_time(V, INT_MAX);
  get_root_order(ord);
  parallel_vector<pair<int, label_entry_t> > pdiff_labels(V);
  parallel_vector<int> pdiff_nxt_que(V), pdiff_touched_vs(V);

  // Compute the labels
  rep (source_i, V) {
    int s = ord[source_i];
    vector<int> crr_que;
    crr_que.push_back(s);
    crr_time[s] = 0;

    size_t num_labels_added = 0;

    for (int d = 0; !crr_que.empty(); ++d) {
      int crr_que_size = crr_que.size();

#pragma omp parallel for schedule(guided, 1)
      rep (que_i, crr_que_size) {
        int v = crr_que[que_i];
        int t = crr_time[v];

        _mm_prefetch(adj[v].data(), _MM_HINT_T0);

        // Prune
        const vector<label_entry_t> &s1 = labels[s];
        const vector<label_entry_t> &s2 = labels[v];

        for (int i1 = 0, i2 = 0;;) {
          /***/if (s1[i1].t > t || s1[i1].d > d) ++i1;
          else if (s2[i2].t > t || s2[i2].d > d) ++i2;
          else if (s1[i1].v < s2[i2].v) ++i1;
          else if (s1[i1].v > s2[i2].v) ++i2;
          else {
            int v = s1[i1].v;
            if (v == INT_MAX) break;

            int q = int(s1[i1].d) + int(s2[i2].d);
            if (q <= d) goto prune;
            ++i1;
            ++i2;
          }
        }

        // Label
        pdiff_labels.push_back(make_pair(v, ((label_entry_t){source_i, d, t})));

        // Traverse
        rep (adj_i, adj[v].size()) {
          const edge_t &e = adj[v][adj_i];
          int tv = e.v;
          int tt = max(t, e.t);

          if (tt < crr_time[tv] && tt < nxt_time[tv]) {
            for (;;) {
              int prv_tt = nxt_time[tv];
              if (prv_tt <= tt) break;

              if (__sync_bool_compare_and_swap(&nxt_time[tv], prv_tt, tt)) {
                if (prv_tt == INT_MAX) {
                  pdiff_nxt_que.push_back(tv);
                  if (crr_time[tv] == INT_MAX) pdiff_touched_vs.push_back(tv);
                }
                break;
              }
            }
          }
        }

     prune:
        {}
      }
#pragma omp flush

      crr_que.clear();

      rep (i, get_max_threads()) {
        rep (j, pdiff_labels.n[i]) {
          labels[pdiff_labels.v[i][j].first].back() = pdiff_labels.v[i][j].second;
          labels[pdiff_labels.v[i][j].first].push_back(((label_entry_t){INT_MAX, 0, 0}));
          num_labels_added++;
        }
        rep (j, pdiff_nxt_que.n[i]) {
          int v = pdiff_nxt_que.v[i][j];
          crr_time[v] = nxt_time[v];
          nxt_time[v] = INT_MAX;
          crr_que.push_back(v);
        }
      }

      pdiff_labels.clear();
      pdiff_nxt_que.clear();

#pragma omp flush
    }

    rep (i, get_max_threads()) {
      rep (j, pdiff_touched_vs.n[i]) {
        int v = pdiff_touched_vs.v[i][j];
        crr_time[v] = INT_MAX;
      }
    }
    pdiff_touched_vs.clear();
  }
}

void historical_pruned_landmark_labeling::get_root_order(vector<int> &ord) {
  vector<pair<pair<int, int>, int> > deg(V);
  rep (v, V) deg[v] = make_pair(make_pair(adj[v].size(), rand()), v);
  sort(deg.begin(), deg.end());
  reverse(deg.begin(), deg.end());

  ord.resize(V);
  rep (i, V) ord[i] = deg[i].second;
}

void historical_pruned_landmark_labeling::get_label(
    int v, vector<label_entry_t> &label) {
  if (v < 0 || V <= v) label.clear();
  else label = labels[v];
}

void historical_pruned_landmark_labeling::get_index(
    vector<vector<label_entry_t>> &index) {
  index = labels;
}

double historical_pruned_landmark_labeling::get_average_label_size() {
  size_t n = 0;
  rep (v, V) n += labels[v].size() - 1;  // -1 for the sentinels
  return n / (double)V;
}

size_t historical_pruned_landmark_labeling::get_index_size() {
  size_t n = 0;
  rep (v, V) n += labels[v].size() - 1;  // -1 for the sentinels
  return n * (sizeof(int32_t) + sizeof(int8_t) + sizeof(int32_t));
}

void historical_pruned_landmark_labeling::partial_bfs(int bfs_i, int sv, int sd, int st) {
  if ((int)que.size() < V) {
    que.resize(V * 2 + 10);
    vis.resize(V * 2 + 10);
    root_label.resize(V * 2 + 10);
  }

  int r = ord[bfs_i];
  const vector<label_entry_t> &idx_r = labels[r];
  for (int i = int(idx_r.size()) - 1; i >= 0; --i) {
    if (idx_r[i].v == INT_MAX) continue;
    root_label[idx_r[i].v] = idx_r[i].d;
  }

  int que_h = 0, que_t = 0;
  que[que_t++] = make_pair(sv, sd);
  vis[sv] = true;

  while (que_h < que_t) {
    int v = que[que_h].first;
    int d = que[que_h].second;
    ++que_h;

    // Puruning test & new label
    {
      vector<label_entry_t> &idx_v = labels[v];

      int i = 0;
      for (; idx_v[i].v <= bfs_i; ++i) {
        label_entry_t &l = idx_v[i];
        if (root_label[l.v] != -1 && root_label[l.v] + l.d <= d) goto prune;
        if (l.v == bfs_i) break;
      }

      for (int j = int(idx_v.size()) - 1; j - 1 >= i; --j) idx_v[j] = idx_v[j - 1];
      idx_v.push_back(((label_entry_t){INT_MAX, 0, 0}));
      idx_v[i] = (label_entry_t){bfs_i, d, st};
    }

    rep (i, adj[v].size()) {
      int w = adj[v][i].v;
      if (!vis[w]) {
        que[que_t++] = make_pair(w, d + 1);
        vis[w] = true;
      }
    }

 prune:;
  }

  rep (i, que_t) vis[que[i].first] = false;
  rep (i, idx_r.size()) {
    if (idx_r[i].v == INT_MAX) continue;
    root_label[idx_r[i].v] = -1;
  }
}

void historical_pruned_landmark_labeling::insert_edge(int u, int v, int t) {
  if (max(u, v) >= V) {
    int x = V;
    V = max(u, v) + 1;
    for (; x < V; ++x) {
      labels.emplace_back();
      labels.back().push_back((label_entry_t){x, 0, 0});
      labels.back().push_back((label_entry_t){INT_MAX, 0, 0});
      ord.push_back(x);
    }
    adj.resize(V);
  }
  adj[u].push_back((edge_t){v, t});
  adj[v].push_back((edge_t){u, t});
  const vector<label_entry_t> &idx_u = labels[u];
  const vector<label_entry_t> &idx_v = labels[v];

  int iu = 0, iv = 0, prv_v = -1;
  for (;;) {
    label_entry_t lu = idx_u[iu];
    label_entry_t lv = idx_v[iv];

    /***/if (iu > 0 && idx_u[iu].v == prv_v) ++iu;
    else if (iv > 0 && idx_v[iv].v == prv_v) ++iv;
    else if (lu.v < lv.v) {  // u -> v
      partial_bfs(prv_v = lu.v, v, lu.d + 1, t);
      ++iu;
    } else if (lu.v > lv.v) {  // v -> u
      partial_bfs(prv_v = lv.v, u, lv.d + 1, t);
      ++iv;
    } else {  // u <-> v
      if (lu.v == INT_MAX) break;  // sentinel
      if (lu.d + 1 < lv.d) partial_bfs(lu.v, v, lu.d + 1, t);
      if (lv.d + 1 < lu.d) partial_bfs(lv.v, u, lv.d + 1, t);
      prv_v = lu.v;
      ++iu;
      ++iv;
    }
  }
}

int historical_pruned_landmark_labeling::query_snapshot(int v, int w, int t) {
  if (v < 0 || w < 0 || V <= v || V <= w) return -1;

  const vector<label_entry_t> &s1 = labels[v];
  const vector<label_entry_t> &s2 = labels[w];
  int d = INT_MAX;

  size_t i1 = 0, i2 = 0;
  while (i1 < s1.size() && i2 < s2.size()) {
    /***/if (s1[i1].t > t) ++i1;
    else if (s2[i2].t > t) ++i2;
    else if (s1[i1].v < s2[i2].v) ++i1;
    else if (s1[i1].v > s2[i2].v) ++i2;
    else {
      int v = s1[i1].v;
      if (v == INT_MAX) break;

      d = min(d, int(s1[i1].d) + int(s2[i2].d));
      for (++i1; s1[i1].v == v; ++i1);
      for (++i2; s2[i2].v == v; ++i2);
    }
  }
  return d == INT_MAX ? -1 : d;
}

void historical_pruned_landmark_labeling::query_change_points(int v, int w,
                                                              vector<pair<int, int>> &cp) {
  cp.clear();
  cp.emplace_back(0, INT_MAX);
  if (v < 0 || w < 0 || V <= v || V <= w) return;

  vector<label_entry_t> &s1 = labels[v];
  vector<label_entry_t> &s2 = labels[w];
  if (s1.back().v != INT_MAX) s1.push_back(((label_entry_t){INT_MAX, 0, 0}));
  if (s2.back().v != INT_MAX) s2.push_back(((label_entry_t){INT_MAX, 0, 0}));

  size_t i1 = 0, i2 = 0;
  for (;;) {
    /***/if (s1[i1].v < s2[i2].v) ++i1;
    else if (s1[i1].v > s2[i2].v) ++i2;
    else {
      int x = s1[i1].v;
      if (x == INT_MAX) break;  // Sentinel

      while (s1[i1].v == x && s2[i2].v == x) {
        cp.emplace_back(max(s1[i1].t, s2[i2].t), s1[i1].d + s2[i2].d);
        if (s2[i2 + 1].v != x || (s1[i1 + 1].v == x && s1[i1].t > s2[i2].t)) ++i1;
        else                                                                 ++i2;
      }
    }
  }

  sort(cp.begin(), cp.end());
  int j = 1;
  for (int i = 1; i < (int)cp.size(); ++i) {
    if (cp[j - 1].second > cp[i].second) cp[j++] = cp[i];
  }
  cp.resize(j);
  for (auto &p : cp) {
    if (p.second == INT_MAX) p.second = -1;
  }
}

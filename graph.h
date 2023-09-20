#include <iostream>
#include <string>
#include <vector>
#include <queue>

#include "pvector.h"
#include "sliding_queue.h"
#include "bitmap.h"

#include "util.h"
#include "cpp-btree-1.0.1/btree_set.h"

using namespace std;

typedef pair<long int, long int> edge;
typedef btree::btree_set<long int> neighbors;

class Graph {
private:
	long int N; // number of nodes/vertices
	long int E; // number of edges
	vector<neighbors> Parents;
	vector<neighbors> Neighbors;
public:
	Graph(int N) {
		this->N = N; // number of nodes/vertices
		this->E = 0;
		for (long int i = 0; i < N; i++) {
			Neighbors.push_back({});
			Parents.push_back({});
		}
	}
	void addEdge (edge edge) {
		if (this->Neighbors[edge.first].find(edge.second)
			== this->Neighbors[edge.first].end()) {
			this->Neighbors[edge.first].insert(edge.second);
			this->Parents[edge.second].insert(edge.first);
			this->E++;
		}
	}
	void bfs (long int source);
	void printNeighbors () {
		for (auto &i : this->Neighbors) {
			for (auto &j : i)
				cout << j << " ";
			cout << endl;
		}
	}
	pvector<long int> bfs_gap (long int source, int alpha, int beta);
	void QueueToBitmap(const SlidingQueue<long int> &q, Bitmap &bm);
	void BitmapToQueue(const Bitmap &bm, SlidingQueue<long int> &q); 
	long int BUStep(pvector<long int> &parent, Bitmap &front, Bitmap &next);
	long int TDStep(pvector<long int> &parent, SlidingQueue<long int> &q);
	bool BFSVerifier(long int source, const pvector<long int> &parent);
};

void Graph::bfs (long int source) {
	// printNeighbors();
	queue<long int> q;
	vector<bool> visited(this->N, false);
	q.push(source); visited[source] = true;
	cout << "BFS: ";
	while (!q.empty()) {
		cout << q.front() << " ";
		for (auto &i : Neighbors[q.front()]) {
			if (!visited[i]) {
				q.push(i);
				visited[i] = true;
			}
		}
		q.pop();
	}
	cout << endl;
}

void Graph::QueueToBitmap(const SlidingQueue<long int> &q, Bitmap &bm) {
	#pragma omp parallel for
	for (auto q_iter = q.begin(); q_iter < q.end(); q_iter++) {
		long int u = *q_iter;
		bm.set_bit_atomic(u);
	}
}

void Graph::BitmapToQueue(const Bitmap &bm,
                   SlidingQueue<long int> &q) {
  #pragma omp parallel
  {
    QueueBuffer<long int> lqueue(q);
    #pragma omp for nowait
    for (long int n=0; n < this->N; n++)
      if (bm.get_bit(n))
        lqueue.push_back(n);
    lqueue.flush();
  }
  q.slide_window();
}

long int Graph::BUStep(pvector<long int> &parent, Bitmap &front, Bitmap &next) {
  long int awake_count = 0;
  next.reset();
  #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
  for (long int u=0; u < this->N; u++) {
    if (parent[u] < 0) {
      for (long int v : this->Parents[u]) {
        if (front.get_bit(v)) {
          parent[u] = v;
          awake_count++;
          next.set_bit(u);
          break;
        }
      }
    }
  }
  return awake_count;
}

long int Graph::TDStep(pvector<long int> &parent, SlidingQueue<long int> &q) {
  long int scout_count = 0;
  #pragma omp parallel
  {
    QueueBuffer<long int> lqueue(q);
    #pragma omp for reduction(+ : scout_count) nowait
    for (auto q_iter = q.begin(); q_iter < q.end(); q_iter++) {
      long int u = *q_iter;
      for (long int v : this->Neighbors[u]) {
        long int curr_val = parent[v];
        if (curr_val < 0) {
          if (compare_and_swap(parent[v], curr_val, u)) {
            lqueue.push_back(v);
            scout_count += -curr_val;
          }
        }
      }
    }
    lqueue.flush();
  }
  return scout_count;
}

pvector<long int> Graph::bfs_gap (long int source, 
	int alpha = 15, 
	int beta = 18) {
	vector<double> timestamps1 = {};
	timestamps1.push_back(get_wall_time());
	pvector<long int> parent(this->N);
	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
	    parent[i] = this->Neighbors[i].size() != 0 
			? -this->Neighbors[i].size() 
			: -1;
	}
	parent[source] = source;
	SlidingQueue<long int> q(this->N);
	q.push_back(source);
	q.slide_window();
	Bitmap curr(this->N);
	curr.reset();
    Bitmap front(this->N);
	front.reset();
	long int edges_to_check = this->E;
	long int scout_count = this->Neighbors[source].size();
	while (!q.empty()) {
		if (scout_count > edges_to_check / alpha) {
			long int awake_count, old_awake_count;
			QueueToBitmap(q, front);
			awake_count = q.size();
			q.slide_window();
			do {
				old_awake_count = awake_count;
				awake_count = BUStep(parent, front, curr);
				front.swap(curr);
			} while ((awake_count >= old_awake_count) 
				|| (awake_count > this->N / beta));
			BitmapToQueue(front, q);
			scout_count = 1;
		} else {
			edges_to_check -= scout_count;
			scout_count = TDStep(parent, q);
			q.slide_window();
		}
	}
	timestamps1.push_back(get_wall_time());
	cout << "Real wall time for bfs_gap after the initial round: "
	 << timestamps1.back() - timestamps1.end()[-2] << endl;
	#pragma omp parallel for
	for (long int n = 0; n < this->N; n++)
		if (parent[n] < -1) 
	  		parent[n] = -1;
	return parent;
}

// BFS verifier does a serial BFS from same source and asserts:
// - parent[source] = source
// - parent[v] = u  =>  depth[v] = depth[u] + 1 (except for source)
// - parent[v] = u  => there is edge from u to v
// - all vertices reachable from source have a parent
bool Graph::BFSVerifier(long int source,
                 const pvector<long int> &parent) {
  pvector<int> depth(this->N, -1);
  depth[source] = 0;
  vector<long int> to_visit;
  to_visit.reserve(this->N);
  to_visit.push_back(source);
  for (auto it = to_visit.begin(); it != to_visit.end(); it++) {
    long int u = *it;
    for (long int v : this->Neighbors[u]) {
      if (depth[v] == -1) {
        depth[v] = depth[u] + 1;
        to_visit.push_back(v);
      }
    }
  }
  for (long int u = 0; u < this->N; u++) {
    if ((depth[u] != -1) && (parent[u] != -1)) {
      if (u == source) {
        if (!((parent[u] == u) && (depth[u] == 0))) {
          cout << "Source wrong" << endl;
          return false;
        }
        continue;
      }
      bool parent_found = false;
      for (long int v : this->Parents[u]) {
        if (v == parent[u]) {
          if (depth[v] != depth[u] - 1) {
            cout << "Wrong depths for " << u << " & " << v << endl;
            return false;
          }
          parent_found = true;
          break;
        }
      }
      if (!parent_found) {
        cout << "Couldn't find edge from " << parent[u] << " to " << u << endl;
        return false;
      }
    } else if (depth[u] != parent[u]) {
      cout << "Reachability mismatch" << endl;
      return false;
    }
  }
  return true;
}

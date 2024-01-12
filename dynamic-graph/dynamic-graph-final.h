#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <unordered_map>

#include "../src-gap/pvector.h"
#include "../src-gap/sliding_queue.h"
#include "../src-gap/bitmap.h"

#include "../util.h"
#include "../tlx/btree_set.hpp"

#include <type_traits>
using namespace std;

typedef pair<long int, long int> edge;
typedef tlx::btree_set<long int> neighbors;

class DynamicGraph {
private:
	long int N; // number of nodes/vertices
	long int E; // number of edges
	vector<neighbors> Parents, Neighbors;
	pvector<long int> depths;
	vector<edge> new_edges;
public:
	DynamicGraph (int N) {
		this->N = N; // number of nodes/vertices
		this->E = 0;
		for (long int i = 0; i < this->N; i++) {
			neighbors n;
      this->Neighbors.push_back(n);
			this->Parents.push_back(n);
			this->depths.reserve(this->N);
			this->depths.fill(-1);
		}
	}
	void buildGraph (vector<edge>& edge_list, 
		vector<edge>& parents_list, 
		long int start_line, 
		long int num_lines);

	void mergeEdges (vector<edge>& edge_list, 
		vector<edge>& parents_list, 
		const pvector<long int>& parents, 
		long int start_line, 
		long int num_lines);


	pvector<long int> bfs_gap (long int source);
	void bfs_gap_incremental_use_depth (long int source, 
		pvector<long int>& parent);
	void bfs_gap_incremental_use_parent (long int source, 
		pvector<long int>& parent);
	pvector<long int> traditional_bfs (long int source);
	long int get_depth (long int curr_val, 
		long int source, 
		const pvector<long int>& parent);

	bool BFSVerifier(long int source, 
		const pvector<long int> &parent);

	long int TDStep (pvector<long int> &parent, 
		SlidingQueue<long int> &q);
	long int TDStepUseParent (pvector<long int> &parent, 
		SlidingQueue<long int> &q, long int source);
	long int TDStepUseDepth (pvector<long int> &parent, 
		SlidingQueue<long int> &q, long int source);
};

void DynamicGraph::buildGraph (vector<edge>& edge_list, 
	vector<edge>& parents_list, 
	long int start_line,
	long int num_lines) {
	this->E = 0;
	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
		this->Neighbors[i].clear();
		this->Parents[i].clear();
	}

	vector<vector<long int> > tmp_p, tmp_e;
	vector<long int> tmp = {};
	for (long int i = 0; i < this->N; i++) {
		tmp_p.push_back(tmp);
		tmp_e.push_back(tmp);
	}

	for (long int i = start_line; i < num_lines; i++) {
		tmp_p[parents_list[i].first].push_back(parents_list[i].second);
		tmp_e[edge_list[i].first].push_back(edge_list[i].second);
	}

	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
		std::sort (tmp_p[i].begin(), tmp_p[i].end());
		std::sort (tmp_e[i].begin(), tmp_e[i].end());
	}

	double timer_start = 0, timer_end = 0;
	timer_start = get_wall_time();
	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
		this->Neighbors[i].bulk_load(tmp_e[i].begin(), tmp_e[i].end());
		this->Parents[i].bulk_load(tmp_p[i].begin(), tmp_p[i].end());
	}
	timer_end = get_wall_time();
	cout << "Wall time to update edges: " 
		<< (timer_end - timer_start) << endl;
	this->E = num_lines;
}

void DynamicGraph::mergeEdges (vector<edge>& edge_list, 
	vector<edge>& parents_list,
	const pvector<long int>& parent, 
	long int start_line,
	long int num_lines) {
	this->new_edges.clear();
	vector<vector<long int> > tmp_p, tmp_e;
	vector<long int> tmp = {};
	for (long int i = 0; i < this->N; i++) {
		tmp_p.push_back(tmp);
		tmp_e.push_back(tmp);
	}

	for (long int i = start_line; i < num_lines; i++) {
		tmp_p[parents_list[i].first].push_back(parents_list[i].second);
		tmp_e[edge_list[i].first].push_back(edge_list[i].second);
		this->new_edges.push_back(edge_list[i]);
	}

	double timer_start = 0, timer_end = 0;
	timer_start = get_wall_time();
	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
		for (auto &j : tmp_e[i]) {
			this->Neighbors[i].insert(j);
		}
		for (auto &j : tmp_p[i]) {
			this->Parents[i].insert(j);
		}
	}
	timer_end = get_wall_time();
	cout << "Wall time to update/merge edges: " 
		<< (timer_end - timer_start) << endl;
	this->E = num_lines;
}

long int DynamicGraph::TDStep (pvector<long int> &parent, 
	SlidingQueue<long int> &q) {
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

pvector<long int> DynamicGraph::bfs_gap (long int source) {
	pvector<long int> parent(this->N);
	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
	    parent[i] = this->Neighbors[i].size() != 0 
			? -this->Neighbors[i].size() 
			: -1;
	}
	vector<double> timestamps = {};
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
	timestamps.push_back(get_wall_time());
	while (!q.empty()) {
		edges_to_check -= scout_count;
		scout_count = TDStep(parent, q);
		q.slide_window();
	}
	timestamps.push_back(get_wall_time());
	cout << "Real wall time for bfs_gap: "
	 << timestamps.back() - timestamps.end()[-2] << endl;
	#pragma omp parallel for
	for (long int n = 0; n < this->N; n++)
		if (parent[n] < -1) 
	  		parent[n] = -1;
	cout << "BFSVerifier output: " 
		 << boolalpha << BFSVerifier (source, parent) << endl;
	return parent;
}

bool DynamicGraph::BFSVerifier(long int source, 
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

long int DynamicGraph::get_depth (long int curr_val, 
	long int source, 
	const pvector<long int>& parent) {
	long int depth = 0;
	while (curr_val != source) {
		if (curr_val < 0)
			return this->N;
		curr_val = parent[curr_val];
		depth++;
	}
	return depth;
}

/*long int DynamicGraph::TDStepUseParent 
	(pvector<long int> &parent, 
	SlidingQueue<long int> &q,
	long int source) {
	long int scout_count = 0;
	#pragma omp parallel
	{
	QueueBuffer<long int> lqueue(q);
	#pragma omp for reduction(+ : scout_count) nowait
	for (auto q_iter = q.begin(); q_iter < q.end(); q_iter++) {
		long int u = *q_iter;
		for (long int v : this->Neighbors[u]) {
			long int curr_val = parent[v];
			if (curr_val < 0 ||
				(curr_val > 0 &&  
				(get_depth(curr_val, source, parent) > 
				get_depth(u, source, parent)))) {
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

void DynamicGraph::bfs_gap_incremental_use_parent 
	(long int source, 
	pvector<long int>& parent) {
	SlidingQueue<long int> q(this->N);
	double timer_start, timer_end;
	timer_start = get_wall_time();
	#pragma omp parallel for
	for (auto &j : this->new_edges) {
		if (parent[j.first] < 0)
      continue;
		else if (parent[j.first] > 0 
			&& parent[j.second] < 0) {
			q.push_back(j.second);
			q.slide_window();
			parent[j.second] = j.first;
			// compare_and_swap(parent[j.second], -1, j.first);
			Bitmap curr(this->N);
			curr.reset();
			Bitmap front(this->N);
			front.reset();
			long int edges_to_check = this->E;
			long int scout_count = this->Neighbors[*q.begin()].size();
			while (!q.empty()) {
				edges_to_check -= scout_count;
				scout_count = TDStepUseParent (parent, q, source);
				q.slide_window();
			}
		}
		else if (parent[j.first] > 0 
			&& parent[j.second] > 0 
			&& get_depth(j.first, source, parent) < 
			get_depth(parent[j.second], source, parent)) {
			compare_and_swap(parent[j.second], parent[j.second], j.first);
		}
	}
	timer_end = get_wall_time();
	cout << "Real wall time to update bfs_gap: "
	 << timer_end - timer_start << endl;
	#pragma omp parallel for
	for (long int n = 0; n < this->N; n++)
		if (parent[n] < -1) 
	  		parent[n] = -1;
	cout << "BFSVerifier output: " 
	 	 << boolalpha << BFSVerifier (source, parent) << endl;
}*/

long int DynamicGraph::TDStepUseDepth 
	(pvector<long int> &parent, 
	SlidingQueue<long int> &q,
	long int source) {
	long int scout_count = 0;
	#pragma omp parallel
	{
	QueueBuffer<long int> lqueue(q);
	#pragma omp for reduction(+ : scout_count) nowait
	for (auto q_iter = q.begin(); q_iter < q.end(); q_iter++) {
		long int u = *q_iter;
		for (long int v : this->Neighbors[u]) {
			long int curr_val = parent[v];
			if (curr_val < 0 ||
				(curr_val > 0 &&
				(this->depths[v] > this->depths[u] + 1))) {
				if (compare_and_swap(parent[v], curr_val, u)) {
					this->depths[v] = this->depths[u] + 1;
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

/*void DynamicGraph::bfs_gap_incremental_use_depth
	(long int source, 
	pvector<long int>& parent) {
	SlidingQueue<long int> q(this->N);
	double timer_start, timer_end;
	timer_start = get_wall_time();
	#pragma omp parallel for
	for (auto &j : this->new_edges) {
		if (parent[j.first] < 0)
      continue;
		else if ((parent[j.first] > 0 
			&& parent[j.second] < 0) ||
			(parent[j.first] > 0 
			&& parent[j.second] > 0 
			&& this->depths[j.first] < this->depths[parent[j.second]])) {
			q.push_back(j.second);
			q.slide_window();
			parent[j.second] = j.first;
			this->depths[j.second] = this->depths[j.first] + 1;
			Bitmap curr(this->N);
			curr.reset();
			Bitmap front(this->N);
			front.reset();
			long int edges_to_check = this->E;
			long int scout_count = this->Neighbors[*q.begin()].size();
			while (!q.empty()) {
				edges_to_check -= scout_count;
				scout_count = TDStepUseDepth (parent, q, source);
				q.slide_window();
			}
		}
	}
	timer_end = get_wall_time();
	cout << "Real wall time to update bfs_gap: "
	 << timer_end - timer_start << endl;

	#pragma omp parallel for
	for (long int n = 0; n < this->N; n++)
		if (parent[n] < 0) 
	  		parent[n] = -1;

	bool bfs_verifier = BFSVerifier (source, parent);
	cout << "BFSVerifier output: " << boolalpha << bfs_verifier << endl;

	if (bfs_verifier != true) {
		pvector<long int> parent_gap_bfs = bfs_gap (source);
		for (int i = 0; i < this->N; i++) {
			long int depth_of_incremental_bfs 
				= get_depth(parent[i], source, parent);
			long int depth_of_static_bfs 
				= get_depth(parent_gap_bfs[i], source, parent_gap_bfs);		
	 		if (depth_of_incremental_bfs != depth_of_static_bfs) {
	 			cout << i << " " << parent[i] << " " 
	 				<< depth_of_incremental_bfs << " " 
	 				<< parent_gap_bfs[i] << " " 
	 				<< depth_of_static_bfs << endl;
	 		}
	 	}
 	}
}*/

void DynamicGraph::bfs_gap_incremental_use_depth
	(long int source, 
	pvector<long int>& parent) {
	double timer_start, timer_end;
	timer_start = get_wall_time();
	#pragma omp parallel for
	for (auto &j : this->new_edges) {
		if (parent[j.first] < 0)
      continue;
		else if ((parent[j.first] > 0 
			&& parent[j.second] < 0) ||
			(parent[j.first] > 0 
			&& parent[j.second] > 0 
			&& this->depths[j.first] < this->depths[parent[j.second]])) {
			SlidingQueue<long int> q(this->N);
			q.push_back(j.second);
			q.slide_window();
			parent[j.second] = j.first;
			this->depths[j.second] = this->depths[j.first] + 1;
			Bitmap curr(this->N);
			curr.reset();
			Bitmap front(this->N);
			front.reset();
			long int edges_to_check = this->E;
			long int scout_count = this->Neighbors[*q.begin()].size();
			while (!q.empty()) {
				edges_to_check -= scout_count;
				scout_count = TDStepUseDepth (parent, q, source);
				q.slide_window();
			}
		}
	}
	timer_end = get_wall_time();
	cout << "Real wall time to update bfs_gap: "
	 << timer_end - timer_start << endl;

	#pragma omp parallel for
	for (long int n = 0; n < this->N; n++)
		if (parent[n] < 0) 
	  		parent[n] = -1;

	bool bfs_verifier = BFSVerifier (source, parent);
	cout << "BFSVerifier output: " << boolalpha << bfs_verifier << endl;

	if (bfs_verifier != true) {
		pvector<long int> parent_gap_bfs = bfs_gap (source);
		for (int i = 0; i < this->N; i++) {
			long int depth_of_incremental_bfs 
				= get_depth(parent[i], source, parent);
			long int depth_of_static_bfs 
				= get_depth(parent_gap_bfs[i], source, parent_gap_bfs);		
	 		if (depth_of_incremental_bfs != depth_of_static_bfs) {
	 			cout << i << " " << parent[i] << " " 
	 				<< depth_of_incremental_bfs << " " 
	 				<< parent_gap_bfs[i] << " " 
	 				<< depth_of_static_bfs << endl;
	 		}
	 	}
 	}
}



pvector<long int> DynamicGraph::traditional_bfs (long int source) {
	std::queue<long int> q;
	pvector<long int> parent(this->N);
	for (long int i = 0; i < this->N; i++) {
		parent[i] = -1;
	}	
	parent[source] = source; 
	this->depths[source] = 0;
	q.push(source);
	double timer_start, timer_end;
	timer_start = get_wall_time();
	while (!q.empty()) {
		for (auto &j : this->Neighbors[q.front()]) {
			if (parent[j] < 0) {
				parent[j] = q.front(); 
				this->depths[j] = this->depths[q.front()] + 1;
				q.push(j);
			}
		}
		q.pop();
	}
	timer_end = get_wall_time();
	cout << "Real wall time for traditional_bfs: "
	 << timer_end - timer_start << endl;
	
	cout << "BFSVerifier output: " 
	 	 << boolalpha << BFSVerifier (source, parent) << endl;
	return parent;
}


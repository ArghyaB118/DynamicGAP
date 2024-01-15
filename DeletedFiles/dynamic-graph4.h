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
bool verbose = false;


class DynGraph {
private:
	long int N; // number of nodes/vertices
	long int E; // number of edges
	vector<neighbors> Parents, Neighbors;
	vector<vector<long int> > p, e;
	unordered_map<long int, long int> depths;
	pvector<long int> depths_bfs_tree;
	vector<vector<long int> > bfs_vector;
	vector<edge> potential_sources;
public:
	DynGraph (int N) {
		this->N = N; // number of nodes/vertices
		this->E = 0;
		for (long int i = 0; i < this->N; i++) {
			neighbors n;
      		this->Neighbors.push_back(n);
			this->Parents.push_back(n);
			vector<long int> tmp = {};
			this->p.push_back(tmp);
			this->e.push_back(tmp);
			this->bfs_vector.push_back(tmp);
			this->depths_bfs_tree.push_back(-1);
		}
	}
	void bfsFromScratch (long int source);
	void bfs (long int source);
	void bfs2 (long int source);
	
	void printNeighbors () {
		for (auto &i : this->Neighbors) {
			for (auto &j : i)
				cout << j << " ";
			cout << endl;
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

	void bfsDelta (vector<edge>& edge_list, vector<edge>& parents_list, long int start_line, long int num_lines);	
	void formBFS ();

	long int TDStep (pvector<long int> &parent, SlidingQueue<long int> &q);
	long int TDStepMod (pvector<long int> &parent, SlidingQueue<long int> &q, long int source);
	
	pvector<long int> bfs_gap (long int source);
	void printBFS (pvector<long int> &parent);
	bool BFSVerifier(long int source,
                 const pvector<long int> &parent);

	long int get_depth (long int curr_val, long int source, const pvector<long int>& parent);

	void bfs_gap_incremental (long int source, pvector<long int>& parent);
	void bfs_gap_incremental2 (long int source, pvector<long int>& parent);
	void bfs_gap_incremental3 (long int source, pvector<long int>& parent);

	pvector<long int> bfs_scratch (long int source);

	void get_potential_sources (vector<edge>& edge_list, 
		long int start_line,
		long int num_lines);
};

void DynGraph::bfsDelta (vector<edge>& edge_list, 
	vector<edge>& parents_list, 
	long int start_line,
	long int num_lines) {
	for (long int i = start_line; i < num_lines; i++) {
		if (this->depths.find(edge_list[i].first) 
			!= this->depths.end()) {
			if (this->depths.find(edge_list[i].second) 
				!= this->depths.end() &&
				this->depths[edge_list[i].second] 
				<= this->depths[edge_list[i].first] + 1)
				continue;
			else {
				this->depths[edge_list[i].second] 
					= this->depths[edge_list[i].first] + 1;
				queue<long int> q;
				q.push(edge_list[i].second);
				while (!q.empty()) {
					for (auto &j : this->Neighbors[q.front()]) {
						if (this->depths[q.front()] + 1 
							< this->depths[j]) {
							q.push(j); 
							this->depths[j] = this->depths[q.front()] + 1;
						}
					}
					q.pop();
				}
			}				
		}
	}
}

void DynGraph::buildGraph (vector<edge>& edge_list, 
	vector<edge>& parents_list, 
	long int start_line,
	long int num_lines) {
	this->E = 0;
	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
		this->Neighbors[i].clear();
		this->Parents[i].clear();
	}
	for (long int i = start_line; i < num_lines; i++) {
		this->p[parents_list[i].first].push_back(parents_list[i].second);
		this->e[edge_list[i].first].push_back(edge_list[i].second);
	}
	double timer_start = 0, timer_end = 0;
	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
		std::sort (p[i].begin(), p[i].end());
		std::sort (e[i].begin(), e[i].end());
	}
	timer_start = get_wall_time();
	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
		this->Neighbors[i].bulk_load(e[i].begin(), e[i].end());
		this->Parents[i].bulk_load(p[i].begin(), p[i].end());
	}
	timer_end = get_wall_time();
	cout << "Wall time to update edges: " 
		<< (timer_end - timer_start) << endl;
	this->E = num_lines;

	/*if (start_line != 0) {
		timer_start = get_wall_time();
		bfsDelta (edge_list, parents_list, start_line, num_lines);
		timer_end = get_wall_time();
		cout << "Wall time to preprocess: " 
			<< (timer_end - timer_start) << endl;
	}
	formBFS();*/
}

void DynGraph::get_potential_sources (vector<edge>& edge_list,  
	long int start_line,
	long int num_lines) {
	this->potential_sources.clear();
	for (long int i = start_line; i < num_lines; i++) {
		this->potential_sources.push_back(edge_list[i]);
	}
}

void DynGraph::mergeEdges (vector<edge>& edge_list, 
	vector<edge>& parents_list,
	const pvector<long int>& parent, 
	long int start_line,
	long int num_lines) {
	this->potential_sources.clear();
	vector<vector<long int> > tmp_p, tmp_e;
	vector<long int> tmp = {};
	for (long int i = 0; i < this->N; i++) {
		tmp_p.push_back(tmp);
		tmp_e.push_back(tmp);
	}
	for (long int i = start_line; i < num_lines; i++) {
		this->p[parents_list[i].first].push_back(parents_list[i].second);
		this->e[edge_list[i].first].push_back(edge_list[i].second);
		tmp_p[parents_list[i].first].push_back(parents_list[i].second);
		tmp_e[edge_list[i].first].push_back(edge_list[i].second);
		if (parent[edge_list[i].first] > 0)
			this->potential_sources.push_back(edge_list[i]);
	}
	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
		sort (this->p[i].begin(), this->p[i].end());
		sort (this->e[i].begin(), this->e[i].end());
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

void DynGraph::bfsFromScratch (long int source) {
	queue<long int> q;
	q.push(source);
	this->depths[source] = 0;

	if (verbose) { cout << "BFS: "; }
	while (!q.empty()) {
		if (verbose) { cout << q.front() << " "; }
		for (auto &i : this->Neighbors[q.front()]) {
			if (this->depths.find(i) == this->depths.end()) {
				q.push(i);
				this->depths[i] = this->depths[q.front()] + 1;
			}
		}
		q.pop();
	}
	if (verbose) { cout << endl; }
}

void DynGraph::bfs (long int source) {
	queue<long int> q;
	q.push(source);
	this->depths[source] = 0;

	while (!q.empty()) {
		long int qSize = q.size();
		#pragma omp parallel for
		for (long int k = 0; k < qSize; k++) {
			long int tmp;
			#pragma omp critical
	    	{
	      		tmp = q.front();
	      		q.pop();
				for (auto &i : this->Neighbors[tmp]) {
					if (this->depths.find(i) == this->depths.end()) {
						q.push(i);
						this->depths[i] = this->depths[tmp] + 1;
					}
				}
			}
		}
	}
}

void DynGraph::bfs2 (long int source) {
	queue<long int> q;
	q.push(source);
	this->depths[source] = 0;

	while (!q.empty()) {
		long int qSize = q.size();
		for (long int k = 0; k < qSize; k++) {
			long int tmp;
			{
	      		tmp = q.front();
	      		q.pop();
	      	}
			for (auto &i : this->Neighbors[tmp]) {
				if (this->depths.find(i) == this->depths.end()) {
					q.push(i);
					// this->depths[i] = this->depths[tmp] + 1;
				}
			}
		}
	}
}


void DynGraph::formBFS () {
	double timer_start = 0, timer_end = 0;
	timer_start = get_wall_time();
	for (auto &i : this->depths) {
		bfs_vector[i.second].push_back(i.first);
	}
	timer_end = get_wall_time();
	cout << "Wall time to form BFS: " 
		<< (timer_end - timer_start) << endl;
}

long int DynGraph::TDStep (pvector<long int> &parent, 
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

pvector<long int> DynGraph::bfs_gap (long int source) {
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
	double timer_start, timer_end;
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
	printBFS (parent);
	cout << "BFSVerifier output: " 
		 << boolalpha << BFSVerifier (source, parent) << endl;
	return parent;
}

void DynGraph::printBFS (pvector<long int> &parent) {
	for (auto &i : parent) {
		cout << i << " ";
	}
	cout << endl;
}

bool DynGraph::BFSVerifier(long int source,
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

long int DynGraph::get_depth (long int curr_val, 
	long int source, 
	const pvector<long int>& parent) {
	long int depth = 0;
	while (curr_val != source) {
		curr_val = parent[curr_val];
		if (curr_val < 0)
			return this->N;
		depth++;
	}
	return depth;
}

long int DynGraph::TDStepMod (pvector<long int> &parent, 
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

void DynGraph::bfs_gap_incremental 
	(long int source, 
	pvector<long int>& parent) {
	SlidingQueue<long int> q(this->N);
	double timer_start, timer_end;
	timer_start = get_wall_time();
	for (auto &j : this->potential_sources) {
		if (parent[j.first] < 0)
      continue;
		else if (parent[j.first] > 0 
			&& parent[j.second] < 0) {
			cout << "a" << q.size() << endl;
			q.push_back(j.second);
			q.slide_window();
			cout << "b" << q.size() << endl;
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
				scout_count = TDStepMod (parent, q, source);
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
	printBFS (parent);
	cout << "BFSVerifier output: " 
	 	 << boolalpha << BFSVerifier (source, parent) << endl;
}

void DynGraph::bfs_gap_incremental2 
	(long int source, 
	pvector<long int>& parent) {
	SlidingQueue<long int> q(this->N);
	for (auto &j : this->potential_sources) {
		if (parent[j.first] < 0)
      continue;
		else if (parent[j.first] > 0 
			&& parent[j.second] < 0) {
			q.push_back(j.second);
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
				scout_count = TDStepMod (parent, q, source);
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
	cout << q.size() << endl;
	q.slide_window();
	cout << q.size() << endl;
	/*Bitmap curr(this->N);
	curr.reset();
	Bitmap front(this->N);
	front.reset();
	long int edges_to_check = this->E;
	long int scout_count = this->Neighbors[*q.begin()].size();*/
	double timer_start, timer_end;
	timer_start = get_wall_time();
	/*while (!q.empty()) {
		edges_to_check -= scout_count;
		scout_count = TDStepMod (parent, q, source);
		q.slide_window();
	}*/
	timer_end = get_wall_time();
	cout << "Real wall time to update bfs_gap: "
	 << timer_start - timer_end << endl;
	#pragma omp parallel for
	for (long int n = 0; n < this->N; n++)
		if (parent[n] < -1) 
	  		parent[n] = -1;
	printBFS (parent);
	cout << "BFSVerifier output: " 
	 	 << boolalpha << BFSVerifier (source, parent) << endl;
}

void DynGraph::bfs_gap_incremental3 
	(long int source, 
	pvector<long int>& parent) {
	std::queue<long int> q;
	for (auto &j : this->potential_sources) {
		if (parent[j.first] < 0)
      continue;
		else if (parent[j.first] > 0 
			&& parent[j.second] < 0) {
			q.push(j.second);
			parent[j.second] = j.first;
			while (!q.empty()) {
				for (auto &k : Neighbors[q.front()]) {
					if (parent[k] < 0) {
						parent[k] = q.front();
						q.push(k);
					}
					else if (get_depth(q.front(), source, parent) 
						< get_depth(parent[k], source, parent)) {
						parent[k] = q.front();
					}
				}
				q.pop();
			}
		}
		else if (parent[j.first] > 0 
			&& parent[j.second] > 0 
			&& get_depth(j.first, source, parent) < 
			get_depth(parent[j.second], source, parent)) {
			parent[j.second] = j.first;
		}
	}
	#pragma omp parallel for
	for (long int n = 0; n < this->N; n++)
		if (parent[n] < 0) 
	  		parent[n] = -1;
	// printBFS (parent);
	cout << "BFSVerifier output: " 
	 	 << boolalpha << BFSVerifier (source, parent) << endl;
}

void DynGraph::bfs_gap_incremental4 
	(long int source, 
	pvector<long int>& parent) {
	std::queue<long int> q;
	for (auto &j : this->potential_sources) {
		if (parent[j.first] < 0)
      continue;
		else if (parent[j.first] > 0 
			&& parent[j.second] < 0) {
			q.push(j.second);
			parent[j.second] = j.first;
			while (!q.empty()) {
				for (auto &k : Neighbors[q.front()]) {
					if (parent[k] < 0) {
						parent[k] = q.front(); 
						this->depths_bfs_tree[k] = this->depths_bfs_tree[q.front()] + 1;
						q.push(k);
					}
					else if (this->depths_bfs_tree[q.front()] 
						< this->depths_bfs_tree[parent[k]]) {
						parent[k] = q.front(); 
						this->depths_bfs_tree[parent[k]] = this->depths_bfs_tree[q.front()];
					}
				}
				q.pop();
			}
		}
		else if (parent[j.first] > 0 
			&& parent[j.second] > 0 
			&& get_depth(j.first, source, parent) < 
			get_depth(parent[j.second], source, parent)) {
			parent[j.second] = j.first;
		}
	}
	#pragma omp parallel for
	for (long int n = 0; n < this->N; n++)
		if (parent[n] < 0) 
	  		parent[n] = -1;
	// printBFS (parent);
	cout << "BFSVerifier output: " 
	 	 << boolalpha << BFSVerifier (source, parent) << endl;
}

pvector<long int> DynGraph::bfs_scratch (long int source) {
	std::queue<long int> q;
	pvector<long int> parent(this->N);
	for (long int i = 0; i < this->N; i++) {
		parent[i] = -1;
	}	
	parent[source] = source; this->depths_bfs_tree[source] = 0;
	q.push(source);
	while (!q.empty()) {
		for (auto &j : this->Neighbors[q.front()]) {
			if (parent[j] < 0) {
				parent[j] = q.front(); this->depths_bfs_tree[j] = this->depths_bfs_tree[q.front()] + 1;
				q.push(j);
			}
		}
		q.pop();
	}
	cout << "BFSVerifier output: " 
	 	 << boolalpha << BFSVerifier (source, parent) << endl;
	return parent;
}


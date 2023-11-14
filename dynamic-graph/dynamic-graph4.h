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
	vector<vector<long int> > bfs_vector;
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
		}
	}
	void bfsFromScratch (long int source);
	void bfs (long int source);
	void printNeighbors () {
		for (auto &i : this->Neighbors) {
			for (auto &j : i)
				cout << j << " ";
			cout << endl;
		}
	}
	void buildGraph (vector<edge>& edge_list, vector<edge>& parents_list, long int start_line, long int num_lines);	
	void bfsDelta (vector<edge>& edge_list, vector<edge>& parents_list, long int start_line, long int num_lines);	
	void formBFS ();
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

	if (start_line != 0) {
		timer_start = get_wall_time();
		bfsDelta (edge_list, parents_list, start_line, num_lines);
		timer_end = get_wall_time();
		cout << "Wall time to preprocess: " 
			<< (timer_end - timer_start) << endl;
	}
	formBFS();
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


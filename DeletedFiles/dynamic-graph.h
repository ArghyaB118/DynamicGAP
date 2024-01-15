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
bool verbose = true;

struct bfsNode {
	// long int value;
	vector<long int> children;
	long int parent;
	long int depth;
};

class DynGraph {
private:
	long int N; // number of nodes/vertices
	long int E; // number of edges
	vector<neighbors> Parents, Neighbors;
	vector<vector<long int> > p, e;
	vector<vector<long int>> bfs_tree;
	vector<bool> visited;
	// std::unordered_map<long int, bfsNode> bfs_tree;
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
			bfs_tree.push_back(tmp);
			visited.push_back(false);
		}
	}
	void bfsRecursive (long int source, long int depth);
	void printNeighbors () {
		for (auto &i : this->Neighbors) {
			for (auto &j : i)
				cout << j << " ";
			cout << endl;
		}
	}
	void buildGraph (vector<edge>& edge_list, vector<edge>& parents_list, long int start_line, long int num_lines);	
};

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
	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
		std::sort (p[i].begin(), p[i].end());
		std::sort (e[i].begin(), e[i].end());
	}
	double timer_start = 0, timer_end = 0;
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
}

void DynGraph::bfsRecursive (long int source, long int depth) {
	this->bfs_tree[depth].push_back(source);
	this->visited[source] = true;
	for (auto& i : this->Neighbors[source]) {
		if (!visited[i])
			bfsRecursive (i, depth + 1);
	}
}

/*void DynGraph::bfsFromScratch (long int source) {
	// printNeighbors();
	queue<long int> q;
	vector<bool> visited(this->N, false);
	q.push(source); visited[source] = true;
	
	bfsNode tmp_bfs_node; 
	tmp_bfs_node.parent = -1; 
	tmp_bfs_node.depth = 0; 
	tmp_bfs_node.children = {};
	this->bfs_tree[source] = tmp_bfs_node;
	tmp_bfs_node.children.clear();
	
	if (verbose) { cout << "BFS: "; }
	while (!q.empty()) {
		if (verbose) { cout << q.front() << " "; }
		for (auto &i : Neighbors[q.front()]) {
			if (!visited[i]) {
				q.push(i);
				visited[i] = true;
				tmp_bfs_node.parent = q.front();
				tmp_bfs_node.depth = bfs_tree[q.front()].depth + 1;
				tmp_bfs_node.children.push_back(i);
				this->bfs_tree[i] = tmp_bfs_node;
				tmp_bfs_node.children.clear();
			}
		}
		q.pop();
	}
	if (verbose) { cout << endl; }
}*/


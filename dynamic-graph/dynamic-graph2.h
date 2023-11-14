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

struct bfs_node {
	long int val;
	struct bfs_node *parent;
	vector<bfs_node *> children;
};

/*
Container:
1. quick lookup in the bfs tree, and find pointers to any particular node.
for that, we keep a hashmap.
2. pointer to parent, to check depth quickly.
3. for bfs, we keep pointers to the children.

Algorithm: incoming edge has -
1. parent with less depth, reconnect.
2. new child, make a new node, start bfs there.
	a. during bfs, form new nodes as necessary.
	b. ignore or reconnect existing nodes based on the depth.
3. new parent, ignore.
*/

class DynGraph {
private:
	long int N; // number of nodes/vertices
	long int E; // number of edges
	vector<neighbors> Parents, Neighbors;
	vector<vector<long int> > p, e;
	unordered_map<long int, struct bfs_node*> node_pointers;
	struct bfs_node *root = new bfs_node;

	bfs_node *new_bfs_node (long int val, bfs_node *parent);
	bfs_node *get_bfs_node (long int val);
	long int get_depth (bfs_node *node);
	void set_bfs_node (long int val, bfs_node *node, bfs_node *parent);

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
		}
	}
	void bfsFromScratch (long int source);
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

	timer_start = get_wall_time();
	if (start_line != 0) {
		for (long int i = start_line; i < num_lines; i++) {
			if (this->node_pointers.find(edge_list[i].first) 
				== this->node_pointers.end())
				continue;
			else if (this->node_pointers.find(edge_list[i].first) 
				!= this->node_pointers.end() &&
				this->node_pointers.find(edge_list[i].second) 
				!= this->node_pointers.end()) {
				bfs_node *node = get_bfs_node(edge_list[i].second);
				bfs_node *current_parent = node->parent;
				bfs_node *potential_parent = get_bfs_node(edge_list[i].first);

				if (get_depth(potential_parent) + 1 < get_depth(node)) {
					auto it = find(current_parent->children.begin(), 
						current_parent->children.end(),
		                node);
					current_parent->children.erase(it);
					node->parent = potential_parent;
					potential_parent->children.push_back(node);
				}
			}
			else if (this->node_pointers.find(edge_list[i].first) 
				!= this->node_pointers.end() &&
				this->node_pointers.find(edge_list[i].second) 
				== this->node_pointers.end()) {
				bfs_node *node = new_bfs_node(edge_list[i].second, 
					get_bfs_node(edge_list[i].first));
				queue<long int> q; q.push(edge_list[i].second);
				while (!q.empty()) {
					// #pragma omp parallel for
					for (auto &j : Neighbors[q.front()]) {
						bfs_node *potential_parent = get_bfs_node(q.front());

						if (this->node_pointers.find(j) 
							== this->node_pointers.end()) {
							new_bfs_node(j, potential_parent);
							q.push(j);
						}
						else if (get_depth(potential_parent) + 1 
								< get_depth(get_bfs_node(j))) {
							bfs_node *node = get_bfs_node(j);
							auto it = find(node->parent->children.begin(), 
								node->parent->children.end(),
		                		node);
							node->parent->children.erase(it);
							node->parent = potential_parent;
							potential_parent->children.push_back(node);
						}
					}
					q.pop();
				}
			}
		}
	}
	timer_end = get_wall_time();
	cout << "Wall time to preprocess: " 
		<< (timer_end - timer_start) << endl;
}

bfs_node* DynGraph::new_bfs_node (long int val, bfs_node *parent) {
	bfs_node *tmp = new bfs_node;
	tmp->val = val;
	tmp->parent = parent;
	tmp->children = {};
	this->node_pointers[val] = tmp;
	if (parent != NULL)
		parent->children.push_back(tmp);
	return tmp;
}

void DynGraph::set_bfs_node (long int val, bfs_node *node, bfs_node *parent) {
	node->val = val;
	node->parent = parent;
	node->children = {};
	node_pointers[val] = node;
	if (parent != NULL)
		parent->children.push_back(node);
}

bfs_node* DynGraph::get_bfs_node (long int val) {
	return node_pointers[val];
}

long int DynGraph::get_depth (bfs_node *node) {
	long int depth = -1;
	while (node != NULL) {
		node = node->parent;
		depth++;
	}
	return depth;
}

void DynGraph::bfsFromScratch (long int source) {
	queue<long int> q;
	vector<bool> visited(this->N, false);
	q.push(source); visited[source] = true;
	this->root = new_bfs_node(source, NULL);

	if (verbose) { cout << "BFS: "; }
	while (!q.empty()) {
		if (verbose) { cout << q.front() << " "; }
		for (auto &i : Neighbors[q.front()]) {
			if (!visited[i]) {
				q.push(i);
				visited[i] = true;
				new_bfs_node (i, get_bfs_node(q.front()));
			}
		}
		q.pop();
	}
	if (verbose) { cout << endl; }
}


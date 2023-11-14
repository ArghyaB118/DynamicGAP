#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;
unordered_map<long int, struct bfs_node*> node_pointers;

struct bfs_node {
	long int val;
	struct bfs_node *parent;
	vector<bfs_node *> children;
};

bfs_node *new_bfs_node (long int val, bfs_node *parent) {
	bfs_node *tmp = new bfs_node;
	tmp->val = val;
	tmp->parent = parent;
	tmp->children = {};
	node_pointers[val] = tmp;
	return tmp;
}

void set_bfs_node (long int val, bfs_node *node, bfs_node *parent) {
	node->val = val;
	node->parent = parent;
	node->children = {};
	node_pointers[val] = node;
}

bfs_node *get_bfs_node (long int val) {
	return node_pointers[val];
}

long int get_depth (bfs_node *node) {
	long int depth = -1;
	while (node != NULL) {
		node = node->parent;
		depth++;
	}
	return depth;
}

int main () {
	// struct bfs_node *root = new_bfs_node(0, NULL);
	struct bfs_node *root = new bfs_node;
	set_bfs_node (0, root, NULL);
	cout << root->val << endl;
	root->children.push_back(new_bfs_node(1, root));
	root->children.push_back(new_bfs_node(2, get_bfs_node(0)));
	bfs_node *node = get_bfs_node(2);
	node->children.push_back(new_bfs_node(3, node));
	node->children.push_back(new_bfs_node(4, node));
	for (auto i : root->children)
		cout << i->val << " ";
	cout << endl;

	for (auto i : node->children)
		cout << i->val << " ";
	cout << endl;

	cout << get_depth(get_bfs_node(0)) << endl;
	cout << get_depth(get_bfs_node(2)) << endl;
	cout << get_depth(get_bfs_node(4)) << endl;
	
	return 0;
}
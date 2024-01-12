#include <iostream>
#include <vector>
#include <string>
#include <fstream>
using namespace std;

void buildGraph(int m, int n, int &e) {
	ofstream graph_file;
	graph_file.open ("synthetic_graph.el");
	int i = 1;
	while (2*i + 1 <= m) {
		graph_file << i << " " << 2*i << endl;
		graph_file << i << " " << 2*i + 1 << endl;
		i++;
		e += 2;
	}
	graph_file.close();
}

void connectLeaves (int m, int n, int& e) {
	ofstream graph_file;
	graph_file.open ("synthetic_graph.el", std::ios::app);
	int i = 1;
	while (2*i + 1 <= m) {
		i++;
	}
	while (i < m) {
		graph_file << i << " " << i + 1 << endl;
		i++;
		e++;
	}
	graph_file.close();
}

void changeAroundHead (int m, int n, int& e) {
	ofstream graph_file;
	graph_file.open ("synthetic_graph.el", std::ios::app);
	graph_file << 1 << " " << 4 << endl;
	graph_file << 1 << " " << 5 << endl;
	graph_file << 1 << " " << 6 << endl;
	graph_file << 1 << " " << 7 << endl;
	e += 4;
	graph_file.close();
}

int main () {
	// number of nodes
	int m = 3072627; // size of orkut
	// n-ary tree
	int n = 2;
	// number of edges
	int e = 0;
	
	buildGraph (m, n, e);
	cout << "Total edge created in the beginning: " << e << endl;
	connectLeaves (m, n, e);
	cout << "Total edge created after joining leaves: " << e << endl;
	
	changeAroundHead (m, n, e);
	cout << "Total edge created after maing the change around head: " << e << endl;
	return 0;
}
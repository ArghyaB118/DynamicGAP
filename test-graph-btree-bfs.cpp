#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <unordered_set>
#include <utility>
#include <time.h>
#include <sys/time.h>
#include "cpp-btree-1.0.1/btree_set.h"
using namespace std;

typedef pair<long int, long int> edge;
typedef btree::btree_set<edge> edgeList;


// Bad idea:  will throw this out later.
static long int bfs_source = 0;

// Note: can't make it member function because
// reference to non-static member function must be called.
// inline defined, not in use anymore.
bool matchesSource (const edge& edge) {
	return (edge.first == bfs_source);
}

class Graph {
private:
	edgeList edge_list;
public:
	void addEdge(edge edge) {
		this->edge_list.insert(edge);
	}
	long int num_edges() {
		return this->edge_list.size();
	}
	void printEdges() {
		for (auto &e : this->edge_list)
			cout << e.first << " " << e.second << endl;
	}
	void bfs (long int source);
};

void Graph::bfs (long int source) {
	queue<long int> q;
	unordered_set<long int> visited;
	q.push(source);
	while (!q.empty()) {
		bfs_source = q.front();
		q.pop(); cout << bfs_source << endl;
		auto it = find_if(this->edge_list.begin(), 
			this->edge_list.end(), 
			[](edge& e) { return e.first == bfs_source; });
		// auto it = find_if(this->edge_list.begin(), 
		//	this->edge_list.end(), matchesSource);
		while (it != this->edge_list.end()) {
			if (it->first != bfs_source)
				break;
			if (visited.find(it->second) == visited.end()) {
				q.push(it->second);
				visited.insert(it->second);
			}
			it++;
		}
	}
}

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

void readEdges (ifstream& myfile, 
	Graph& g, // Graph& g, // edgeList& edge_list, 
	long int num_lines) {
	string myline;
	if (myfile.is_open()) {	
		while (myfile && num_lines > 0) {
			getline(myfile, myline);
			vector<string> tmp; 
			tmp.push_back("");
			for (char &c : myline) {
				if (c != ' ')
					tmp.back().push_back(c);
				else
					tmp.push_back("");
			}
			if (tmp[0] != "" && tmp[1] != "") {
				edge tmp_edge = {stol(tmp[0]), stol(tmp[1])};
				g.addEdge(tmp_edge);
				// edge_list.insert(tmp_edge);
			}	
			myline.clear();
			num_lines--;
		}
	}
	return;
}

void run_benchmark (string filename, 
	long int num_edges_init_round, 
	long int num_edges_each_round,
	long int num_rounds) {
	// read the file name from command line input.
	ifstream myfile (filename);
	Graph g;
	vector<double> timestamps = {};
	timestamps.push_back(get_wall_time());
	readEdges(myfile, g, num_edges_init_round);
	timestamps.push_back(get_wall_time());
	cout << "Wall time to read edges in the initial round: "
	 << timestamps.back() - timestamps.end()[-2] << endl;
	// cout << edge_list.size() << endl;
	cout << "Number of edges after initial round: " 
	 << g.num_edges() << endl;
	g.bfs(0); // I know that the first entry is [0, 13]
	timestamps.push_back(get_wall_time());
	cout << "Wall time for bfs after the initial round: "
	 << timestamps.back() - timestamps.end()[-2] << endl;
	for (long int i = 0; i < num_rounds; i++) {
		timestamps.push_back(get_wall_time());
		readEdges(myfile, g, num_edges_each_round);
		timestamps.push_back(get_wall_time());
		cout << "Wall time to read edges in the " 
		 << i + 1 << "-th round: "
	 	 << timestamps.back() - timestamps.end()[-2] << endl;
		// cout << edge_list.size() << endl;
		cout << "Number of edges after " << i + 1
		 << "-th round: " << g.num_edges() << endl;
		g.bfs(0);
		timestamps.push_back(get_wall_time());
		cout << "Wall time for bfs after " << i + 1
		 << "-th round: "
	 	 << timestamps.back() - timestamps.end()[-2] << endl;
	}
	myfile.close();
	g.printEdges();
}

bool sanity_check (string filename, 
	long int num_edges_init_round, 
	long int num_edges_each_round,
	long int num_rounds) {
	long int num_lines = 0;
	ifstream myfile (filename);
	if (myfile.is_open()) {	
		while (myfile) {
			string tmp;
			getline(myfile, tmp);
			num_lines++;
		}
	}
	if (num_lines >= num_edges_init_round + 
		num_edges_each_round * num_rounds)
		return true;
	return false;
}

int main (int argc, char *argv[]) {
	if (argc < 5) {
		cout << "Correct usage: ./a.out "
		 << "<filename (file in .el format)> "
		 << "<num_edges_init_round> " 
		 <<	"<num_edges_each_round> " 
		 <<	"<num_rounds>" << endl;
		 return 0;
	}
	string filename = argv[1];
	long int num_edges_init_round = stol(argv[2]);
	long int num_edges_each_round = stol(argv[3]);
	long int num_rounds = stol(argv[4]);
	if (!sanity_check(filename, 
		num_edges_init_round, 
		num_edges_each_round,
		num_rounds))
		return 0;
	run_benchmark(filename, 
		num_edges_init_round, 
		num_edges_each_round,
		num_rounds);
	return 0;
}
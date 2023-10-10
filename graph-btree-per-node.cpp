#include <iostream>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <utility>
#include <cstdlib>

#include "graph.h"
#include "sanity-check.h"

using namespace std;

bool verify = false;

vector<edge> readEdges (ifstream& myfile, 
	long int num_lines) {
	vector<edge> edge_list = {};
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
				edge_list.push_back(tmp_edge);
			}
			myline.clear();
			num_lines--;
		}
	}
	return edge_list;
}

void updateEdges (vector<edge>& edge_list, 
	Graph& g, // edgeList& edge_list, 
	long int start_line,
	long int num_lines) {
	double timer_start = 0, timer_end = 0;
	double update_time = 0;
	#pragma omp for
	for (auto i = 0; i < num_lines; i++) {
		timer_start = get_wall_time();
		g.addEdge(edge_list[start_line + i]);
		timer_end = get_wall_time();
		update_time += (timer_end - timer_start);
	}
	cout << "Wall time to update edges: " << update_time << endl;
	return;
}

void run_benchmark (long int source,
	string filename, 
	long int num_edges_init_round, 
	long int num_edges_each_round,
	long int num_rounds) {
	// read the file name from command line input.
	ifstream myfile (filename);
	long int num_nodes = countNodes(filename, 
		num_edges_init_round, 
		num_edges_each_round, 
		num_rounds);
	cout << "Number of nodes: " 
	 << num_nodes << endl;
	vector<edge> edge_list = 
	 readEdges(myfile, num_edges_init_round 
	 	+ num_edges_each_round * num_rounds);
	Graph g(num_nodes);
	vector<double> timestamps = {};
	timestamps.push_back(get_wall_time());
	updateEdges(edge_list, g, 0, num_edges_init_round);
	timestamps.push_back(get_wall_time());
	cout << "Wall time to read edges in the initial round: "
	 << timestamps.back() - timestamps.end()[-2] << endl;
	if (verify) {
		g.bfs(source);
		timestamps.push_back(get_wall_time());
		cout << "Wall time for bfs after the initial round: "
		 << timestamps.back() - timestamps.end()[-2] << endl;
	}
	pvector<long int> parent = g.bfs_gap(source); 
	timestamps.push_back(get_wall_time());
	cout << "Wall time for bfs_gap after the initial round: "
	 << timestamps.back() - timestamps.end()[-2] << endl;
	if (verify) {
		cout << "BFSVerifier output: " 
		 << boolalpha << g.BFSVerifier(source, parent) << endl;
	}
	for (long int i = 0; i < num_rounds; i++) {
		timestamps.push_back(get_wall_time());
		updateEdges(edge_list, g, 
			num_edges_init_round + i * num_edges_each_round, 
			num_edges_each_round);
		timestamps.push_back(get_wall_time());
		cout << "Wall time to read edges in the " 
		 << i + 1 << "-th round: "
	 	 << timestamps.back() - timestamps.end()[-2] << endl;
		if (verify) {
			g.bfs(source);
			timestamps.push_back(get_wall_time());
			cout << "Wall time for bfs after " << i + 1
			 << "-th round: "
		 	 << timestamps.back() - timestamps.end()[-2] << endl;
		}
	 	parent = g.bfs_gap(source);
	 	timestamps.push_back(get_wall_time());
		cout << "Wall time for bfs_gap after " << i + 1
		 << "-th round: "
	 	 << timestamps.back() - timestamps.end()[-2] << endl;
	 	if (verify) {
		 	cout << "BFSVerifier output: " 
		 	 << boolalpha << g.BFSVerifier(source, parent) << endl;
	 	}
	}
	myfile.close();
}

int main (int argc, char *argv[]) {
	if (argc < 7) {
		cout << "Correct usage: ./a.out "
		 << "<-random> OR <-deterministic>"
		 << "<number of points> OR <source>"
		 << "<filename (file in .el format)> "
		 << "<num_edges_init_round> " 
		 <<	"<num_edges_each_round> " 
		 <<	"<num_rounds>" << endl;
		 return 0;
	}
	string experiment_type = argv[1];
	long int experiment_setup = stol(argv[2]);
	string filename = argv[3];
	long int num_edges_init_round = stol(argv[4]);
	long int num_edges_each_round = stol(argv[5]);
	long int num_rounds = stol(argv[6]);
	if (!sanity_check(filename, 
		num_edges_init_round, 
		num_edges_each_round,
		num_rounds))
		return 0;
	if (experiment_type == "-deterministic") {
		run_benchmark(experiment_setup,
			filename, 
			num_edges_init_round, 
			num_edges_each_round,
			num_rounds);
	}
	else if (experiment_type == "-random") {
		for (long int i = 0; i < experiment_setup; i++) {
			long int random_source = 
				rand() % countNodes(filename, num_edges_init_round, 
							num_edges_each_round, num_rounds);
			run_benchmark(random_source,
				filename, 
				num_edges_init_round, 
				num_edges_each_round,
				num_rounds);
		}
	} 
	return 0;
}
#include <iostream>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <utility>
#include <cstdlib>

#include "dynamic-graph/dynamic-graph4.h"
#include "sanity-check.h"

using namespace std;

bool verify = true;

void readEdgesBulkLoad (ifstream& myfile,
	vector<edge>& edge_list,
	vector<edge>& parents_list,
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
				edge tmp_parent = {stol(tmp[1]), stol(tmp[0])};
				edge_list.push_back(tmp_edge);
				parents_list.push_back(tmp_parent);
			}
			myline.clear();
			num_lines--;
		}
	}
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
	vector<edge> edge_list = {};
	edge_list.reserve (num_edges_init_round 
		+ num_edges_each_round * num_rounds);
	vector<edge> parents_list = {};
	parents_list.reserve (num_edges_init_round 
		+ num_edges_each_round * num_rounds);
	readEdgesBulkLoad (myfile, edge_list, parents_list,
		num_edges_init_round 
	 	+ num_edges_each_round * num_rounds);
	
	DynGraph g(num_nodes);
	vector<double> timestamps = {};
	timestamps.push_back(get_wall_time());
	g.buildGraph (edge_list, parents_list, 0, num_edges_init_round);
	timestamps.push_back(get_wall_time());
	cout << "Wall time to read edges in the initial round: "
	 << timestamps.back() - timestamps.end()[-2] << endl;
	if (verify) {
		g.bfsFromScratch (source);
		timestamps.push_back(get_wall_time());
		cout << "Wall time for bfs after the initial round: "
		 << timestamps.back() - timestamps.end()[-2] << endl;
	}
	for (long int i = 0; i < num_rounds; i++) {
		timestamps.push_back(get_wall_time());
		g.buildGraph (edge_list, parents_list,
			num_edges_init_round + i * num_edges_each_round,
			num_edges_init_round + (i + 1) * num_edges_each_round);
		timestamps.push_back(get_wall_time());
		cout << "Wall time to read edges in the " 
		 << i + 1 << "-th round: "
	 	 << timestamps.back() - timestamps.end()[-2] << endl;
		if (verify) {
			g.bfsFromScratch (source);
			timestamps.push_back(get_wall_time());
			cout << "Wall time for bfs after " << i + 1
			 << "-th round: "
		 	 << timestamps.back() - timestamps.end()[-2] << endl;
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

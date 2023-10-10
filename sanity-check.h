#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>

using namespace std;

typedef pair<long int, long int> p;

p countNodesAndEdgesDynamically (string filename,
	long int num_edges_init_round, 
	long int num_edges_each_round,
	long int num_rounds) {
	unordered_map<long int, vector<long int>> neighbors;
	ifstream myfile (filename);
	string myline;
	long int num_lines = num_edges_init_round + 
	 num_edges_each_round * num_rounds;
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
				if (neighbors[stol(tmp[0])].size() > 0) {
					neighbors[stol(tmp[0])].push_back(stol(tmp[1]));
				}
				else {
					neighbors[stol(tmp[0])] = {stol(tmp[1])};
				}
			}
			myline.clear();
			num_lines--;
		}
	}
	myfile.close();
	long int node_count = 0, edge_count = 0;
	for (auto & i : neighbors) {
		node_count++;
		for (auto& j : i.second) {
			edge_count++;
		}
	}
	return {node_count, edge_count};
}

// This is used as there may be some numbered nodes missing 
// between 0 and the max_num.
long int countNodes (string filename,
	long int num_edges_init_round, 
	long int num_edges_each_round,
	long int num_rounds) {
	ifstream myfile (filename);
	string myline;
	long int max_num = 0;
	long int num_lines = num_edges_init_round + 
	 num_edges_each_round * num_rounds;
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
			if (tmp[0] != "" && stol(tmp[0]) > max_num)
				max_num = stol(tmp[0]);
			if (tmp[1] != "" && stol(tmp[1]) > max_num)
				max_num = stol(tmp[1]);
			myline.clear();
			num_lines--;
		}
	}
	myfile.close();
	return max_num + 1;
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
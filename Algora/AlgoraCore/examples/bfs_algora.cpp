/**
 * Copyright (C) 2013 - 2020 : Kathrin Hanauer
 *
 * This file is part of Algora.
 *
 * Algora is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Algora is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Algora.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact information:
 *   http://algora.xaikal.org
 */

#include "graph.incidencelist/incidencelistgraph.h"
#include "algorithm.basic.traversal/breadthfirstsearch.h"
#include "property/fastpropertymap.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include </home/arghya/DynamicGAP/util.h>

using namespace std;
using namespace Algora;

long int countNodes (string filename, long int num_lines) {
	ifstream myfile (filename);
	string myline;
	long int max_num = 0;
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

std::map<int, Vertex*> makeVertices (string filename, 
	IncidenceListGraph& g, long int num_lines) {
	std::map<int, Vertex*> vertices;
	long int count_nodes = countNodes(filename, num_lines);
	for (auto i = 0; i < count_nodes; i++) {
		auto u = g.addVertex();
		u->setName(std::to_string(i));
		vertices[i] = u;
	}
	return vertices;
}

void readEdges (string filename, 
	IncidenceListGraph& g, 
	std::map<int, Vertex*> &vertices, 
	long int num_lines) {
	ifstream myfile (filename);
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
			if (tmp[0] != "" && tmp[1] != "")
				g.addArc(vertices[stoi(tmp[0])], vertices[stoi(tmp[1])]);
			myline.clear();
			num_lines--;
		}
	}
}

int main(int argc, char *argv[]) {
	string filename = argv[1];
	long int num_edges_init_round = stol(argv[2]);
	long int source = stol(argv[3]);
	
	long int num_lines = num_edges_init_round;
	// First, we create a graph.
	// Calling setName() is optional, we only use it here to get a nicer output.
	vector<double> timestamps = {};
	timestamps.push_back(get_wall_time());
	IncidenceListGraph g;
	std::map<int, Vertex*> vertices = makeVertices(filename, g, num_lines);
	readEdges(filename, g, vertices, num_lines);
	
	// Configure the breadth-first search algorithm:
	BreadthFirstSearch<FastPropertyMap> bfs;
	bfs.setGraph(&g);
	bfs.setStartVertex(vertices[source]);
	// We want it to store the BFS number of each vertex in this map:
	FastPropertyMap<DiGraph::size_type> bfsNum(0);
	bfs.useModifiableProperty(&bfsNum);

	// We want to be notified whenever the algorithm discovers a new tree arc:
	bfs.onTreeArcDiscover([](Arc *a) { std::cout << "Discovered tree arc " << a << std::endl; });

	// prepare() checks the algorithm's preconditions:
	if (!bfs.prepare()) {
		std::cerr << "Could not prepare BFS algorithm. Please check the preconditions." << std::endl;
		return 1;
	}
	timestamps.push_back(get_wall_time());
	cout << "Wall time to read edges in the initial round: "
	 << timestamps.back() - timestamps.end()[-2] << endl;
	
	// run() runs the algorithm, whereas deliver() only returns the result computed during run()
	std::cout << "Running breadth-first search..." << std::endl;
	bfs.run();
	timestamps.push_back(get_wall_time());
	cout << "Wall time for bfs after the initial round: "
	 << timestamps.back() - timestamps.end()[-2] << endl;
	auto numReachedVertices = bfs.deliver();
	std::cout << "Done, reached " << numReachedVertices << " vertices." << std::endl;

	// Finally, output the BFS number for each vertex:
	g.mapVertices([&bfsNum](Vertex *v) { std::cout << "BFS number of " << v << ": " << bfsNum(v) << std::endl; });

	return 0;
}

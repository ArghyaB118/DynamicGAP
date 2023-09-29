#include <iostream>
#include <vector>

#include "util.h"
using namespace std;

int main () {
	vector<double> timestamps = {};
	system ("g++ -std=c++11 -O3 -Wall -fopenmp src-gap/bfs.cc");
	timestamps.push_back(get_wall_time());
	system ("./a.out -f graphs/slashdot.el -n 1");
	timestamps.push_back(get_wall_time());
	cout << "Wall time for reading the graph and running bfs in GAP: "
	 << timestamps.back() - timestamps.end()[-2] << endl;
	return 0;
}
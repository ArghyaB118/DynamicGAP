#include <iostream>
#include <vector>

#include "util.h"
using namespace std;

int main () {
	system ("vmtouch -vt graphs/slashdot.el");
	vector<double> timestamps = {};
	timestamps.push_back(get_wall_time());
	system ("./src-gap/bfs -f graphs/slashdot.el -n 1");
	timestamps.push_back(get_wall_time());
	cout << "Wall time for reading the graph and running bfs in GAP: "
	 << timestamps.back() - timestamps.end()[-2] << endl;
	return 0;
}
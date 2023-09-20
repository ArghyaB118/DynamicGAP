#include <iostream>
using namespace std;

class Graph {
private:
	int n = 0;
public:
	int getN () { return this->n; }
	void setN (int n) {
		this->n = n;
	}
};

int main () {
	Graph g;
	cout << g.getN() << endl;
	g.setN(5);
	cout << g.getN() << endl;
	return 0;
}


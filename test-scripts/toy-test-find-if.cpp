#include <iostream>
#include <vector>
using namespace std;

int main() {
	vector<int> arr = {12, 3, 5, 10, 6, 7, 2, 4};
	auto it = find_if(arr.begin(), arr.end(), 
		[](int& i) { return (i % 2 == 1); });
	while (it != arr.end()) {
		if (*it == 6)
			break;
		cout << *it << endl;
		it++;
	}
	return 0;
}
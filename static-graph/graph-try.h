#include <iostream>
#include <string>
#include <vector>
#include <queue>

#include "../src-gap/pvector.h"
#include "../src-gap/sliding_queue.h"
#include "../src-gap/bitmap.h"

#include "../util.h"
#include "../cpp-btree-1.0.1/btree_set.h"

using namespace std;

typedef pair<long int, long int> edge;
typedef btree::btree_set<long int> neighbors;

class Graph {
private:
	long int N; // number of nodes/vertices
	long int E; // number of edges
	vector<neighbors> Parents;
	vector<neighbors> Neighbors;
	vector<float> property;
	vector<bool> affected;
public:
	Graph(int N) {
		this->N = N; // number of nodes/vertices
		this->E = 0;
		for (long int i = 0; i < N; i++) {
			Neighbors.push_back({});
			Parents.push_back({});
		}
	}
	void addEdge (edge edge) {
		if (this->Neighbors[edge.first].find(edge.second)
			== this->Neighbors[edge.first].end()) {
			this->Neighbors[edge.first].insert(edge.second);
			this->Parents[edge.second].insert(edge.first);
			this->E++;
		}
	}
	void bfs (long int source);
	void printNeighbors () {
		for (auto &i : this->Neighbors) {
			for (auto &j : i)
				cout << j << " ";
			cout << endl;
		}
	}
	pvector<long int> bfs_gap (long int source, int alpha, int beta);
	void QueueToBitmap(const SlidingQueue<long int> &q, Bitmap &bm);
	void BitmapToQueue(const Bitmap &bm, SlidingQueue<long int> &q); 
	long int BUStep(pvector<long int> &parent, Bitmap &front, Bitmap &next);
	long int TDStep(pvector<long int> &parent, SlidingQueue<long int> &q);
	bool BFSVerifier(long int source, const pvector<long int> &parent);
	void BFSStartFromScratch(long int source);
	void dynBFSAlg(long int source);
	void BFSIter0(SlidingQueue<long int>& queue);
};

void Graph::bfs (long int source) {
	// printNeighbors();
	queue<long int> q;
	vector<bool> visited(this->N, false);
	q.push(source); visited[source] = true;
	cout << "BFS: ";
	while (!q.empty()) {
		cout << q.front() << " ";
		for (auto &i : Neighbors[q.front()]) {
			if (!visited[i]) {
				q.push(i);
				visited[i] = true;
			}
		}
		q.pop();
	}
	cout << endl;
}

void Graph::QueueToBitmap(const SlidingQueue<long int> &q, Bitmap &bm) {
	#pragma omp parallel for
	for (auto q_iter = q.begin(); q_iter < q.end(); q_iter++) {
		long int u = *q_iter;
		bm.set_bit_atomic(u);
	}
}

void Graph::BitmapToQueue(const Bitmap &bm,
                   SlidingQueue<long int> &q) {
  #pragma omp parallel
  {
    QueueBuffer<long int> lqueue(q);
    #pragma omp for nowait
    for (long int n=0; n < this->N; n++)
      if (bm.get_bit(n))
        lqueue.push_back(n);
    lqueue.flush();
  }
  q.slide_window();
}

long int Graph::BUStep(pvector<long int> &parent, Bitmap &front, Bitmap &next) {
  long int awake_count = 0;
  next.reset();
  #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
  for (long int u=0; u < this->N; u++) {
    if (parent[u] < 0) {
      for (long int v : this->Parents[u]) {
        if (front.get_bit(v)) {
          parent[u] = v;
          awake_count++;
          next.set_bit(u);
          break;
        }
      }
    }
  }
  return awake_count;
}

long int Graph::TDStep(pvector<long int> &parent, SlidingQueue<long int> &q) {
  long int scout_count = 0;
  #pragma omp parallel
  {
    QueueBuffer<long int> lqueue(q);
    #pragma omp for reduction(+ : scout_count) nowait
    for (auto q_iter = q.begin(); q_iter < q.end(); q_iter++) {
      long int u = *q_iter;
      for (long int v : this->Neighbors[u]) {
        long int curr_val = parent[v];
        if (curr_val < 0) {
          if (compare_and_swap(parent[v], curr_val, u)) {
            lqueue.push_back(v);
            scout_count += -curr_val;
          }
        }
      }
    }
    lqueue.flush();
  }
  return scout_count;
}

pvector<long int> Graph::bfs_gap (long int source, 
	int alpha = 15, 
	int beta = 18) {
	pvector<long int> parent(this->N);
	#pragma omp parallel for
	for (long int i = 0; i < this->N; i++) {
	    parent[i] = this->Neighbors[i].size() != 0 
			? -this->Neighbors[i].size() 
			: -1;
	}
	vector<double> timestamps1 = {};
	parent[source] = source;
	SlidingQueue<long int> q(this->N);
	q.push_back(source);
	q.slide_window();
	Bitmap curr(this->N);
	curr.reset();
  Bitmap front(this->N);
	front.reset();
	long int edges_to_check = this->E;
	long int scout_count = this->Neighbors[source].size();
	timestamps1.push_back(get_wall_time());
	while (!q.empty()) {
		if (scout_count > edges_to_check / alpha) {
			long int awake_count, old_awake_count;
			QueueToBitmap(q, front);
			awake_count = q.size();
			q.slide_window();
			do {
				old_awake_count = awake_count;
				awake_count = BUStep(parent, front, curr);
				front.swap(curr);
			} while ((awake_count >= old_awake_count) 
				|| (awake_count > this->N / beta));
			BitmapToQueue(front, q);
			scout_count = 1;
		} else {
			edges_to_check -= scout_count;
			scout_count = TDStep(parent, q);
			q.slide_window();
		}
	}
	timestamps1.push_back(get_wall_time());
	cout << "Real wall time for bfs_gap: "
	 << timestamps1.back() - timestamps1.end()[-2] << endl;
	#pragma omp parallel for
	for (long int n = 0; n < this->N; n++)
		if (parent[n] < -1) 
	  		parent[n] = -1;
	return parent;
}

void Graph::BFSIter0(SlidingQueue<long int>& queue){  
    pvector<bool> visited(this->N, false);     
  
    #pragma omp parallel     
    {
        QueueBuffer<long int> lqueue(queue);
        #pragma omp for schedule(dynamic, 64)
        for(long int n=0; n < this->N; n++){
            if(this->affected[n]){
                float old_depth = this->property[n];
                float new_depth = std::numeric_limits<float>::max();

                // pull new depth from incoming neighbors
                for(auto v: this->Parents[n]){
                    if (this->property[v] != -1) {
                        new_depth = std::min(new_depth, this->property[v] + 1);
                    }
                }
                
                // trigger happens if it is:
                // 1) brand new vertex with old_prop = -1 and we found a new valid min depth 
                // 2) already existing vertex and we found a new depth smaller than old depth 
                bool trigger = (
                ((new_depth < old_depth) || (old_depth == -1)) 
                && (new_depth != std::numeric_limits<float>::max())                 
                );               

                /*if(trigger){                                                 
                    ds->property[n] = new_depth; 
                    for(auto v: out_neigh(n, dataStruc, ds, directed)){
                        float curr_depth = ds->property[v];
                        float updated_depth = ds->property[n] + 1;                        
                        if((updated_depth < curr_depth) || (curr_depth == -1)){   
                            if(compare_and_swap(ds->property[v], curr_depth, updated_depth)){                                                              
                                lqueue.push_back(v); 
                            }
                        }
                    }
                }*/

                // Note: above is commented and included this new thing. 
                // Above was leading to vertices being queued redundantly
                // Above assumes updated_depth < curr_depth only once. 
                // This is not true in dynamic case because we start from affected vertices
                // whose depths are not all necessary the same.
                // In static version, the above works because static version starts from the source 
                // and we know that updated_depth < curr_depth only once. 

                if(trigger){
                    this->property[n] = new_depth; 
                    for(auto v: this->Neighbors[n]){
                        float curr_depth = this->property[v];
                        float updated_depth = this->property[n] + 1;
                        if((updated_depth < curr_depth) || (curr_depth == -1)){
                            bool curr_val = visited[v];
                            if(!curr_val){
                                if(compare_and_swap(visited[v], curr_val, true))
                                    lqueue.push_back(v);
                            }
                            while(!compare_and_swap(this->property[v], curr_depth, updated_depth)){
                                curr_depth = this->property[v];
                                if(curr_depth <= updated_depth){
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        lqueue.flush();
    }   
}


void Graph::dynBFSAlg(long int source){
    std::cout <<"Running dynamic BFS " << std::endl;
    SlidingQueue<long int> queue(this->N);
    if(this->property[source] == -1) this->property[source] = 0;

    BFSIter0(queue);
    queue.slide_window();   
    
    while(!queue.empty()){             
        //std::cout << "Queue not empty, Queue size: " << queue.size() << std::endl;
        pvector<bool> visited(this->N, false); 

        #pragma omp parallel
        {
            QueueBuffer<long int> lqueue(queue);
            #pragma omp for schedule(dynamic, 64)
            for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++){
                long int n = *q_iter;                        
                for(auto v: this->Neighbors[n]){
                    float curr_depth = this->property[v];
                    float new_depth = this->property[n] + 1;
                    if((new_depth < curr_depth) || (curr_depth == -1)){
                        bool curr_val = visited[v];
                        if(!curr_val){
                            if(compare_and_swap(visited[v], curr_val, true))
                                    lqueue.push_back(v);
                        }

                        while(!compare_and_swap(this->property[v], curr_depth, new_depth)){
                            curr_depth = this->property[v];
                            if(curr_depth <= new_depth){
                                break;
                            }
                        }
                    }               
                }
            }
            lqueue.flush();
        }
        queue.slide_window();               
    }    

    // clear affected array to get ready for the next update round
    #pragma omp parallel for schedule(dynamic, 64)
    for(long int i = 0; i < this->N; i++){
        this->affected[i] = false;
    }
}  

void Graph::BFSStartFromScratch(long int source) {  
    //std::cout << "Source " << source << std::endl;
    std::cout << "Running BFS from scratch" << std::endl;

		#pragma omp parallel for 
    for(long int n = 0; n < this->N; n++)
        this->property.push_back(-1);

    this->property[source] = 0;    
		cout << "here" << endl;
    
    SlidingQueue<long int> queue(this->N);   
    queue.push_back(source);
    queue.slide_window();  

    cout << "here" << endl;
    while(!queue.empty()){       
        //std::cout << "Queue not empty, Queue size: " << queue.size() << std::endl;         
        #pragma omp parallel
        {             
            QueueBuffer<long int> lqueue(queue);
            #pragma omp for 
            for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++){
                long int u = *q_iter;
                for(auto v: this->Neighbors[u]){
                    float curr_depth = this->property[v];
                    float new_depth = this->property[u] + 1;
                    if(curr_depth < 0){
                        if(compare_and_swap(this->property[v], curr_depth, new_depth)){
                            lqueue.push_back(v);
                        }
                    }
                }
            }
            lqueue.flush();
        }
        queue.slide_window();        
    }
}

// BFS verifier does a serial BFS from same source and asserts:
// - parent[source] = source
// - parent[v] = u  =>  depth[v] = depth[u] + 1 (except for source)
// - parent[v] = u  => there is edge from u to v
// - all vertices reachable from source have a parent
bool Graph::BFSVerifier(long int source,
                 const pvector<long int> &parent) {
  pvector<int> depth(this->N, -1);
  depth[source] = 0;
  vector<long int> to_visit;
  to_visit.reserve(this->N);
  to_visit.push_back(source);
  for (auto it = to_visit.begin(); it != to_visit.end(); it++) {
    long int u = *it;
    for (long int v : this->Neighbors[u]) {
      if (depth[v] == -1) {
        depth[v] = depth[u] + 1;
        to_visit.push_back(v);
      }
    }
  }
  for (long int u = 0; u < this->N; u++) {
    if ((depth[u] != -1) && (parent[u] != -1)) {
      if (u == source) {
        if (!((parent[u] == u) && (depth[u] == 0))) {
          cout << "Source wrong" << endl;
          return false;
        }
        continue;
      }
      bool parent_found = false;
      for (long int v : this->Parents[u]) {
        if (v == parent[u]) {
          if (depth[v] != depth[u] - 1) {
            cout << "Wrong depths for " << u << " & " << v << endl;
            return false;
          }
          parent_found = true;
          break;
        }
      }
      if (!parent_found) {
        cout << "Couldn't find edge from " << parent[u] << " to " << u << endl;
        return false;
      }
    } else if (depth[u] != parent[u]) {
      cout << "Reachability mismatch" << endl;
      return false;
    }
  }
  return true;
}

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <climits>
#include <algorithm>
using namespace std;

#define MILLION 1000000
#define TEN_MILLION 10000000
#define HUNDRED_MILLION 100000000

string input, output;
int64_t in_lines, num_edges;
int32_t num_vertices;

pair<int32_t, int32_t> *edge_list;
int64_t *offsets;
int32_t *raw_edges;

int main(int argc, char *argv[]) {
  if(argc < 3) exit(-1);

  input = argv[1];          // input file
  output = argv[2];  // output file

  std::ifstream infile(input);
  if (!infile.is_open()) {
    std::cout << "Couldn't open file " << input << std::endl;
    std::exit(-2);
  }

  int64_t line = 0;
  int64_t el = 0, ol = 0;
  int64_t u_, v_;
  int32_t u, v;
  int32_t mx_u = -(INT_MAX - 1), mn_u = INT_MAX;
  string line_str;
  while (infile >> line_str) {
    if(line == 0) {
      cout << "Dataset format: " << line_str << endl;
    }
    else if (line == 1) {
      num_vertices = stoi(line_str);
      cout << "num-vertices: " << num_vertices << endl;
      offsets = (int64_t *) malloc(sizeof(int64_t) * num_vertices);
    }
    else if (line == 2) {
      num_edges = stoll(line_str);
      cout << "num-edges: " << num_edges << endl;

      in_lines = num_edges * 2;
      raw_edges = (int32_t *) malloc(sizeof(int32_t) * num_edges);
      edge_list = (pair<int32_t, int32_t> *) malloc(sizeof(pair<int32_t, int32_t>) * in_lines);
    }
    else {
      v = stoi(line_str);
      if(line < num_vertices+3) {
        offsets[ol++] = v;
      }
      else {
        raw_edges[el++] = v;
      }
    }
    line += 1;
  }

  cout << "[DONE] reading the dataset." << endl;

  int64_t st, nd;
  int32_t max_degree = 0, max_vid = 0;
  line = 0;
  for(u=0; u<num_vertices; u+=1) {
    st = offsets[u];
    nd = (u == num_vertices - 1) ? num_edges : offsets[u + 1];
    if((nd - st) > max_degree) {
      max_degree = (nd - st);
      max_vid = u;
    }

    for(int64_t i=st; i<nd; i+=1) {
      v = raw_edges[i];
      edge_list[line++] = make_pair(u, v);
      edge_list[line++] = make_pair(v, u);
    }
    mx_u = max(u, mx_u);
    mn_u = min(u, mn_u);

    mx_u = max(v, mx_u);
    mn_u = min(v, mn_u);
  }

  cout << "[DONE] preparing the un-directed edge-list." << endl;
  cout << "MAX vertex-id: " << mx_u << ", MIN vertex-id: " << mn_u << endl;
  cout << "MAX degree: " << max_degree << ", max-degree vertex-id: " << max_vid << endl;

  free(raw_edges);
  cout << "[DONE] raw edges freed." << endl;

  sort(edge_list, edge_list+in_lines);
  cout << "[DONE] undirected edgelist sorted." << endl;

  int64_t out_line = 0;
  int32_t last_u = -1, last_v = -1;
  for(int64_t i=0; i<in_lines; i+=1) {
    // removing redundant edges
    if(edge_list[i].first == last_u && edge_list[i].second == last_v) continue;
    // removing self loops
    // if(edge_list[i].first == edge_list[i].second) continue;
    edge_list[out_line++] = edge_list[i];
    last_u = edge_list[i].first;
    last_v = edge_list[i].second;
  }
  cout << "[DONE] unique edge-list prepared. Total unique edges: " << out_line << endl;

  // obtain a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  shuffle(edge_list, edge_list + (out_line-1), std::default_random_engine(seed));
  cout << "[DONE] unique edge-list shuffled." << endl;

  ofstream out_file(output);
  if(!out_file) {
    cout << "Cannot open the unique-single-output file!" << endl;
    exit(-3);
  }

  for(int64_t i=0; i<out_line; i+=1) {
    out_file << edge_list[i].first << " " << edge_list[i].second << "\n";
  }
  cout << "[DONE] writing the unique-single-datasets." << endl;
  out_file.flush();
  out_file.close();

  return 0;
}

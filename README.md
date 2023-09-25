# Dynamic GAP
Question: What is GAP?
GAP is a static graph processing benchmark.

## GAP Benchmark Suite included.
We merged the GAP benchmark BFS implementation.
In order to use the BFS from GAp, we used the following header files from GAP implementations:
1. sliding_queue
2. bitmap
3. pvector
4. platform_atomics
For details, please visit https://github.com/sbeamer/gapbs/tree/master
Reference: Scott Beamer, Krste Asanović, and David Patterson. 
"Direction-Optimizing Breadth-First Search." 
International Conference on High Performance Computing, Networking, Storage and Analysis (SC), Salt Lake City, Utah, November 2012.

btree implementation is done by Google as part of the absl implementation.
Instructions: https://code.google.com/archive/p/cpp-btree/wikis/UsageInstructions.wiki
Number of edges will be limited to ~2B if we take 'int' to contain edge_list.
Hence, we used 'long int'.
Assumption: The nodes are marked from 0 onwards.

## Using graph converters from .adj to .el
g++ -o adj_to_el adj_to_el.cc -std=c++11 -O3
./adj_to_el INPUT_FILE OUTPUT_FILE

## Sanity check for GAP benchmark
The runtime matches between Dynamic GAP benchmark and the native GAP benchmark. Note that, Dynamic GAP benchmark uses individual btrees to store the in- and out-neighbors, whereas GAP uses CSR representation of a graph.
1. GAP's `src/` is clones in `src-gap/`.
2. Build `bfs.cc` by running `g++ -std=c++11 -O3 -Wall -fopenmp src-gap/bfs.cc`. Skip the option `-fopenmp` if not suppported on mac.
3. Run `./a.out -f graphs/slashdot.el`

## Including AlgoraCore
[Project page](https://libalgora.gitlab.io/#algora)
[AlgoraCore](https://gitlab.com/libAlgora/AlgoraCore) and [GitHub](https://github.com/libAlgora/AlgoraCore/tree/master)
[AlgoraDyn](https://gitlab.com/libAlgora/AlgoraDyn)


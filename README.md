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

Run on slashdot:
1. Build: `g++ -std=c++11 -fopenmp graph-btree-per-node.cpp`.
2. Run: `./a.out -deterministic 15462 graphs/slashdot.el 938360 0 0`.

## Sanity check for GAP benchmark
The runtime matches between Dynamic GAP benchmark and the native GAP benchmark. Note that, Dynamic GAP benchmark uses individual btrees to store the in- and out-neighbors, whereas GAP uses CSR representation of a graph.
1. GAP's `src/` is clones in `src-gap/`.
2. Build `bfs.cc` by running `g++ -std=c++11 -O3 -Wall -fopenmp src-gap/bfs.cc`. Skip the option `-fopenmp` if not suppported on mac.
3. Run `./a.out -f graphs/slashdot.el`

P.S. Build `g++ -std=c++11 test-sanity-GAP.cpp` and run `./a.out >> tmp`.
P.P.S. Run `make gap` to build all files in `src-gap`. Please remember to do `make clean` before pushing on git.

## Including AlgoraCore
[Project page](https://libalgora.gitlab.io/#algora) \
[AlgoraCore](https://gitlab.com/libAlgora/AlgoraCore) and [GitHub](https://github.com/libAlgora/AlgoraCore/tree/master) \
[AlgoraDyn](https://gitlab.com/libAlgora/AlgoraDyn)

1. Go to `Algora/AlgoraCore/`.
2. Run `./easyCompile` if not already run.
3. Go to `Algora/AlgoraCore/examples`.
4. Build `g++ -std=c++17 -Wall -I../src/ -L../build/Release/ bfs_algora.cpp -lAlgoraCore`.
5. Run `./a.out /home/arghya/DynamicGAP/graphs/slashdot.el 938360 144 > tmp`.

## Including Kickstarter
[Project page](https://github.com/pdclab/graphbolt/)
1. `sudo vim /etc/apt/sources.list` and Add `deb http://dk.archive.ubuntu.com/ubuntu/ xenial main` and `deb http://dk.archive.ubuntu.com/ubuntu/ xenial universe`.
2. `sudo apt update`.
3. `sudo apt install g++-5 gcc-5`.
4. `sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 50`.
5. `sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 60`.
6. `sudo update-alternatives --config gcc`.
7. `sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 5`.
8. `sudo update-alternatives --config g++`.
9. `g++ -v`.
10. `sudo chmod +x ./graphbolt/install_mimalloc.sh && ./graphbolt/install_mimalloc.sh`.
11. `export LD_PRELOAD=/home/arghya/DynamicGAP/graphbolt/lib/mimalloc/out/release/libmimalloc.so`.
12. `cd ./graphbolt/apps && make -j && cd ../..`.
13. `cd ./graphbolt/tools/converters/ && make && cd ../../..`.
14. `sudo find / -name libmimalloc.so`.
15. `echo $LD_LIBRARY_PATH`.
16. `LD_LIBRARY_PATH=/usr/local/lib`.
17. `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/arghya/DynamicGAP/graphbolt/lib/mimalloc/out/release/`.
18. `export LD_LIBRARY_PATH`.

Reference: [stack overflow](https://stackoverflow.com/questions/67280779/cilk-h-no-such-file-or-directory), 
[ask ubuntu] (https://askubuntu.com/questions/1235819/ubuntu-20-04-gcc-version-lower-than-gcc-7), and
[stack overflow](https://stackoverflow.com/questions/480764/linux-error-while-loading-shared-libraries-cannot-open-shared-object-file-no-s).


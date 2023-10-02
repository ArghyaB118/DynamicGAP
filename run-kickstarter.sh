#!/bin/bash

#sed -n '1,937300p' ~/DynamicGAP/graphs/slashdot.el > tmp1.el
sed -n '1,10000000p' ~/DynamicGAP/graphs/er_graph.el > tmp1.el
./graphbolt/tools/converters/SNAPtoAdjConverter tmp1.el tmp1.adj
sed -n '10000001,100000000p' ~/DynamicGAP/graphs/er_graph.el > tmp2.txt
sed -e 's/^/a /' tmp2.txt > tmp3.txt

#./graphbolt/apps/BFS -source 15462 -numberOfUpdateBatches 10 -nEdges 100 -streamPath tmp3.txt -outputFile bfs_output tmp1.adj
./graphbolt/apps/BFS -source 2001011 -numberOfUpdateBatches 2 -nEdges 10000000 -streamPath tmp3.txt -outputFile bfs_output tmp1.adj

rm tmp1.el tmp1.adj tmp2.txt tmp3.txt

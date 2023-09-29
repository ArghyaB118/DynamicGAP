#!/bin/bash

sed -n '1,937300p' ~/DynamicGAP/graphs/slashdot.el > tmp1.el
./graphbolt/tools/converters/SNAPtoAdjConverter tmp1.el tmp1.adj
sed -n '937301,938360p' ~/DynamicGAP/graphs/slashdot.el > tmp2.txt
sed -e 's/^/a /' tmp2.txt > tmp3.txt

./graphbolt/apps/BFS -source 15462 -numberOfUpdateBatches 10 -nEdges 100 -streamPath tmp3.txt -outputFile bfs_output tmp1.adj

rm tmp1.el tmp1.adj tmp2.txt tmp3.txt

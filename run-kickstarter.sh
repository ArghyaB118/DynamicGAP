#!/bin/bash
now=$(date)
echo $now

while getopts s:n:e:p: flag
do
	case "${flag}" in
        s) source=${OPTARG};;
        n) numberOfUpdateBatches=${OPTARG};;
        e) nEdges=${OPTARG};;
		p) streamPath=${OPTARG};;
    esac
done

echo $streamPath

#sed -n '1,185083p' $streamPath > tmp.el
#./graphbolt/tools/converters/SNAPtoAdjConverter tmp.el tmp.adj
#sed -n '5084,117370166p' $streamPath > tmp.txt
#sed -e 's/^/a /' tmp.txt > tmp2.txt
#./graphbolt/apps/BFS -source $source -numberOfUpdateBatches 1 -nEdges 117185083 -streamPath tmp2.txt -outputFile bfs_output tmp.adj

#sed -n '1,117370166p' $streamPath > tmp.el
#./graphbolt/tools/converters/SNAPtoAdjConverter tmp.el tmp.adj
#sed -n '117370167,234370166p' $streamPath > tmp.txt
#sed -e 's/^/a /' tmp.txt > tmp2.txt
#./graphbolt/apps/BFS -source $source -numberOfUpdateBatches $numberOfUpdateBatches -nEdges $nEdges -streamPath tmp2.txt -outputFile bfs_output tmp.adj

head -n "$(($(wc -l<$streamPath)/2))" $streamPath > tmp.el
./graphbolt/tools/converters/SNAPtoAdjConverter tmp.el tmp.adj
tail -n "$(($(wc -l<$streamPath)/2))" $streamPath > tmp.txt
sed -e 's/^/a /' tmp.txt > tmp2.txt
vmtouch -vt tmp.adj tmp2.txt
./graphbolt/apps/BFS -source $source -numberOfUpdateBatches $numberOfUpdateBatches -nEdges $nEdges -streamPath tmp2.txt -outputFile bfs_output tmp.adj

rm tmp.el tmp.adj tmp.txt tmp2.txt bfs_output*

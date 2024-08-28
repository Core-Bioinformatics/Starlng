
#!/bin/bash

dts_name=$1
ncells=$2
ncores=$3
ncores_graph_learn=$4

# Run the benchmarking script
bash memory_profiling.sh 5 ${dts_name}/${dts_name}_${ncells}_memory_${ncores}_cores_${ncores_graph_learn}_for_graph_learn.txt &
sleep 15
Rscript 2_benchmark_specific_ncells.R $dts_name $ncells $ncores $ncores_graph_learn
touch trigger.txt
sleep 15


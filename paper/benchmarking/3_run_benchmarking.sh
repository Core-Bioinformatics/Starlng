
#!/bin/bash

dts_name=$1
# ncores=$2
# ncores_graph_learn=$3

# Run the benchmarking script
for ncells in 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000
do
    bash memory_profiling.sh 5 ${dts_name}/${dts_name}_${ncells}_memory_${ncores}_cores_${ncores_graph_learn}_for_graph_learn.txt &
    sleep 15
    Rscript 2_benchmark_specific_ncells.R $dts_name $ncells 30 30
    touch trigger.txt
    sleep 15
done


#!/bin/bash

dts_name=$1
ncells=100000

# Run the benchmarking script
for ncores in 25 20 15 10 5 1
do
    bash memory_profiling.sh 5 ${dts_name}/${dts_name}_${ncells}_memory_${ncores}_cores_${ncores}_for_graph_learn.txt &
    sleep 15
    Rscript 2_benchmark_specific_ncells.R $dts_name $ncells $ncores $ncores
    echo "create trigger"
    touch trigger.txt
    sleep 15
done

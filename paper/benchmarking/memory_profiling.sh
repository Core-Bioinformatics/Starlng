#!/bin/bash

# two arguments: frequency (in seconds) and output file
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <frequency> <output_file>"
    exit 1
fi

freq=$1
output=$2

while true; do
    total_mem=$(cat /proc/meminfo | grep MemTotal | tr -d -c 0-9)
    free_mem=$(cat /proc/meminfo | grep MemFree | tr -d -c 0-9)
    used_mem=$(($total_mem - $free_mem))

    echo "$used_mem" >> $output
    sleep $freq

    if [ -f "trigger.txt" ]; then
        rm trigger.txt
	echo "done memory"
        break
    fi
done

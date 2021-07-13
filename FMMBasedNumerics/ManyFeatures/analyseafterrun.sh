#!/bin/bash
echo "Error messages":
cat LOGFILES/*err
echo "Number of times successfully finished":
cat LOGFILES/*out | grep "FMM ended." | wc -l
echo "Max memory used; used / max; be aware of potential 1000 / 1024 issue":
cd LOGFILES
sh makesummary.sh
cd ..
cat LOGFILES/summary_*.log | grep batch | awk '{ print $4, $3 }' | sort -h | tail -n 2

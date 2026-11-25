#!/bin/bash
cd /home/user/Mathathon/sphere_packing_research/rounds/R9_trace_lp_lattices
while pgrep -f run2.sh > /dev/null; do sleep 1; done
for n in 7 9 11 13; do timeout 25 python3 lp_scan.py $n gm 24 10 1 data/lpscan_n${n}_gm.json; done
timeout 40 python3 lp_scan.py 14 gm 30 20 1 data/lpscan_n14_gm.json
timeout 80 python3 lp_scan.py 22 gm 10 50 1 data/lpscan_n22_gm.json
timeout 80 python3 lp_scan.py 24 gm 6 50 1 data/lpscan_n24_gm.json
echo P4 DONE

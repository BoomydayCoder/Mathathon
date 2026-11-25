#!/bin/bash
cd /home/user/Mathathon/sphere_packing_research/rounds/R9_trace_lp_lattices
timeout 75 python3 lp_scan.py 26 gm 6 50 1 data/lpscan_n26_gm.json
timeout 85 python3 lp_scan.py 28 gm 4 55 1 data/lpscan_n28_gm.json
echo P5 DONE

#!/bin/bash
cd /home/user/Mathathon/sphere_packing_research/rounds/R9_trace_lp_lattices
timeout 70 python3 lp_scan.py 16 gm 30 45 1 data/lpscan_n16_gm.json
timeout 80 python3 lp_scan.py 18 gm 20 55 1 data/lpscan_n18_gm.json
timeout 90 python3 lp_scan.py 20 gm 12 60 1 data/lpscan_n20_gm.json
echo P3 DONE

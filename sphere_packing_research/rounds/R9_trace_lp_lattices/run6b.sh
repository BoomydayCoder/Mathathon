#!/bin/bash
cd /home/user/Mathathon/sphere_packing_research/rounds/R9_trace_lp_lattices
timeout 40 python3 lp_scan.py 16 gauss 14 18 1 data/lpscan_n16_gauss.json
timeout 50 python3 lp_scan.py 20 gauss 8 25 1 data/lpscan_n20_gauss.json
echo P6b DONE

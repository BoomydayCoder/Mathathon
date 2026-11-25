#!/bin/bash
cd /home/user/Mathathon/sphere_packing_research/rounds/R9_trace_lp_lattices
timeout 140 python3 tracelp.py 10 gm 24 110 1 24 data/n10_gm.json
timeout 115 python3 tracelp.py 10 gauss 20 85 1 20 data/n10_gauss.json
timeout 200 python3 tracelp.py 14 gm 10 160 1 2 data/n14_gm.json
echo P2 DONE

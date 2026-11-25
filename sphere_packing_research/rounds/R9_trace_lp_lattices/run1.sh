#!/bin/bash
cd /home/user/Mathathon/sphere_packing_research/rounds/R9_trace_lp_lattices
timeout 40 python3 tracelp.py 6 gm 24 20 1 24 data/n6_gm.json
timeout 40 python3 tracelp.py 6 gauss 24 20 1 24 data/n6_gauss.json
timeout 80 python3 tracelp.py 8 gm 24 50 1 24 data/n8_gm.json
timeout 80 python3 tracelp.py 8 gauss 24 50 1 24 data/n8_gauss.json
timeout 190 python3 tracelp.py 12 gm 20 150 1 6 data/n12_gm.json
timeout 120 python3 tracelp.py 12 gauss 12 90 1 4 data/n12_gauss.json
echo P1 DONE

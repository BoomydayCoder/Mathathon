#!/bin/bash
cd /home/user/Mathathon/sphere_packing_research/rounds/R9_trace_lp_lattices
GMP=1000003 timeout 30 python3 lp_scan.py 16 gm 14 14 1 data/lpscan_n16_gmP1e6.json
GMP=1000003 timeout 35 python3 lp_scan.py 20 gm 7 18 1 data/lpscan_n20_gmP1e6.json
echo P7 DONE

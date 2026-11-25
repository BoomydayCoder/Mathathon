#!/bin/bash
cd /home/user/Mathathon/sphere_packing_research/rounds/R9_trace_lp_lattices
GMP=90000049 timeout 40 python3 lp_scan.py 16 gm 14 18 1 data/lpscan_n16_gmP9e7.json
GMP=90000049 timeout 50 python3 lp_scan.py 20 gm 8 25 1 data/lpscan_n20_gmP9e7.json
echo P6a DONE

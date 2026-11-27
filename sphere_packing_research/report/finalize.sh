#!/bin/bash
# Final assembly: research log tail (from its source pieces) + report PDF.  Idempotent; does not stamp any time
# (the recorded end of the timed loop, 10:39:49 UTC, is fixed text in log_footer.md).
set -e
cd "$(dirname "$0")"
python3 assemble_log.py
./build.sh
cp main.pdf ../sphere_packing_report.pdf
pdfinfo ../sphere_packing_report.pdf | grep Pages

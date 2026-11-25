#!/bin/bash
# Final assembly: append the closing summary to the log (once), rebuild the report, copy the PDF next to the log.
set -e
cd "$(dirname "$0")"
LOG=../research_log.md
if ! grep -q "^## 종료 정리" "$LOG"; then
  cat final_summary.md >> "$LOG"
  printf '\n**기록 종료 시각**: %s UTC (시작 06:24:15 UTC 로부터 %s분 경과).\n' "$(date -u +%H:%M:%S)" "$(( ( $(date -u +%s) - 1795587855 ) / 60 ))" >> "$LOG"
fi
./build.sh
cp main.pdf ../sphere_packing_report.pdf
pdfinfo ../sphere_packing_report.pdf | grep Pages

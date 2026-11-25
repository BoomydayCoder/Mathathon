#!/bin/bash
# Final assembly (idempotent): replace the closing summary at the end of the log by the current final_summary.md,
# rebuild the report, copy the PDF next to the log.
set -e
cd "$(dirname "$0")"
python3 - <<'PY'
import time, calendar
LOG = '../research_log.md'
t = open(LOG, encoding='utf-8').read()
marker = '\n---\n\n## 종료 정리'
i = t.find(marker)
if i >= 0:
    t = t[:i]
t = t.rstrip('\n') + '\n' + open('final_summary.md', encoding='utf-8').read().rstrip('\n') + '\n'
now = time.gmtime()
elapsed = (calendar.timegm(now) - 1795587855) // 60
t += '\n**기록 종료 시각**: %s UTC (시작 06:24:15 UTC 로부터 %d분 경과; 5시간 제한 이내).\n' % (time.strftime('%H:%M:%S', now), elapsed)
open(LOG, 'w', encoding='utf-8').write(t)
print('log finalized, elapsed', elapsed, 'min')
PY
./build.sh
cp main.pdf ../sphere_packing_report.pdf
pdfinfo ../sphere_packing_report.pdf | grep Pages

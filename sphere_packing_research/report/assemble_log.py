#!/usr/bin/env python3
"""Re-assemble the tail of ../research_log.md from its source pieces (idempotent; never stamps a time).

  body (rounds 0-19, untouched)  +  log_r20_r21.md (post-record rounds)  +  final_summary.md  +  log_footer.md
"""
import os
here = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(here, '..', 'research_log.md')
rd = lambda f: open(os.path.join(here, f), encoding='utf-8').read()
t = open(LOG, encoding='utf-8').read()
cuts = [i for i in (t.find('\n---\n\n## 20회차 (기록 종료 후'), t.find('\n---\n\n## 종료 정리')) if i >= 0]
assert cuts, 'markers not found'
body = t[:min(cuts)].rstrip('\n') + '\n'
old = ('- **읽는 법**: 이 로그는 종료된 기록이다(0–19회차, 기록 종료 10:39:49 UTC, 경과 255분).')
new = ('- **읽는 법**: 이 로그는 시간 제한이 있던 연구 루프의 기록(0–19회차, 기록 종료 10:39:49 UTC, 경과 255분)에, '
       '그 뒤 사용자의 추가 요청으로 수행한 20·21회차(11월 25일, 27일; 5시간 루프 밖)를 덧붙인 것이다.')
if old in body:
    body = body.replace(old, new)
assert new in body, 'header line not found'
out = body + rd('log_r20_r21.md').rstrip('\n') + '\n' + rd('final_summary.md').rstrip('\n') + '\n\n' + rd('log_footer.md').rstrip('\n') + '\n'
open(LOG, 'w', encoding='utf-8').write(out)
print('research_log.md assembled:', out.count('\n'), 'lines')

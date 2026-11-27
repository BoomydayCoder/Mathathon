# 고차원 격자 구 채우기 하한 $c\,n^2\log\log n\cdot2^{-n}$ 을 넘어서려는 시도 — 연구 기록

Klartag (2025, arXiv:2504.05042)의 격자 구 채우기 하한 $\delta_n\ge c\,n^2\,2^{-n}$ 과 그 후속 결과 $c\,n^2\log\log n\cdot2^{-n}$
(Abuya–Gargava–Zhao, Gargava–Gröbner–Simkin, 2026)을 개선하려고 한 자율 연구 세션의 기록이다.

**결론부터: 차수 $n^2\log\log n$ 을 개선하는 정리는 얻지 못했다.** 얻은 것은 방법의 한계에 대한 정리, 명시적 상수의 개선
($n^2$ 차수 안에서), 유일하게 살아남은 후보 족에 대한 부분 결과와 계산, 그리고 실패한 접근의 기록이다.

> **검증 수준에 대한 주의.** 이 저장소의 "증명됨" 은 모두 *AI 에이전트가 쓴 증명을 다른 AI 에이전트가 반박하려고 검토한*
> 수준이며 **사람의 검증은 전혀 거치지 않았다.** 기록 종료 후의 20·21회차는 그 AI 심사조차 없다.

## 읽을 파일

| 파일 | 내용 |
|---|---|
| [`sphere_packing_report.pdf`](sphere_packing_report.pdf) | 정리된 보고서(143쪽): 초록, 결과 요약표(§2), 본문(§3–§9), 증명 전문 부록(A–D), 연구 로그 전문(E), 파일 색인(F), 회차 노트 원문(G) |
| [`research_log.md`](research_log.md) | 회차별 연구 로그(0–21회차)와 맨 끝의 "종료 정리"(확실한 결과 / 남은 빈틈 / 다음 과제) |
| `README.md` | 이 파일 |

보고서를 처음 읽는다면 초록 → §2 결과 요약표 → §9 "실패한 접근, 남은 빈틈, 다음 연구 과제" 순서를 권한다. 20·21회차의 내용은 §7.10, §7.11 에 있다.

## 진행 경과

- **시간 제한 루프 (0–19회차)**: 2026-11-25 06:24:15 UTC 시작, 10:39:49 UTC 기록 종료(255분; 5시간 제한 이내). 회차마다
  "과제 선택 → 증명 시도 → 반례·수치 검증 → 기록" 을 반복했고, 대부분의 회차는 증명/계산 담당 에이전트와 그것을 반박하도록
  지시받은 심사 에이전트의 쌍으로 구성했다.
- **기록 종료 후의 추가 회차** (사용자의 추가 질문에 따른 것; 5시간 루프에 포함되지 않음)
  - 20회차 (11-25 18:25–19:55 UTC): "$(p+1)/2$ 가 소수인 $p$ 가 무한히 많으면 안정한 $p$ 도 무한히 많은가" 의 증명 시도.
    하위 에이전트 하나, 심사 없음, 오케스트레이터의 부분 점검.
  - 21회차 (11-27): 가정 "$(p+1)/2$ 가 소수" 를 약화할 수 있는가. 오케스트레이터의 직접 유도와 GAP 계산, 정식 심사 없음.
  - 11-27: 두 회차를 세 문서에 반영한 뒤 비평 에이전트 하나가 원자료와 대조했고, 지적 19건을 반영해 개정했다.

## 결과 요약

표기: $V_n$ = (공부피 1 인 격자에서 원점 외의 격자점이 없는 가장 큰 열린 공의 부피)의 상한 $=2^n\delta_n$.
알려진 최선: $V_n\ge c\,n^2\log\log n$; 명시적 상수는 Rezvyakova 의 $V_n\ge(\tfrac1{2e}-o(1))n^2$.

### 얻은 것

1. **이득 지표와 자유 작용의 한계** [증명됨; Zassenhaus 구조정리 인용]. 대칭을 쓰는 Klartag 형 논증의 이득은
   $g=|G|/(2\dim_{\mathbb Q}A)=\tfrac12\sum_g\rho(g)^2$ 이고, 고정점 없는 작용에서는 $g=O(\log\log|G|)$ 다.
   (방법의 회계량에 대한 정리이지 모든 논증에 대한 장벽 정리는 아니다.)
2. **명시적 상수 (차수는 $n^2$)** [증명됨·인용 정리 전제]. Hurwitz 가군: $V_n\ge(\tfrac3{2e}-o(1))n^2$ ($4\mid n$);
   icosian 가군: $V_n\ge(\tfrac{15}{4e}-o(1))n^2$ ($8\mid n$).
3. **비자유 작용의 틀 정리와 Leech 텐서족** [증명됨·인용 정리 전제 + 계산기]. 고정된 격자 $\Lambda_W$ 가 극대 order 위에서
   "경계 안정" ($R_U\le1$)이면 텐서족 $\Lambda_W\otimes_{\mathcal O}L$ 에서 $V_n\ge(1-o(1))\tfrac{|G|}{4e\dim_{\mathbb Q}A}n^2$.
   Leech 격자의 $SL_2(13)$ 대칭에 적용하면 **$V_n\ge(\tfrac{91}{12e}-o(1))n^2\approx2.79\,n^2$ ($24\mid n$)** — Rezvyakova 상수의
   15.2배. 순수 점근 진술이다(증명 그대로는 $n\gtrsim10^{14}$ 에서야 비자명).
4. **유일하게 살아남은 후보: $SL_2(p)$ Weil 격자 $W_p$** ($p\equiv1\bmod4$, $(p+1)/2$ 소수; 랭크 $2(p-1)$, 이득 $\approx p$).
   - $p=13$ (Leech): 정확히 경계 안정 [증명됨(계산기)]. $p=37$: **불안정** [증명됨(계산기; 정수 인증서)].
     $p=61$: 수행한 모든 구조적 시험을 통과했으나 미결정.
   - 표수 $p$ 의 곡선 $y^2=x^p-x$ 의 Mordell–Weil 격자와의 사전, 몫곡선 사전 [증명됨·인용 전제].
   - (20회차, 미심사) $p\equiv5\bmod8$ 에서 $W_p$ 의 족은 원분체 $\mathbb Q(\zeta_{4p})$ 의 명시적인 ideal 격자로 주어진다:
     짝 unimodular 인 것은 $L_1=(\mathfrak P_1^{(m+1)/2}\mathfrak P_2^{-(m-1)/2},\ \mathrm{Tr}(x\bar x)/2p)$ 와
     $L_0=(\mathfrak P_1^{-(3m-1)/2}\mathfrak P_2^{-(m-1)/2},\ \mathrm{Tr}(\varepsilon x\bar x/2\sqrt p))$ 둘($m=(p-1)/4$; 둘의 동형은 미증명).
     벡터 수준의 안정성 조건은 $L_1$ 에서 "모든 $x\ne0$ 에서 $T_+(x)T_-(x)\ge p(p-1)^2$", $L_0$ 에서 "$T_+T_-\ge(p-1)^2$" 과
     동치이고, 뒤의 것은 $L_0$ 의 지수 $p^{(p-3)/2}$ 인 부분격자 $\mathbb Z[\zeta_{4p}]$ 위에서 성립한다. 정확한 Gram 행렬로 랭크 24
     부분격자를 $p<8900$ 의 85개 소수에서 전수 탐색: 위반은 $p=37$ 뿐.
5. **모형과 소차원 실험** [모형에서 증명 / 수치관찰]. Poisson 모형에서 대각합-LP 최적 타원체의 부피 $\sim\chi^2_{2N}$;
   등변 Klartag 과정은 $E_8,K_{12},BW_{16},\Lambda_{20}$, Leech 를 재발견(새 기록 없음).

### 얻지 못한 것

- 차수 개선. 조건부 경로($V_n\gtrsim n^{5/2}/\sqrt{\log n}$, 희소 차원열)는 (i) 무한히 많은 $p$ 에서 $W_p$ 의 안정성,
  (ii) $p$ 와 랭크가 함께 커질 때의 균등한 정리, (iii) 그런 소수가 무한히 많다는 것을 모두 요구하며 셋 다 미증명이다.
- "안정한 $p$ 가 무한히 많다" (20회차): 증명도 반증도 되지 않았다. 증명된 하한은 필요한 값의 약 $1/p$ 이고, 비주 직선과
  고차원 부분공간에 대해서는 진전이 없다. $p\equiv1\bmod8$ 은 모형조차 없다.
- 가정 "$(p+1)/2$ 가 소수" 의 약화 (21회차; unimodular 격자에 대한 필요조건): $SL_2$ 족 안에서는 불가능하다 —
  $(p+1)/2$ 가 합성수이면(서로 다른 두 소수의 제곱으로 동시에 나누어지는 경우만 빼고) 어떤 충실한 유리 표현에서도 비분리
  토러스가 이득을 주지 못하고, 홀 이산계열의 합에서는 이득이 약 2 이하다. 그 밖의 군에서는 방법을 찾지 못했다: GAP 지표표
  라이브러리 2313개의 스캔에서 필요조건을 통과하며 icosian 값(7.5)을 넘는 것은 $SL_2(13)$ 뿐이었다(분류 정리는 아니다). 다만 그런 소수는 드물지 않다(2000 미만에 34개). 상수 개선만이라면 안정한 소수
  하나로 충분하다.

## 디렉터리

```
sphere_packing_research/
├── README.md
├── research_log.md              연구 로그 (0–21회차 + 종료 정리)
├── sphere_packing_report.pdf    보고서
├── report/                      보고서의 LaTeX 원본과 빌드 스크립트
│   ├── main.tex, sec_*.tex, preamble.tex
│   ├── sec_post.tex             20·21회차 절
│   ├── final_summary.md, log_r20_r21.md, log_footer.md   로그 끝부분의 원본 조각
│   ├── assemble_log.py          로그 끝부분 재조립 (시간을 새로 찍지 않는다)
│   └── build.sh, build_notes.sh, finalize.sh
└── rounds/                      회차별 코드·노트·심사 보고
    ├── R1_* … R19_classnumber/          시간 제한 루프
    ├── R20_infinitely_many_stable/      20회차 (proof.md, ideal_model.py, scan_exact.py, orchestrator_check/)
    └── R21_weaker_conditions/           21회차 (notes.md, scan.g, scan3.g, scan4.g, sl2_family.py, composite_M*.py)
```

`papers/` (내려받은 논문)는 저장소에 넣지 않았다.

## 재현 방법

필요한 것: Python 3 (numpy, scipy, sympy), gcc. 보고서 빌드에는 XeLaTeX 과 pandoc, 21회차 스캔에는 GAP 4.12 + CTblLib.

```bash
# 보고서와 로그 끝부분 다시 만들기
cd report && ./finalize.sh
```

```bash
# 구조적 부분격자 스캔 (13회차): p = 37 의 위반(최소 Gamma/theta = 7/9)과 p = 61 의 통과
cd rounds/R10_weil_lattices/orchestrator && python3 scan.py 37 61
```

```bash
# 20회차: ideal 격자 모형의 독립 점검 (랭크 24 격자의 근; 인자 없으면 p = 13, 37)
cd rounds/R20_infinitely_many_stable && python3 orchestrator_check/chk.py 13 37 61
```

```bash
# 21회차: SL_2(q) 족, 합성수 (p+1)/2, 지표표 스캔
cd rounds/R21_weaker_conditions && python3 sl2_family.py && python3 composite_M.py && gap -q -b scan4.g
```

주의: `fp_enum.c` (R2, R8, R10/A)는 좌표 버퍼가 64개라 랭크 64 이하에서만 써야 한다. 더 큰 랭크에는
`rounds/R20_infinitely_many_stable/fp_enum_big.c` 나 `rounds/R12_p61/lat.c` 를 쓴다.

## 지위 표기

[증명됨] / [증명됨(계산기)] / [증명됨·인용 정리 전제] / [모형에서 증명] / [조건부] / [추측] / [휴리스틱] / [수치관찰] / [실패].
열람하지 못한 문헌(Gross 1990, Nebe 1998, Dummigan, Tiep, Weil 의 환원 이론 등)에 의존하는 곳은 보고서 본문에 표시했다.
새로움은 어디서도 단정하지 않는다(제한된 문헌 검색에서 "찾지 못함" 수준).

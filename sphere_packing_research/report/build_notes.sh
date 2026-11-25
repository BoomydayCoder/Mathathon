#!/bin/bash
# Appendix with the key per-round working notes (plain-text English notes are included verbatim).
set -e
cd "$(dirname "$0")"
out=app_notes.tex
: > $out
add() {
  f="$1"; title="$2"
  [ -f "$f" ] || return 0
  echo "\\subsection{$title}" >> $out
  echo "{\\footnotesize\\noindent 원본: \\path{${f#../}}\\par}" >> $out
  mkdir -p notes_tmp
  g="notes_tmp/$(echo "${f#../rounds/}" | tr '/' '_')"
  sed -e 's/𝒫/P/g; s/∖/\\/g; s/⟹/==>/g; s/⟺/<=>/g; s/∤/ does not divide /g' "$f" > "$g"
  echo "\\VerbatimInput[breaklines=true,breakanywhere=true,fontsize=\\scriptsize,formatcom=\\notesfont]{$g}" >> $out
  echo "" >> $out
}
add ../rounds/R3_G_klartag/notes.md "3회차: 이득 항등식과 자유 작용의 한계 (증명 노트)"
add ../rounds/R3_non_free/notes.md "3회차: 비자유 작용, 텐서 안정성, no-go (증명 노트)"
add ../rounds/R3_process/notes.md "3회차: 과정의 한계, Lemma N (증명 노트)"
add ../rounds/R4_perceptron_theory/referee/notes.md "4회차: 정확한 \$\\chi^2\$ 법칙 등 (심사자 노트)"
add ../rounds/R5_leech13/notes.md "5회차: Leech 격자의 계산 (노트)"
add ../rounds/R5_sl2p_theory/notes.md "5회차: 일반 \$p\$ 의 보조정리 (증명 노트)"
add ../rounds/R7_gluing/notes.md "7회차: 원분 접착 (노트)"
add ../rounds/R10_weil_lattices/A/notes.md "10회차: Weil 격자 --- 빌더 A 노트"
add ../rounds/R10_weil_lattices/B/notes.md "10회차: Weil 격자 --- 빌더 B 노트"
add ../rounds/R11_mordell_weil/T1/notes.md "11회차: Mordell--Weil 재정식화의 증명 (T1 노트)"
add ../rounds/R11_mordell_weil/T2/notes.md "11회차: 조건부 주장에 대한 반박 시도 (T2 노트)"
add ../rounds/R12_p61/notes.md "12회차: \$p=61\$ 의 체계적 시험 (노트)"
add ../rounds/R12_p61/verify/notes.md "12회차: 독립 검증 노트"
add ../rounds/R14_leech_tensor/novelty_check.md "14회차: 문헌 점검 보고"
add ../rounds/R14_leech_tensor/referee_A.md "14회차: 심사 보고 A"
add ../rounds/R14_leech_tensor/referee_B.md "14회차: 심사 보고 B"
add ../rounds/R15_rogers/referee.md "15회차: 심사 보고 A"
add ../rounds/R15_rogers/referee_B.md "15회차: 심사 보고 B"
add ../rounds/R16_quotient_curve/referee.md "16회차: 몫곡선 사전의 심사 보고"
add ../rounds/R17_gargava_check/notes.md "17회차: Gargava 정리 16 의 점검 노트"
add ../rounds/R17_gargava_check/referee.md "17회차: 심사 보고"
add ../rounds/R18_integration/report.md "18회차: 14·15회차 결합 점검 보고 (통합 진술 포함)"
add ../rounds/R19_classnumber/referee.md "19회차: 심사 보고 A"
add ../rounds/R19_classnumber/referee_B.md "19회차: 심사 보고 B (독립 심사)"

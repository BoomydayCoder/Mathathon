"""Static figures for the report (print, light surface). Palette: validated categorical slots 1-3
(#257de0, #f47230, #00b894) of the reference palette; text in text inks; direct labels + legend."""
import json, math, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
fm.fontManager.addfont('/usr/share/fonts/truetype/nanum/NanumGothic.ttf')
plt.rcParams.update({'font.family': 'NanumGothic', 'axes.unicode_minus': False, 'font.size': 9,
                     'axes.edgecolor': '#999999', 'axes.linewidth': 0.6, 'text.color': '#191919',
                     'axes.labelcolor': '#191919', 'xtick.color': '#666666', 'ytick.color': '#666666',
                     'mathtext.fontset': 'cm'})
C1, C2, C3 = '#257de0', '#f47230', '#00b894'
INK, INK2, GRID = '#191919', '#666666', '#e6e6e6'

# ---------- Figure 1: V/N versus dimension on genuine random lattices ----------
nK = [2, 3, 4, 5, 6, 7, 8, 9, 10]
vK = [1.21, 0.99, 0.93, 0.97, 1.08, 1.26, 1.24, 1.01, 1.03]
nL = [6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24, 28]
vL = [1.12, 1.33, 1.64, 1.166, 1.17, 1.162, 1.18, 1.195, 1.198, 1.276, 1.350, 1.441, 1.545, 1.550, 1.752]
fig, ax = plt.subplots(figsize=(6.3, 3.3), dpi=200)
ax.set_facecolor('white'); fig.patch.set_facecolor('white')
ax.grid(axis='y', color=GRID, linewidth=0.6); ax.set_axisbelow(True)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.axhline(2.0, color=INK2, linewidth=0.9, linestyle=(0, (4, 3)))
ax.axhline(1.0, color=INK2, linewidth=0.9, linestyle=(0, (1, 2)))
ax.text(2.1, 2.03, 'Poisson 모형의 대각합-LP 값 $V/N=2$', color=INK2, fontsize=8, va='bottom')
ax.text(11.5, 0.97, 'Poisson 모형의 sticky 과정 $V/N\\approx1$', color=INK2, fontsize=8, va='top')
ax.plot(nL, vL, color=C1, linewidth=2, marker='o', markersize=4.5, markeredgecolor='white', markeredgewidth=0.8, label='대각합-LP 꼭짓점 (9회차)')
ax.plot(nK, vK, color=C2, linewidth=2, marker='s', markersize=4.5, markeredgecolor='white', markeredgewidth=0.8, label='Klartag 과정을 동결까지 (2회차)')
ax.text(28.3, 1.752, '대각합-LP', color=INK, fontsize=8.5, va='center')
ax.text(10.35, 1.06, 'Klartag 과정', color=INK, fontsize=8.5, va='bottom')
ax.annotate('$E_8$', (8, 1.64), xytext=(9.3, 1.74), fontsize=8.5, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, linewidth=0.6))
ax.set_xlim(1.5, 31.5); ax.set_ylim(0.8, 2.2)
ax.set_xlabel('차원 $n$'); ax.set_ylabel('$V/N$   ($N=n(n+1)/2$)')
ax.set_title('실제 랜덤 격자에서의 부피: 접근 가능한 차원에서는 모형값 2 에 못 미친다', fontsize=9.5, loc='left', color=INK)
ax.legend(loc='lower right', bbox_to_anchor=(1.0, 0.16), frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig('figs/fig_VN.pdf'); fig.savefig('figs/fig_VN.png'); plt.close(fig)

# ---------- Figure 2: profile of the twisted minimum over the circle of invariant forms ----------
d = json.load(open('../rounds/R10_weil_lattices/orchestrator/figdata_tcircle.json'))
fig, ax = plt.subplots(figsize=(6.3, 3.3), dpi=200)
ax.set_facecolor('white'); fig.patch.set_facecolor('white')
ax.grid(axis='y', color=GRID, linewidth=0.6); ax.set_axisbelow(True)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.axhline(1.0, color=INK2, linewidth=0.9, linestyle=(0, (4, 3)))
ax.text(-0.99, 1.012, '안정성 문턱', color=INK2, fontsize=8, va='bottom')
series = [('13', C1, '$p=13$: Leech 전체 (랭크 24)', 'o'), ('37', C2, '$p=37$: $C_3$-고정 부분격자 (랭크 24)', 's'), ('61', C3, '$p=61$: $C_5$-고정 부분격자 (랭크 24)', '^')]
for k, c, lab, mk in series:
    ts, v = d[k]
    ax.plot(ts, v, color=c, linewidth=2, label=lab)
ts, v = d['37']; i = int(np.nanargmin(v))
ax.annotate('위반: $\\sqrt{7/9}=0.882$', (ts[i], v[i]), xytext=(-0.99, 0.795), fontsize=8.5, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, linewidth=0.6))
ax.text(0.60, 1.055, '$p=13$', color=INK, fontsize=8.5)
ax.text(-0.56, 1.56, '$p=61$', color=INK, fontsize=8.5)
ax.text(0.74, 0.915, '$p=37$', color=INK, fontsize=8.5)
ax.set_xlim(-1.02, 1.02); ax.set_ylim(0.76, 1.72)
ax.set_xlabel('$t/\\log\\varepsilon_p$   ($t=0$: Mordell–Weil 점,  $t=\\pm\\log\\varepsilon_p$: unimodular 점)')
ax.set_ylabel('$\\min_w Q_t(w)\\,/\\,2\\sqrt{\\theta_pD}$')
ax.set_title('$G$-불변 형식의 원 위에서의 최소값: 1 아래로 내려가면 텐서 불안정', fontsize=9.5, loc='left', color=INK)
ax.legend(loc='lower center', frameon=False, fontsize=8, ncol=1, bbox_to_anchor=(0.5, 0.0))
fig.tight_layout(); fig.savefig('figs/fig_tcircle.pdf'); fig.savefig('figs/fig_tcircle.png'); plt.close(fig)
print('figures written')

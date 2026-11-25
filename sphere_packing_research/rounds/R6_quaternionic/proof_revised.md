# R6 (REVISED). Klartag's stochastic ellipsoid inside the invariant forms of a division algebra: explicit constants

> **STATUS BOX** (revision of `proof.md` after two referee reports; both verdicts: *correct modulo minor fixes*, no mathematical error found; the reviser re-derived Lemma 0.2, Lemma 1.3, Prop. 2.4, Lemma 3.3, Lemma 3.4, Prop. 3.5, Lemma 4.2 and the bookkeeping of §5 by hand once more).
>
> * **Theorem A (Hurwitz, $n=4k$, constant $\frac3{2e}$): PROVED, with cited inputs (L1)-(L3).**
> * **Theorem B (icosians, $n=8k$, constant $\frac{15}{4e}$): PROVED, with cited inputs (L1)-(L3)** (see the caveat on (L1) for $D_\mathbb R=\mathbb H\times\mathbb H$ in Remark 4.1a).
> * **Theorem C (general template): PROVED (with cited inputs (L1)-(L3)) in the case $D_\mathbb R\cong\mathbb K^s$** (all archimedean completions are division algebras) -- this is **Theorem C′** below; it contains A and B and was checked by both referees. **General case ($D_{\mathbb R}$ with matrix blocks $M_m(\mathbb K)$, $m>1$): PROVED UP TO THE FOLLOWING GAP:** the two missing ingredients (formula/minimum of $q$; a Wishart-type angular lemma) are supplied with proofs in the new §7 of this revision (Lemma 7.1, Lemma 7.2), sanity-checked numerically (`check_blocks.py`), but §7 was written during the revision and has **not been refereed**; it must be regarded as unverified until an independent check. No hypothesis "positive involution" is needed in §7 (finite groups are unitarisable at the real places).
>
> **Cited inputs.** (L1) Gargava, arXiv:2107.04844, Theorem 16 (Siegel mean value theorem for $g\mathcal O^k$, $k\ge2$, $\mathcal O$ an order [G, Def. 1] in a $\mathbb Q$-division algebra, $f\in C_c$): used as a black box; the quoted statement was compared with the text of [G] by both referees and by the reviser; its proof was NOT re-verified. All three theorems are conditional on it. (L2) Klartag, arXiv:2504.05042: the deterministic Lemma 2.1 and Lemma A.1 are used as stated; the proofs of Lemma 2.2, Prop. 2.3, Cor. 2.4, Lemma 3.1, Cor. 3.2, Lemma 3.3, Prop. 4.1, Prop. 4.2 are used with $F_A$ replaced by the smaller space (13$_\mathcal Q$); the properties of $F_A,\pi_A$ that these proofs use are listed and verified in §1-§3 (tag [V]); referee 1 re-read each of these proofs with the substitution and confirmed. (L3) Standard tools: Itô formula, Lévy characterisation, Dambis-Dubins-Schwarz, reflection principle, the GOE operator-norm tail bound (Vershynin, HDP, Cor. 4.4.8) as used in [K, Cor. 3.2].
>
> **Nature of the statements.** Purely asymptotic along $n=4k$ (A), $n=8k$ (B), $n=k\delta$ (C′), $k\to\infty$ with $(D,G)$ fixed; the $o(1)$ is not effective and decays slowly (§6.6). Nothing is claimed for other $n$. "Factor $3$ / $7.5$ over Rezvyakova" compares leading constants only. The bounds do NOT improve the order of magnitude: they are weaker than the $c\,n^2\log\log n$ bounds of Abuya-Gargava-Zhao and Gargava-Groebner-Simkin for large $n$ (those have an unspecified small $c$); what is gained is an explicit leading constant.

Status tags: **[P]** proved here in full; **[V]** "verbatim": the cited proof of Klartag applies word for word once the listed
properties are checked (they are checked here); **[L]** quoted from the literature with exact reference and hypotheses; **[GAP]** not proved.
References: [K] Klartag, arXiv:2504.05042v2 (numbering of v2); [R] Rezvyakova, arXiv:2608.22590; [G] Gargava, *Lattice packings
through division algebras*, arXiv:2107.04844; [V] Vance, arXiv:1105.3779. Numerical sanity checks: `check_q.py` (same directory).

Normalisation: for a lattice $L\subset\mathbb R^n$, $V(L)=\sup\{\mathrm{Vol}(B)/\mathrm{covol}(L)\}$ over open origin-centred balls $B$ with
$B\cap L=\{0\}$; $V_n=\sup_L V(L)$; packing density $=V/2^n$. $\kappa_n=\mathrm{Vol}(B^n)$.

## Main results

Let $\mathbb K\in\{\mathbb R,\mathbb C,\mathbb H\}$, $d=\dim_{\mathbb R}\mathbb K\in\{1,2,4\}$.

**Standing data (D).** $D$ is a division algebra of dimension $\delta$ over $\mathbb Q$ such that $D_{\mathbb R}:=D\otimes_{\mathbb Q}\mathbb R\cong\mathbb K^{s}$
as $\mathbb R$-algebras ($\delta=sd$), and $G\subset D^\times$ is a finite subgroup, $w:=|G|$, whose $\mathbb Q$-span is $D$.
Put $\mathcal O:=\mathbb Z[G]$ (the $\mathbb Z$-span of $G$), and
$$D_k:=s\Big(k+\frac{d\,k(k-1)}2\Big)=\dim_{\mathbb R}\{(H_v)_{v\le s}:H_v\in \mathrm{Herm}_k(\mathbb K)\}.$$

**Theorem C′ (template, product-of-division-algebras case) [P modulo L].** For data (D), as $k\to\infty$ ($D,G$ fixed), with $n=k\delta$,
$$V_{n}\ \ge\ (1-o(1))\,\frac{w\,\delta\,k^2}{4e}\;=\;(1-o(1))\,\frac{w\,n^2}{4e\,\delta}\;=\;(1-o(1))\,\frac{w\,D_k}{2e}.$$
Moreover the lattice can be taken to be a $G$-invariant lattice $\lambda g\mathcal O^k$, $g\in GL_k(D_{\mathbb R})$.

**Theorem A (Hurwitz) [P modulo L].** $D=(\tfrac{-1,-1}{\mathbb Q})$, $G=2T$ = the 24 Hurwitz units, $\mathcal O$ = Hurwitz order, $s=1$, $\mathbb K=\mathbb H$,
$\delta=4$, $D_k=k(2k-1)$: for $n=4k\to\infty$,
$$V_n\ \ge\ (1-o(1))\frac{24k^2}{e}=\Big(\frac{3}{2e}-o(1)\Big)n^2 .$$

**Theorem B (icosians) [P modulo L].** $D=(\tfrac{-1,-1}{\mathbb Q(\sqrt5)})$, $G=2I$ (binary icosahedral group, $w=120$), $\mathcal O=\mathbb Z[2I]$ (icosian ring),
$s=2$, $\mathbb K=\mathbb H$, $\delta=8$, $D_k=2k(2k-1)$: for $n=8k\to\infty$,
$$V_n\ \ge\ (1-o(1))\frac{240k^2}{e}=\Big(\frac{15}{4e}-o(1)\Big)n^2 .$$

(Compare: [R] gives $(\frac1{2e}-o(1))n^2$ for all $n$; A and B improve this explicit leading constant by the factors $3$ and $7.5$, along $4\mid n$ resp. $8\mid n$ only, and only asymptotically.)

**Relation to the literature (novelty NOT verified line by line).** The mechanism -- run Klartag's process inside the forms invariant under a finite group acting freely, so that an orbit of $w$ contact points costs one linear constraint -- is that of Abuya-Gargava-Zhao (arXiv:2606.05105) and Gargava-Groebner-Simkin (arXiv:2610.12619) for cyclic groups of roots of unity (cyclotomic fields), where it yields $c\,n^2\log\log n$ with unspecified $c$. What this note adds, to the best of our knowledge: (i) the non-commutative orders (Hurwitz, icosians; random lattices from [G, Thm 16]); (ii) explicit leading constants $\frac{w}{4e\delta}$, via the identity $P(x\otimes x)=L_{xx^*/d}$, i.e. the projected variance $q(x)\ge q_0=\frac1\delta$ (Lemma 0.2), with $q\equiv\frac14$ for Hurwitz; (iii) the angular Lemma 3.4 handling non-constant $q$ ($s\ge2$). Neither referee nor the reviser checked whether those two papers already contain (i)-(iii).

**Literature inputs [L]** (everything else is proved below): (i) Gargava's mean value theorem [G, Thm 16] (Section 4);
(ii) Klartag's deterministic Lemma 2.1 and Lemma A.1 [K], the GOE tail bound used in [K, Cor. 3.2] (Vershynin, HDP, Cor. 4.4.8),
and standard stochastic calculus (Itô formula, Lévy characterisation, Dambis–Dubins–Schwarz, reflection principle) exactly as used in [K];
(iii) for Theorem B: the 120 unit icosians form a group (Conway–Sloane, SPLAG, Ch. 8 §2.1). This is also verified here by EXACT arithmetic in $\mathbb Q(\sqrt5)$ (`check_q.py`, part (2b): the monoid generated by the three generators of 0.1 has exactly 120 elements, all of reduced norm 1, contains $-1,i,j,k$ and an element outside $\mathbb Q\langle1,i,j,k\rangle$; a finite submonoid of a group is a group), so (iii) is in fact **[P, computer-checked]**.

**Theorem C as posed by the orchestrator** (arbitrary $D$ with positive involution, i.e. $D_{\mathbb R}\cong\prod_v M_{m_v}(\mathbb K_v)$ with some $m_v>1$,
e.g. cyclic algebras of degree $\ge3$ from Amitsur's list): not covered by §0-§5; see §7 (new in the revision, **un-refereed**). Under the hypothesis "$G$ finite with $\mathbb Q$-span $D$"
the case $D_\mathbb R\cong \mathbb K^s$ treated here is exactly: $D=\mathbb Q$ ($G=\pm1$, Klartag/Rezvyakova), $D$ a cyclotomic field ($G=\mu_m$), or $D$ a totally
definite quaternion algebra over a totally real field ($G\subset D^1$: $2T,2O,2I$, dicyclic groups). (Reason: if $D\otimes_F F_v$ is a division algebra for every
archimedean place $v$ of the centre $F$, then either $D=F$, and a field spanned by roots of unity is $\mathbb Q$ or CM, or $\dim_FD=4$ and all places are real and ramified.)

---

## 0. Setup

**0.1 Algebra [P].** $\mathcal O=\mathbb Z[G]$ is a subring of $D$ containing $1$, finitely generated as a $\mathbb Z$-module, with $\mathbb Q\mathcal O=D$:
an order in the sense of [G, Def. 1]. Fix the isomorphism $D_{\mathbb R}\cong\mathbb K^s$, $c\mapsto(c_v)_{v\le s}$. Every $u\in G$ has finite order, hence
$|u_v|=1$ for all $v$ ($|\cdot|$ is multiplicative on $\mathbb K$). $D$ is a division algebra, so $u-1\in D^\times$ for $1\ne u\in G$, hence $(u-1)_v\neq0$ for all $v$
(invertible elements of $D$ are invertible in $D_\mathbb R$).

*Hurwitz:* $D=\mathbb Q\langle1,i,j,k\rangle$, $G=\{\pm1,\pm i,\pm j,\pm k,\tfrac12(\pm1\pm i\pm j\pm k)\}$, $\mathbb Z[G]$ = Hurwitz order, $\mathbb Q G=D$ since $1,i,j,k\in G$.
*Icosians:* $K=\mathbb Q(\sqrt5)$, $\tau=\frac{1+\sqrt5}2$, $D=K\langle1,i,j,k\rangle$; $D$ is a division algebra because $x\bar x=x_0^2+\dots+x_3^2\in K\subset\mathbb R$
is nonzero for $x\ne0$. $G=2I$ is generated by $i$, $\tfrac12(1+i+j+k)$, $\tfrac12(\tau+\tau^{-1}i+j)$; it has order 120 and contains $2T\ni 1,i,j,k$.
Its $\mathbb Q$-span $A$ is a subalgebra (hence a division ring) properly containing $\mathbb Q\langle1,i,j,k\rangle$; $D$ is a left $A$-vector space, so
$\dim_{\mathbb Q}A$ divides $8$ and exceeds $4$: $A=D$. $D_{\mathbb R}=\mathbb H\times\mathbb H$ via the two embeddings $\sqrt5\mapsto\pm\sqrt5$ applied to coordinates.

**0.2 The Euclidean space and the action [P].** $E:=D_{\mathbb R}^k=\bigoplus_{v\le s}\mathbb K^k$, elements $x=(x^{(v)})_v$, $x^{(v)}\in\mathbb K^k$ **column** vectors;
$n=\dim_{\mathbb R}E=k\delta$. Inner product $\langle x,y\rangle=\sum_v\mathrm{Re}(x^{(v)*}y^{(v)})$, $|x|^2=\sum_v\sum_{i\le k}|x^{(v)}_i|^2$ (the standard
Euclidean structure of $\mathbb R^n$; $x^*$ = conjugate transpose). **Conventions:** matrices $M=(M_v)\in M_k(D_{\mathbb R})=\prod_vM_k(\mathbb K)$ act on the **left**,
$L_M x=(M_vx^{(v)})_v$; scalars $c\in D_{\mathbb R}$ act on the **right**, $R_cx=(x^{(v)}c_v)_v$. These commute. (This is the convention of [G]: lattices $g\mathcal O^k$ are
right $\mathcal O$-modules; [V] uses the mirror convention, related by $x\mapsto\bar x$.) For $u\in G$, $R_u$ is orthogonal because $|u_v|=1$.
$G$ acts **freely** on $E\setminus\{0\}$: if $xu=x$, $x^{(v)}_i\ne0$, then $x^{(v)}_i(u_v-1)=0$ in the division ring $\mathbb K$, so $(u-1)_v=0$, so $u=1$ by 0.1.
Hence every $G$-orbit in $E\setminus\{0\}$ has exactly $w$ elements.

**0.3 Invariant forms.** Let $\mathrm{Sym}_n$ be the symmetric endomorphisms of $E$ with $\langle A,B\rangle=\mathrm{Tr}[AB]$, and
$$\mathcal Q:=\{A\in\mathrm{Sym}_n:\ R_u^{T}AR_u=A\ \ \forall u\in G\}.$$

**Lemma 0.1 [P].** $\mathcal Q=\{L_H:\ H=(H_v)_v,\ H_v\in\mathrm{Herm}_k(\mathbb K)\}$; $\dim\mathcal Q=D_k$; $\mathrm{Id}\in\mathcal Q$; and for such $H,H'$:
$\langle L_H,L_{H'}\rangle=d\sum_v\mathrm{Re\,tr}(H_vH'_v)$.

*Proof.* $R_u$ is orthogonal, so $R_u^TAR_u=A$ iff $AR_u=R_uA$. The map $c\mapsto R_c$ is $\mathbb R$-linear on $D_{\mathbb R}$ and the $\mathbb R$-span of $G$ is
$D_{\mathbb R}$ (the $\mathbb Q$-span is $D$). So $A\in\mathcal Q$ iff $A$ is symmetric and commutes with all $R_c$, $c\in D_{\mathbb R}$, i.e. is right $D_{\mathbb R}$-linear.
Commuting with the idempotents $e_v$, $A$ preserves each summand $\mathbb K^k$ and is right $\mathbb K$-linear there, hence $A(x)=A(\sum_ie_ix_i)=\sum_iA(e_i)x_i$
is left multiplication by a matrix: $A=L_M$, $M\in\prod_vM_k(\mathbb K)$. Since $(M_vx)^*=x^*M_v^*$, $\langle L_Mx,y\rangle=\langle x,L_{M^*}y\rangle$; as
$M\mapsto L_M$ is injective, $L_M$ is symmetric iff $M=M^*$. Conversely every $L_H$ commutes with every $R_u$. The dimension is $s(k+d\binom k2)$.
Finally $\mathrm{Tr}_{\mathbb R}L_M=d\sum_v\mathrm{Re\,tr}M_v$ (for $c\in\mathbb K$, left multiplication by $c$ on $\mathbb K=\mathbb R^d$ has trace $d\,\mathrm{Re}\,c$;
off-diagonal entries do not contribute). $\square$

**Lemma 0.2 (projection and the constant $q$) [P].** Let $P:\mathrm{Sym}_n\to\mathrm{Sym}_n$ be the orthogonal projection onto $\mathcal Q$. Then
(a) $P(B)=\frac1w\sum_{u\in G}R_uBR_u^T$; in particular $\|P(B)\|_{op}\le\|B\|_{op}$ and $P(xu\otimes xu)=P(x\otimes x)$ for $u \in G$;
(b) $P(x\otimes x)=L_{K}$ with $K_v=\frac1d\,x^{(v)}x^{(v)*}$;
(c) $|P(x\otimes x)|^2=\frac1d\sum_v|x^{(v)}|^4=q(x)|x|^4$, where, with $u_v:=|x^{(v)}|^2/|x|^2$ and $q_0:=\frac1{\delta}$,
$$q(x)=\frac1d\sum_vu_v^2=q_0\Big(1+s\sum_{v}\big(u_v-\tfrac1s\big)^2\Big)\in[q_0,\ sq_0].$$
*Hurwitz ($s=1$, $d=4$): $q\equiv\frac14$. Icosians ($s=2,d=4$): $q=\frac14(u_1^2+u_2^2)$, $q_0=\frac18$.*

*Proof.* (a) $B\mapsto R_uBR_u^T$ is an orthogonal map of $\mathrm{Sym}_n$ with inverse given by $u^{-1}$, so the average is self-adjoint; it is idempotent with
image $\mathcal Q$: it is $P$. The norm bound is the triangle inequality; $xu\otimes xu=R_u(x\otimes x)R_u^T$. (b) $K_v$ is Hermitian, and for $B=L_H\in\mathcal Q$,
by Lemma 0.1 and $\mathrm{Re\,tr}(XY)=\mathrm{Re\,tr}(YX)$: $\langle L_K,L_H\rangle=\sum_v\mathrm{Re\,tr}(x^{(v)}x^{(v)*}H_v)=\sum_v\mathrm{Re}(x^{(v)*}H_vx^{(v)})
=\langle L_Hx,x\rangle=\langle B,x\otimes x\rangle$. (c) $|L_K|^2=d\sum_v d^{-2}\mathrm{Re\,tr}(x^{(v)}(x^{(v)*}x^{(v)})x^{(v)*})=\frac1d\sum_v|x^{(v)}|^4$ since
$x^{(v)*}x^{(v)}=|x^{(v)}|^2$ is a real scalar; and $\sum_vu_v^2=\frac1s+\sum_v(u_v-\frac1s)^2$ as $\sum_v u_v=1$. $\square$

(`check_q.py`: $D_k=1,6,15$ for Hurwitz $k=1,2,3$; $2,12$ for icosians $k=1,2$; $q=0.25$ to 12 digits; $P(x\otimes x)=L_{xx^*/4}$ to $10^{-15}$.)

**0.4 Lattices [P].** A *$G$-lattice* is $L=\lambda\,g\,\mathcal O^k$, $\lambda>0$, $g\in GL_k(D_{\mathbb R})$ ($\mathcal O^k\subset D^k\subset E$). It is a lattice in $E$
and $LG=L$ (right multiplication commutes with $L_g$ and $\mathcal O u=\mathcal O$). If $A\in\mathcal Q$ is positive definite, $\sqrt A\in\mathcal Q$
(a polynomial in $A$), $\sqrt A=L_h$, $h\in GL_k(D_\mathbb R)$, so $\sqrt A\,L$ is again a $G$-lattice.

---

## 1. The process in $\mathcal Q$  (Klartag §2)

Let $L\subset\mathbb R^n$ be a lattice and $\mathcal Q\subseteq\mathrm{Sym}_n$ **any** linear subspace containing $\mathrm{Id}$. For $A\in\mathrm{Sym}_n$ put
$$F_A:=\{B\in\mathcal Q:\ Bx\cdot x=0\ \ \forall x\in\partial E_A\cap L\},\qquad \pi_A:=\text{orthogonal projection of }\mathrm{Sym}_n\text{ onto }F_A .\tag{13$_\mathcal Q$}$$
(Klartag's (13) is the case $\mathcal Q=\mathrm{Sym}_n$.) Let $(W_t)$ be a Dyson Brownian motion in $\mathrm{Sym}_n$, as in [K].

**Proposition 1.1 (= [K, Lemma 2.2, Prop. 2.3, Cor. 2.4] with (13) replaced by (13$_\mathcal Q$)) [V].** Let $a_0>0$ with $a_0\mathrm{Id}$ $L$-free.
There is a continuous adapted process $(A_t)_{t\ge0}$ with (A)–(E) of [K, Prop. 2.3] and $A_t=a_0\mathrm{Id}+\int_0^t\pi_s(dW_s)$, $\pi_s=\pi_{A_s}$. In addition $A_t\in\mathcal Q$ for all $t$.

*Proof.* We list what the proofs in [K] use about $F_A,\pi_A$ and check it for (13$_\mathcal Q$).
* [K, Lemma 2.1] is a deterministic statement about continuous families in $\mathrm{Sym}_n$; no change.
* [K, Lemma 2.2] uses: (i) $\pi_{M_0}$ is an orthogonal projection, not identically zero when $F_{M_0}\ne0$, **hence** there is $x_0$ with $\pi_{M_0}(x_0\otimes x_0)\ne0$
  (true since the $x\otimes x$ span $\mathrm{Sym}_n$); then $\pi_{M_0}(W_t)x_0\cdot x_0=\langle W_t,\pi_{M_0}(x_0\otimes x_0)\rangle$ is a nondegenerate one-dimensional Brownian
  motion, giving (16); (ii) $Bx\cdot x=0$ for $B\in F_{M_0}$, $x\in\partial E_{M_0}\cap L$ (true by definition), giving (11). Nothing else.
* [K, Prop. 2.3] uses in addition: $F_A$ depends on $A$ only through the contact set $\partial E_A\cap L$ (true), the bounds (7),(8) for lattices (unchanged),
  and the strong Markov property of $W$. The recursion stops when $F_{\tau_{i-1}}=\{0\}$. Each step adds at least one contact point (Lemma 2.1(C)), so $M\le\tilde C_L$.
  (In fact $M\le\dim\mathcal Q$: if $x$ is a new contact point at time $\tau_i$, then $A_tx\cdot x>1$ for $\tau_{i-1}\le t<\tau_i$ and $A_{\tau_i}x\cdot x=1$, so $\langle A_{\tau_i}-A_{\tau_{i-1}},x\otimes x\rangle\ne0$ with $A_{\tau_i}-A_{\tau_{i-1}}\in F_{\tau_{i-1}}$; hence $F_{\tau_i}\subsetneq F_{\tau_{i-1}}$, the dimension drops strictly at each step. For $G$-lattices the whole orbit $xG$ becomes contact at the same time, Lemma 1.3.)
* [K, Cor. 2.4] uses only (A),(B) and boundedness of $M$.
Finally $A_t-a_0\mathrm{Id}$ is a sum of elements $\pi_{\tau_{i-1}}(\cdot)\in F_{\tau_{i-1}}\subseteq\mathcal Q$ and $\mathrm{Id}\in\mathcal Q$. $\square$

**Remark 1.2.** $\pi_s=\pi_sP$ and $\pi_s(dW_s)=\pi_s(d(PW)_s)$ where $PW$ is a standard Brownian motion in $\mathcal Q$: this is "Klartag's process run inside $\mathcal Q$".
Writing it with the ambient Dyson motion makes all of [K, §2–3] literally applicable.

From now on $\mathcal Q$ is the space of 0.3 and $L$ is a $G$-lattice.

**Lemma 1.3 (orbit counting) [P].** Almost surely for all $t$: $\partial E_t\cap L$ is $G$-invariant and
$$N_t:=\dim F_t=\mathrm{Tr}[\pi_t]\ \ge\ D_k-\frac{|\partial E_t\cap L|}{w}.$$
*Proof.* $A_t\in\mathcal Q$ gives $A_t(xu)\cdot(xu)=A_tx\cdot x$, and $Lu=L$. $F_t$ is the orthogonal complement in $\mathcal Q$ of $\mathrm{span}\{P(x\otimes x):x\in\partial E_t\cap L\}$
(as $Bx\cdot x=\langle B,P(x\otimes x)\rangle$ for $B\in\mathcal Q$). By Lemma 0.2(a) $P(x\otimes x)$ is constant on $G$-orbits, and by 0.2 all orbits in $L\setminus\{0\}$
have size $w$; so the span has dimension $\le|\partial E_t\cap L|/w$. $\square$

## 2. Volume  (Klartag §3)

**Lemma 2.1 (= [K, Lemma 3.1], Maurey) [V].** $A_t=a_0\mathrm{Id}+\frac12(W^{(1)}_t+W^{(2)}_t)$ with two Dyson Brownian motions in $\mathrm{Sym}_n$.
*Proof (written out).* $\pi_s$ is an adapted (piecewise constant, by Prop. 1.1) orthogonal projection of $\mathrm{Sym}_n$, so $U_s:=2\pi_s-\mathrm{Id}$ is an adapted orthogonal transformation of $\mathrm{Sym}_n$. Put $W^{(1)}=W$, $W^{(2)}_t=\int_0^tU_s(dW_s)$. $W^{(2)}$ is a continuous martingale in $\mathrm{Sym}_n$ with $d[\langle W^{(2)},B\rangle,\langle W^{(2)},B'\rangle]_t=\langle U_s^TB,U_s^TB'\rangle dt=\langle B,B'\rangle dt$, so it is a Dyson Brownian motion by Lévy's characterisation, and $\frac12(W^{(1)}_t+W^{(2)}_t)=\int_0^t\pi_s(dW_s)=A_t-a_0\mathrm{Id}$. Nothing about the range of $\pi_s$ is used, so [K, Cor. 3.2] holds with Klartag's same constants. $\square$

**Corollary 2.2 (= [K, Cor. 3.2]) [V].** For $t>0$, $r\ge\sqrt{tn}$: $\ \mathbb P(\|A_t-a_0\mathrm{Id}\|_{op}\ge C_0r)\le Ce^{-r^2/t}$, with Klartag's universal $C_0,C$.

**Lemma 2.3 (= [K, Lemma 3.3]) [V].** $\mathbb E\log\det A_T\le n\log a_0-\frac12\int_0^T\mathbb E[\|A_t\|_{op}^{-2}N_t]\,dt$, $N_t=\mathrm{Tr}[\pi_t]=\dim F_t$.
*Proof (the points where $\pi_t$ enters, written out).* Itô's formula for $\log\det$ ($\nabla\log\det A=A^{-1}$, $\nabla^2\log\det A(B,B)=-\mathrm{Tr}[A^{-1}BA^{-1}B]$) and $dA_t=\pi_t(dW_t)$ give $d\log\det A_t=\langle\pi_t(A_t^{-1}),dW_t\rangle-\frac12\delta_tdt$, $\delta_t=\sum_\alpha\mathrm{Tr}[A_t^{-1}E_\alpha A_t^{-1}E_\alpha]$ for an orthonormal basis $(E_\alpha)$ of $F_t$. In an orthonormal eigenbasis of $A_t$ (eigenvalues $\lambda_i>0$): $\mathrm{Tr}[A_t^{-1}EA_t^{-1}E]=\sum_{i,j}E_{ij}^2/(\lambda_i\lambda_j)\ge\|A_t\|_{op}^{-2}|E|^2$, so $\delta_t\ge\|A_t\|_{op}^{-2}\dim F_t$ -- valid for an arbitrary orthogonal projection $\pi_t$. The stochastic integral has zero expectation because $|\pi_t(A_t^{-1})|\le|A_t^{-1}|$ and $\mathbb E\int_0^T|A_t^{-1}|^2dt<\infty$ exactly as in [K] ((32) and Cor. 3.2, which do not involve $F_t$). $\square$

**Proposition 2.4 (quantitative Prop. 3.4; cf. [R, Prop. 2]) [P].** Let $T>0$, $\rho\ge0$, $1\le a_0\le1+\frac{3\rho}n$. Then
$$\mathbb E\log\det A_T\le n\log a_0-\frac{D_kT}2+\frac1{2w}\int_0^T\mathbb E|\partial E_t\cap L|\,dt+C'D_k\big(\tfrac{\rho}{n}T+n^{1/2}T^{3/2}+e^{-n}T\big),$$
$C'$ universal. (As $D_k\le n^2$ the error is $\le C'(\rho nT+n^{5/2}T^{3/2}+n^2e^{-n}T)$, the form of [R, Prop. 2]; no upper bound on $T$ is needed.)

*Proof.* Let $S_t=\{\|A_t-a_0\mathrm{Id}\|_{op}\le C_0\sqrt{tn}\}$; $\mathbb P(S_t^c)\le Ce^{-n}$ (Cor. 2.2, $r=\sqrt{tn}$). On $S_t$,
$\|A_t\|_{op}^{-2}\ge(1+\frac{3\rho}n+C_0\sqrt{tn})^{-2}\ge1-\frac{6\rho}n-2C_0\sqrt{tn}$. Since $0\le N_t\le D_k$, pointwise $\|A_t\|_{op}^{-2}N_t\ge N_t-(\frac{6\rho}n+2C_0\sqrt{tn})D_k-\mathbf 1_{S_t^c}D_k$ (on $S_t$ by the previous line, on $S_t^c$ trivially), so
$\mathbb E[\|A_t\|_{op}^{-2}N_t]\ge\mathbb E[N_t]-(\frac{6\rho}n+2C_0\sqrt{tn})D_k-CD_ke^{-n}$. Integrate over $[0,T]$, insert into Lemma 2.3, and use Lemma 1.3. $\square$

## 3. Contact points  (Klartag §4)

**Proposition 3.1 (= [K, Prop. 4.1] with the projected variance) [P].** For fixed $T>0$ and $0\ne x\in L$,
$$\mathbb P(x\in\partial E_T)\le2\,\mathbb P\Big(Z\ge\frac{1}{\sqrt{q(x)T}}\Big(a_0-\frac1{|x|^2}\Big)\Big),\qquad Z\sim N(0,1).$$
*Proof.* As in [K]: $M_t=\langle A_t,x\otimes x\rangle-1$, $dM_t=\langle dW_t,\pi_t(x\otimes x)\rangle$. Since $F_t\subseteq\mathcal Q$, $\pi_t=\pi_tP$, so
$[M]_t=\int_0^t|\pi_sP(x\otimes x)|^2ds\le t\,|P(x\otimes x)|^2=t\,q(x)|x|^4$ by Lemma 0.2(c). The rest of Klartag's proof (Dambis–Dubins–Schwarz and the reflection principle)
goes through with $|x|^4$ replaced by $q(x)|x|^4$: $\mathbb P(x\in\partial E_T)\le2\mathbb P(B_{Tq(x)|x|^4}\ge a_0|x|^2-1)$. $\square$

Let $\Phi$ be as in [K, (48)], $\mathcal D_t=\{x:(a_0-C_0\sqrt{tn})|x|^2<1\}$, and
$$K_t(L):=\sum_{0\ne x\in L\cap\mathcal D_t}\Phi\Big(\frac{1}{\sqrt{q(x)t}}\Big(a_0-\frac1{|x|^2}\Big)\Big).$$

**Proposition 3.2 (= [K, Prop. 4.2]) [V].** $\mathbb E|\partial E_t\cap L|\le2K_t(L)+Ce^{-cn}$. (Proof: verbatim [K, proof of Prop. 4.2]; it uses Prop. 3.1 for $x\in L\cap\mathcal D_t$, the bound (49) $\mathbb P(Z\ge r)\le\Phi(r)$ together with $a_0-|x|^{-2}\ge0$, Cor. 2.2 for $\mathbb P(\overline{E_t}\not\subseteq\mathcal D_t)\le Ce^{-n}$, and the bound [K, (8)] $|\partial E_t\cap L|\le2(2^n-1)$ for lattices. Recall [K, (48)]: $\Phi(r)=\min\{\frac12,\frac{1}{\sqrt{2\pi}r}e^{-r^2/2}\}$, the two bounds used in Lemma 3.3.)

**Lemma 3.3 (integral lemma; = [R, Lemma 1], reproved with uniformity) [P].** For every $\varepsilon>0$ there is $\eta\in(0,\frac14]$ such that for all $n\ge1/\eta$, $S>0$, $Y>1$ with
$$n\sqrt S\ge1/\eta,\qquad n^3S^2\le\eta,\qquad Y\sqrt S\le\eta,$$
we have, with $m=\frac{n+2}2$,
$$\int_0^S\frac{n\sqrt\sigma}2\int_0^Y\Phi(y)\,(1-y\sqrt\sigma)^{-m}\,dy\,d\sigma\ \le\ (1+\varepsilon)\frac8{n^2}e^{n^2S/8}.$$
*Proof.* Put $a=n\sqrt S/2$, $\epsilon_1=Y\sqrt S$, $\beta=1-2mS$; note $mS\le nS=(n^3S^2)^{1/2}n^{-1/2}\le\eta\le\frac14$ (using $n\ge1/\eta\ge4$), hence $\beta\ge\frac12$, $\epsilon_1\le\eta\le\frac14$ (so $y\sqrt\sigma\le\frac12$ below) and $a\ge\frac1{2\eta}$. For $0\le u\le\frac12$, $-\log(1-u)\le u+u^2$, so for
$y\le Y,\sigma\le S$: $(1-y\sqrt\sigma)^{-m}\le\exp(\tfrac n2y\sqrt\sigma+\epsilon_1+my^2S)$.
*Range $0\le y\le1$* ($\Phi\le\frac12$): contribution $\le\frac{n\sqrt S}4\,S\,e^{\epsilon_1+mS}e^{n\sqrt S/2}=\frac8{n^2}e^{n^2S/8}\cdot\frac{a^3}{4}e^{\epsilon_1+mS}e^{a-a^2/2}$, and $a^3e^{a-a^2/2}\to0$.
*Range $1\le y\le Y$* ($\Phi(y)\le e^{-y^2/2}/(\sqrt{2\pi}y)$): integrating first in $\sigma$ (substitute $\sigma=v^2$),
$\int_0^S\frac{n\sqrt\sigma}2e^{ny\sqrt\sigma/2}d\sigma=\int_0^{\sqrt S}nv^2e^{nyv/2}dv\le nS\cdot\frac2{ny}e^{ay}=\frac{2S}ye^{ay}$. Hence the contribution is at most
$$e^{\epsilon_1}\,2S\int_1^Y\frac{1}{\sqrt{2\pi}\,y^2}e^{-\beta y^2/2+ay}dy=e^{\epsilon_1}\,2S\,e^{a^2/(2\beta)}\int_1^Y\frac{1}{\sqrt{2\pi}y^2}e^{-\frac\beta2(y-a/\beta)^2}dy .$$
Here $\frac{a^2}{2\beta}=\frac{n^2S}8+\frac{a^2}{2}\frac{2mS}{\beta}=\frac{n^2S}8+O(n^3S^2)$. Split the last integral at $y_0=a(1-a^{-1/2})$: for $y\ge y_0$, $y^{-2}\le a^{-2}(1-a^{-1/2})^{-2}$ and
the Gaussian integral is $\le\beta^{-1/2}$; for $1\le y<y_0$, $y^{-2}\le1$ and $|y-a/\beta|\ge a^{1/2}$, giving $\le\beta^{-1/2}e^{-\beta a/2}=o(a^{-2})$.
So the contribution is $\le(1+o(1))\frac{2S}{a^2}e^{n^2S/8}=(1+o(1))\frac8{n^2}e^{n^2S/8}$, where every "$o(1)$" is bounded by an explicit function of $\eta$ tending to $0$ as $\eta\to0$ (the tail bound used is $\mathbb P(Z\ge u)\le\frac12e^{-u^2/2}$). The dominant relative error is $(1-a^{-1/2})^{-2}-1=O(a^{-1/2})$: in the application $a\asymp\sqrt{\log n}$, so the lemma as proved gives the factor $1+O((\log n)^{-1/4})$ (a split at $a/2$ with a second-order expansion would give $O(1/\log n)$; not needed, only $o(1)$ is used). $\square$
(`check_q.py` (4): ratio LHS/RHS $=1.045,1.017,1.013$ for $n=10^5,10^7,10^9$, $S=16\log n/n^2$, $Y=2\sqrt n$.)

**Lemma 3.4 (angular lemma) [P].** Let $u=(u_1,\dots,u_s)\sim\mathrm{Dirichlet}(a,\dots,a)$ and $b$ an integer with $0\le b<a$. Then
$$\mathbb E\exp\Big(\frac{bs}{2}\sum_v\big(u_v-\tfrac1s\big)^2\Big)\le\Big(1-\frac ba\Big)^{-bs}.$$
*Proof.* Let $x_v=su_v-1\in(-1,s-1]$, $\sum_vx_v=0$. For $x\in(-1,s-1]$: $\log(1+x)\le x-\frac{x^2}{2s}$ (for $x\le0$: $\log(1+x)\le x-\frac{x^2}2$; for $x\ge0$:
$\log(1+x)\le x-\frac{x^2}{2(1+x)}$ and $1+x\le s$). Summing, $\sum_v\log(su_v)\le-\frac1{2s}\sum_vx_v^2=-\frac s2\sum_v(u_v-\frac1s)^2$, i.e.
$\exp(\frac{bs}2\sum_v(u_v-\frac1s)^2)\le\prod_v(su_v)^{-b}$. By the Dirichlet integral, for $b<a$,
$$\mathbb E\prod_v(su_v)^{-b}=s^{-bs}\frac{\Gamma(a-b)^s\,\Gamma(sa)}{\Gamma(a)^s\,\Gamma(sa-sb)}=\frac{\prod_{i=1}^{bs}(a-\frac is)}{\prod_{j=1}^{b}(a-j)^s}\le\Big(\frac a{a-b}\Big)^{bs}.\qquad\square$$
(For $s=1$ the statement is trivial. Checked numerically in `check_q.py` (3).)

**Proposition 3.5 (shell integral; replaces [K, Lemma 4.3] and [R, Prop. 1]) [P].** Let $\rho=\rho(n)\ge1$, $\rho=o(\sqrt n)$, $\alpha=(1-\frac\rho n)^{-2}$,
$R_t=\{x:\frac1\alpha\le|x|^2<\frac1{\alpha-C_0\sqrt{tn}}\}$, and let $z=z(n)$ satisfy $\log n\le z\le 3\log n$; put $T=\frac{8z}{q_0n^2}$. Then, as $k\to\infty$ ($D,G$ fixed),
$$\frac1{\kappa_n}\int_0^T\!\!\int_{R_t}\Phi\Big(\frac{1}{\sqrt{q(x)t}}\Big(\alpha-\frac1{|x|^2}\Big)\Big)dx\,dt\ \le\ (1+o(1))\,e^{-\rho}\,\frac{8}{q_0n^2}\,e^{z}.$$
*Proof.* $q$ is $0$-homogeneous. In polar coordinates $x=r\omega$ ($\sigma$ = uniform probability on $S^{n-1}$), substituting $y=(q(\omega)t)^{-1/2}(\alpha-r^{-2})$ exactly as in
[K, proof of Lemma 4.3], the inner integral equals
$$\int_{S^{n-1}}\frac{n\sqrt{q(\omega)t}}2\,\alpha^{-\frac{n+2}2}\int_0^{C_0\sqrt{n/q(\omega)}}\Phi(y)\Big(1-\frac{y\sqrt{q(\omega)t}}{\alpha}\Big)^{-\frac{n+2}2}dy\,d\sigma(\omega).$$
Use $\alpha\ge1$, integrate in $t$ and substitute $\sigma=q(\omega)t$: by Lemma 3.3 with $S=q(\omega)T\in[q_0T,sq_0T]$, $Y=C_0\sqrt{n/q(\omega)}$
($n\sqrt S\ge\sqrt{8z}\to\infty$, $n^3S^2\le64s^2z^2/n\to0$, $Y\sqrt S=C_0\sqrt{nT}\to0$, uniformly in $\omega$) the left side is at most
$$(1+o(1))\,\alpha^{-\frac{n+2}2}\,\frac8{n^2}\int_{S^{n-1}}\frac{1}{q(\omega)}e^{n^2q(\omega)T/8}d\sigma(\omega).$$
Now $\alpha^{-(n+2)/2}=(1-\frac\rho n)^{n+2}\le e^{-\rho}$, $q\ge q_0$, and by Lemma 0.2(c) $n^2qT/8=z(1+s\sum_v(u_v-\frac1s)^2)$. For $\omega$ uniform on the sphere,
$(u_v)_v\sim\mathrm{Dirichlet}(\frac{kd}2,\dots,\frac{kd}2)$. Lemma 3.4 with $b=\lceil2z\rceil$, $a=\frac{kd}2$ gives
$\int e^{zs\sum_v(u_v-1/s)^2}d\sigma\le(1-\frac ba)^{-bs}=\exp(O(s\log^2n/k))=1+o(1)$. (For $s=1$ this step is empty: $q\equiv q_0$.) $\square$

## 4. Random $G$-lattices  (Klartag §5)

**Theorem 4.1 (Gargava [G, Theorem 16 and Remark 17]) [L].** Let $D$ be a $\mathbb Q$-division algebra, $\mathcal O\subseteq D$ an order, $k\ge2$, $\mathsf G=SL_k(D_{\mathbb R})$
(unit-norm elements of $M_k(D_\mathbb R)$), $\Gamma=SL_k(\mathcal O)$, $dg$ the $\mathsf G$-invariant probability measure on $\mathsf G/\Gamma$. Then for $f\in C_c(D_{\mathbb R}^k)$:
$\int_{\mathsf G/\Gamma}\sum_{v\in g\mathcal O^k\setminus\{0\}}f(v)\,dg=\int_{D_{\mathbb R}^k}f(x)\,dx$, $dx$ = the Lebesgue measure giving $\mathcal O^k$ covolume $1$.
**Remark 4.1a (how the citation was checked; conventions; caveats).** (i) The statement above is [G, Thm 16] verbatim: $\mathsf G=\{a\in M_k(D_\mathbb R):N(a)=1\}$, $\Gamma=\{a\in M_k(\mathcal O):N(a)=1\}$ [G, §1.2], $\mathsf G/\Gamma$ carries a unique $\mathsf G$-invariant probability measure [G, Prop. 14], "order" is [G, Def. 1] (finitely generated $\mathbb Z$-module, $\mathbb Q\otimes\mathcal O=D$, closed under multiplication, $1\in\mathcal O$ -- satisfied by $\mathbb Z[G]$, 0.1; maximality is not required). (ii) Side conventions: [G] identifies $\mathsf G/\Gamma\simeq\{g\mathcal O^k\}$ (matrices act on the left of $\mathcal O^k$), and in the proof of [G, Thm 21] the finite group $G_0\subseteq\mathcal O$ acts on $D_\mathbb R^2$ by $(v_1,v_2)\mapsto(v_1g^{-1},v_2g^{-1})$, i.e. by right multiplication: exactly the structure of 0.2/0.4. (iii) We use only the displayed identity; Lemma 4.2 derives everything else from it. (iv) **Caveat on the cited proof.** The proof of [G, Thm 16] says that $\mathsf G$ acts transitively on $D_\mathbb R^k\setminus\{0\}$. This is literally true when $D_\mathbb R$ is a division algebra (Hurwitz, $D_\mathbb R=\mathbb H$: Theorem A). For $D_\mathbb R=\mathbb K^s$, $s\ge2$ (icosians; also cyclotomic fields, the main application in [G]), the $\mathsf G$-orbit of any $v\in D^k\setminus\{0\}$ is the open dense conull set $\{x:x^{(v)}\ne0\ \forall v\}$ (all components of $v$ are nonzero since a nonzero coordinate of $v$ is invertible in $D$; and for $k\ge2$ the stabiliser of a nonzero vector in $GL_k(\mathbb K)$ contains elements of every positive norm). The argument of [G] appears to go through with $D_\mathbb R^k\setminus\{0\}$ replaced by this set, but we did NOT re-verify the proof of [G, Thm 16]; Theorem B depends on [G, Thm 16] being correct as stated for $D_\mathbb R=\mathbb H\times\mathbb H$. (v) [Unverified remark, not used] For Theorem A one could presumably use instead Vance's Hlawka-type theorem for Hurwitz lattices [V]; we did not check its exact hypotheses.

**Lemma 4.2 (Siegel formula in Klartag's normalisation) [P from 4.1].** Let $k\ge2$, $c_{\mathcal O}$ = Euclidean covolume of $\mathcal O^k$ in $E$, $\lambda=(\kappa_n/c_{\mathcal O})^{1/n}$,
and let $L=\lambda g\mathcal O^k$ with $g\Gamma$ distributed according to $dg$. Then $L$ is a $G$-lattice, $\mathrm{covol}(L)=\kappa_n$, and for every Borel $\varphi:\mathbb R^n\to[0,\infty]$
$$\mathbb E\sum_{0\ne x\in L}\varphi(x)=\frac1{\kappa_n}\int_{\mathbb R^n}\varphi .\tag{63$_G$}$$
*Proof.* Each $h\in\mathsf G$ preserves Lebesgue measure: for $0\le f\in C_c$ with $\int f>0$, Theorem 4.1 applied to $f\circ L_h$ gives $|\det h|^{-1}\int f\,dx$, while
$\sum_{v\in g\mathcal O^k\setminus0}f(hv)=\sum_{v\in hg\mathcal O^k\setminus0}f(v)$ and invariance of $dg$ give $\int f\,dx$. So $\mathrm{covol}(g\mathcal O^k)=c_{\mathcal O}$,
$\mathrm{covol}(L)=\lambda^nc_{\mathcal O}=\kappa_n$. For $\varphi\in C_c$, $\varphi\ge0$: $\mathbb E\sum_{0\ne x\in L}\varphi(x)=c_{\mathcal O}^{-1}\int\varphi(\lambda x)dx=\kappa_n^{-1}\int\varphi$.
The set function $\mu(B)=\mathbb E|B\cap L\setminus\{0\}|$ is a Borel measure (monotone convergence; $g\mapsto\sum_{\gamma\in\mathcal O^k\setminus0}\mathbf 1_B(\lambda g\gamma)$ is Borel),
finite on compacts, and agrees with $\kappa_n^{-1}\mathrm{Leb}$ on $C_c$; by uniqueness in the Riesz representation theorem $\mu=\kappa_n^{-1}\mathrm{Leb}$, and (63$_G$) follows for
all nonnegative Borel $\varphi$ by monotone approximation. $\square$

**Proposition 4.3 (good lattice; replaces [K, Prop. 5.1], [R, Cor. 1]) [P].** Let $\rho\ge1$, $\rho=o(\sqrt n)$, $a_0=(1-\frac\rho n)^{-2}$, $z$ and $T$ as in Prop. 3.5, $0<p<1-e^{-\rho}$.
For $k$ large there is a $G$-lattice $L$ of covolume $\kappa_n$ with $a_0|x|^2>1$ for all $0\ne x\in L$ **and**
$$\int_0^TK_t(L)\,dt\le\frac{1+o(1)}{p}\,e^{-\rho}\,\frac{8}{q_0n^2}\,e^z .$$
*Proof.* Take $L$ random as in Lemma 4.2. (i) By (63$_G$) with $\varphi=\mathbf 1_{\{|x|\le1-\rho/n\}}$ and Markov: $\mathbb P(\mathcal B_1)\le(1-\frac\rho n)^n\le e^{-\rho}$,
$\mathcal B_1=\{\exists\,0\ne x\in L,\ |x|\le1-\frac\rho n\}$. (Since the count is a multiple of $w$, even $\mathbb P(\mathcal B_1)\le e^{-\rho}/w$; not needed.)
(ii) Let $X(L)=\int_0^T\tilde K_t(L)dt$, $\tilde K_t$ = the sum defining $K_t$ taken over $0\ne x\in L\cap R_t$ (with $a_0=\alpha$). By Tonelli and (63$_G$),
$\mathbb EX$ equals the left side of Prop. 3.5, so $\mathbb EX\le(1+o(1))e^{-\rho}\frac8{q_0n^2}e^z$, and $\mathbb P(\mathcal B_2)\le p$ for $\mathcal B_2=\{X>\mathbb EX/p\}$.
Since $p+e^{-\rho}<1$ **both** requirements hold simultaneously with positive probability. For $L\notin\mathcal B_1\cup\mathcal B_2$: $|x|>1-\frac\rho n$ for all nonzero $x \in L$, hence
$(L\setminus0)\cap R_t=(L\setminus0)\cap\mathcal D_t$ and $\tilde K_t(L)=K_t(L)$, as in [K]. $\square$

## 5. Conclusion: proof of Theorems C′, A, B

Fix $\rho=\log\log n$, $p=1-2e^{-\rho}$, $a_0=(1-\frac\rho n)^{-2}$ (so $1\le a_0\le1+\frac{3\rho}n$ and $n\log a_0=2\rho+o(1)$), and
$$z:=\log\frac{w\,q_0n^2}{4}+\rho,\qquad T:=\frac{8z}{q_0n^2}\quad(\text{so } z=2\log n+O(\rho),\ \ T=O(\log n/n^2)\le n^{-5/3}).$$
Let $L$ be the lattice of Prop. 4.3 and run the process of Prop. 1.1. By Prop. 2.4, Prop. 3.2 and Prop. 4.3 (all error terms are $o(1)$ since $T=O(\log n/n^2)$):
$$\mathbb E\log\det A_T\le2\rho-\frac{D_kT}2+\frac1w\cdot\frac{1+o(1)}p\,e^{-\rho}\frac{8}{q_0n^2}e^z+o(1)=2\rho-\theta z+2+o(1),\qquad\theta:=\frac{4D_k}{q_0n^2}=2+\frac{4-2d}{kd},$$
using $e^z=\frac{wq_0n^2}4e^\rho$. Since $|\theta-2|\le\frac2k$ and $z=O(\log n)$: $\theta z=2z+o(1)$, hence
$$\mathbb E\log\det A_T\le-2\log\frac{w\,q_0n^2}{4}+2+o(1)=-2\log\frac{w\,q_0\,n^2}{4e}+o(1).$$
So with positive probability $\det A_T\le(1+o(1))(\frac{4e}{wq_0n^2})^2$; $A_T$ is positive definite, $L$-free and in $\mathcal Q$ (Prop. 1.1). The ellipsoid $E_{A_T}$ is $L$-free with
$$\frac{\mathrm{Vol}(E_{A_T})}{\mathrm{covol}(L)}=\det(A_T)^{-1/2}\ge(1-o(1))\frac{w\,q_0\,n^2}{4e}=(1-o(1))\frac{w\,\delta k^2}{4e},$$
and $\tilde L=\mu\sqrt{A_T}L$ ($\mu$ normalising the covolume to 1) is a $G$-lattice (0.4) with $V(\tilde L)\ge(1-o(1))\frac{w\delta k^2}{4e}$. As $D_k=\frac{\delta k^2}2(1+O(\frac1k))$ this is
$(1-o(1))\frac{wD_k}{2e}$. **Theorem C′ is proved.** Theorem A: $w=24,\delta=4,n=4k$: $\frac{24\cdot4k^2}{4e}=\frac{24k^2}e=\frac{3}{2e}n^2$. Theorem B: $w=120,\delta=8,n=8k$:
$\frac{120\cdot 8k^2}{4e}=\frac{240k^2}{e}=\frac{15}{4e}n^2$. Sanity: $G=\{\pm1\}$, $D=\mathbb Q$: $\frac{2n^2}{4e}=\frac{n^2}{2e}$, Rezvyakova's constant. $\blacksquare$

## 6. Status, caveats, gaps

1. **Proved here [P]:** 0.1–0.4, Lemma 1.3, Prop. 2.4, Prop. 3.1, Lemma 3.3, Lemma 3.4, Prop. 3.5, Lemma 4.2, Prop. 4.3, §5.
2. **[V] steps** (Prop. 1.1, Lemma 2.1, Cor. 2.2, Lemma 2.3, Prop. 3.2): Klartag's proofs are used verbatim; the only change is the definition (13$_\mathcal Q$) of $F_A$, and the
   properties of $F_A$ used in each proof are listed in Prop. 1.1. A referee should re-read [K, §2–§3, Prop. 4.2] with this substitution; I found no other dependence on $\mathcal Q=\mathrm{Sym}_n$.
3. **[L] inputs:** Gargava's Theorem 16 (the only non-classical input; hypothesis $k\ge2$; we use it with $\mathcal O=\mathbb Z[G]$, which is an order in Gargava's sense, not
   necessarily maximal). The exact definition of the norm $N$ on $M_k(D_\mathbb R)$ in [G] is irrelevant for us: Lemma 4.2 derives measure preservation from the theorem itself.
   Not independently re-proved here.
4. **Asymptotic regime:** $k\to\infty$ with $(D,G)$ fixed; "$n$ large" is used as in [K],[R]. The $o(1)$ is not effective here.
   **[Unverified remark, NOT part of the theorems]** If $(D,G)$ varies with $k$ (e.g. $D=\mathbb Q(\zeta_m)$, $w=m$ even, $\delta=\varphi(m)$, target
   $V_n\ge(1-o(1))\frac m{\varphi(m)}\frac{n^2}{4e}$, which would make the constant in the $n^2\log\log n$ bound explicit), the places where $\delta,s$ enter are:
   $T=8z\delta/n^2$ (needs $T\le n^{-5/3}$, or better: redo Prop. 2.4 with $N_t\le D_k$ instead of $N_t\le n^2$, giving the error $O(z(\rho/n+\sqrt{Tn}))$);
   Lemma 3.3 (needs $n^3(sq_0T)^2=64s^2z^2/n\to0$); Lemma 3.4 (needs $s\log^2n=o(k)$; the truth should be $\log n=o(k)$, since
   $\mathbb E\,s\sum_v(u_v-\frac1s)^2\approx\frac{2}{kd}$); and $|\theta-2|z\to0$. None of this was checked line by line.
5. **General Theorem C: formerly [GAP]; now addressed in §7 (un-refereed).** Original description of the gap: ($D_\mathbb R$ with matrix blocks $M_{m}(\mathbb K)$, $m>1$): then $x^{(v)}\in M_{km\times m}(\mathbb K)$, $P(x\otimes x)\propto(x^{(v)}x^{(v)*})_v$ and
   $q(x)\propto\sum_v\mathrm{tr}[(x^{(v)*}x^{(v)})^2]/|x|^4$ depends on the Gram matrices; one needs (a) the exact formula and its minimum $q_0$, (b) a matrix (Wishart-type) analogue of Lemma 3.4.
   Everything else above is unchanged. (Both (a) and (b) are now supplied in §7, Lemma 7.1 and Lemma 7.2 -- un-refereed.)
6. **Second-order terms:** the proof gives $V\ge(1-o(1))\,w q_0n^2/(4e)$ where the $o(1)$ contains $e^{O(\log n/k)}$ from $\theta\ne2$ and from Lemma 3.4, and $O((\log n)^{-1/4})$ (Lemma 3.3 as proved) and $1/p-1=O(1/\log n)$ terms; referee 1's direct quadrature of the left side of Prop. 3.5 divided by the claimed bound gave $2.6,\ 1.19,\ 1.045,\ 1.02$ at $n=10^3,10^4,10^5,10^6$, so the convergence is slow; these are the error terms
   from Lemma 3.3; no claim is made about the rate.
7. The initial radius is irrelevant for the constant (as in [R]); the orbit structure enters only through Lemma 1.3 (each orbit of $w$ contact points costs one dimension)
   and through $q_0=1/\delta$ (the projected variance of $A_tx\cdot x$), i.e. gain over Klartag $=\frac{w}{2}\cdot\frac{q_0}{1}=\frac{w}{2\delta}$: $3$ for Hurwitz, $7.5$ for icosians.

8. **What the two referees checked independently** (their scripts: `referee1/ref1_check.py`, `referee1/ref1_icos.py`, `referee2/ref2_check.py`): every [V] step against the text of [K] (Lemma 2.1 + Appendix, Lemma 2.2, Prop. 2.3, Cor. 2.4, Lemma 3.1, Cor. 3.2, Lemma 3.3, Prop. 4.1, Prop. 4.2, Lemma A.1); all [P] steps by hand; numerically: $q\equiv\frac14$ for Hurwitz $k=1,2,3$, $\dim\mathcal Q=1,6,15$; icosians: group of order 120, $q=\frac14(u_1^2+u_2^2)\ge\frac18$, $\dim\mathcal Q=2,12$, free action; the constants $\frac3{2e}=0.551819$, $\frac{15}{4e}=1.379548$, $\frac1{2e}=0.183940$ by numerical minimisation of the bookkeeping functional at $n=10^6$. Not checked by anybody in this project: the proof of [G, Thm 16]; novelty relative to arXiv:2606.05105 and arXiv:2610.12619.
9. **Changes in this revision:** new §7 (general Theorem C, un-refereed) and `check_blocks.py`; corrected citation in Prop. 3.2 ([K, (8)] instead of Lemma A.1); status box; explicit proofs for the [V] steps Lemma 2.1 and Lemma 2.3 and the bound $M\le\dim\mathcal Q$ in Prop. 1.1; Prop. 2.4 restated without the unused hypothesis $T\le n^{-5/3}$ and with $D_k$ in the error term; Lemma 3.3 with $\eta\le\frac14$ explicit and the honest rate $(\log n)^{-1/4}$; Remark 4.1a (exact citation, side conventions, caveat on transitivity for $s\ge2$); the Vance alternative downgraded to an unverified remark; novelty paragraph. No theorem statement had to be weakened. Theorem C beyond Theorem C′ rests on the un-refereed §7.

## 7. The general case of Theorem C (matrix blocks) -- NEW IN THE REVISION, NOT REFEREED

**Setting.** $D$ a division algebra of dimension $\delta$ over $\mathbb Q$, $G\subset D^\times$ finite with $\mathbb Q$-span $D$, $w=|G|$, $\mathcal O=\mathbb Z[G]$. $D$ is central simple over a number field $F$, so $D_\mathbb R=\prod_{v\mid\infty}D\otimes_FF_v\cong\prod_{v\le s}M_{m_v}(\mathbb K_v)$, $\mathbb K_v\in\{\mathbb R,\mathbb C,\mathbb H\}$, $d_v=\dim_\mathbb R\mathbb K_v$, $\delta=\sum_vd_vm_v^2$. The image of the finite group $G$ in $GL_{m_v}(\mathbb K_v)$ preserves a positive definite $\mathbb K_v$-Hermitian form (average the standard one), so after conjugating the $v$-th coordinate of the isomorphism we may and do assume that every $u\in G$ maps to a **unitary** matrix $u_v$ for all $v$. $E=D_\mathbb R^k=\bigoplus_vM_{km_v\times m_v}(\mathbb K_v)\ni x=(x_v)_v$ (stack the $k$ blocks vertically), $\langle x,y\rangle=\sum_v\mathrm{Re\,tr}(x_v^*y_v)$, $n=k\delta$; $M_k(D_\mathbb R)=\prod_vM_{km_v}(\mathbb K_v)$ acts on the left, $D_\mathbb R$ on the right. $R_u$ is orthogonal for $u\in G$, and $G$ acts freely on $E\setminus\{0\}$: for $1\ne u\in G$, $u-1\in D^\times$, so $R_u-\mathrm{Id}=R_{u-1}$ is invertible. Put $c_v=d_vm_v$ and
$$D_k:=\sum_v\Big(km_v+\frac{d_v\,km_v(km_v-1)}2\Big)=\frac{\delta k^2}2+O(k).$$

**Lemma 7.1 [P, un-refereed].** (a) $\mathcal Q=\{L_H:H_v\in\mathrm{Herm}_{km_v}(\mathbb K_v)\}$, $\dim\mathcal Q=D_k$, $\langle L_H,L_{H'}\rangle=\sum_vc_v\mathrm{Re\,tr}(H_vH_v')$. (b) $P(x\otimes x)=L_K$, $K_v=\frac1{c_v}x_vx_v^*$. (c) $|P(x\otimes x)|^2=\sum_v\frac1{c_v}\mathrm{tr}[(x_v^*x_v)^2]$. (d) Let $\lambda_{v,i}\ge0$ ($i\le m_v$) be the eigenvalues of the Gram matrix $x_v^*x_v\in\mathrm{Herm}_{m_v}(\mathbb K_v)$, $p_{v,i}:=c_v/\delta$ (so $\sum_{v,i}p_{v,i}=1$) and $t_{v,i}:=\delta\lambda_{v,i}/(c_v|x|^2)$ (so $\sum_{v,i}p_{v,i}t_{v,i}=1$). Then with $q_0:=1/\delta$
$$\frac{q(x)}{q_0}-1=\sum_{v,i}p_{v,i}(t_{v,i}-1)^2\ \ge0,\qquad q_0\le q\le1 .$$
*Proof.* (a) As in Lemma 0.1: $A\in\mathcal Q$ iff $A$ is symmetric and right $D_\mathbb R$-linear; $E$ is a free right $D_\mathbb R$-module of rank $k$, so $A=L_M$, $M\in M_k(D_\mathbb R)$, and the adjoint of $L_M$ is $L_{M^*}$. $L_M$ acts on each of the $m_v$ columns of $x_v$ separately, so $\mathrm{Tr}_\mathbb RL_M=\sum_vm_vd_v\mathrm{Re\,tr}M_v$. (b) $\langle L_K,L_H\rangle=\sum_v\mathrm{Re\,tr}(x_vx_v^*H_v)=\sum_v\mathrm{Re\,tr}(x_v^*H_vx_v)=\langle L_Hx,x\rangle$. (c) $|L_K|^2=\sum_vc_v\cdot c_v^{-2}\mathrm{Re\,tr}(x_vx_v^*x_vx_v^*)$ and cyclicity of $\mathrm{Re\,tr}$. (d) $\sum_{v,i}\lambda_{v,i}=|x|^2$ and $q/q_0=\delta\sum_{v,i}\lambda_{v,i}^2/(c_v|x|^4)=\sum_{v,i}p_{v,i}t_{v,i}^2$; subtract $(\sum p\,t)^2=1$. $q\le1$ since $P$ is a projection. $\square$ (Numerically: `check_blocks.py`, toy case $M_2(\mathbb R)$ with the dihedral group of order 8: $P(x\otimes x)=L_{xx^T/2}$ to $10^{-16}$, $\dim\mathcal Q=3,10,21$ for $k=1,2,3$.)

**Lemma 7.2 (angular lemma, matrix version) [P, un-refereed].** Let $\omega$ be uniform on the unit sphere of $E$, $b\ge0$, and $\kappa:=(\delta+4b\delta+2)/k<1$. Then $\ \mathbb E\exp\big(b(q(\omega)/q_0-1)\big)\le(1-\kappa)^{-2b\delta}$. In particular the left side is $1+o(1)$ if $b=O(\log n)$, $k\to\infty$, $\delta$ fixed.
*Proof.* $p_{v,i}t_{v,i}\le1$ gives $t_{v,i}\le\delta$; by the inequality $\log(1+y)\le y-\frac{y^2}{2\delta}$ on $(-1,\delta-1]$ (proof of Lemma 3.4), $\sum p_{v,i}\log t_{v,i}\le-\frac1{2\delta}(q/q_0-1)$. With $\beta_v:=2bc_v$, $\Delta_v(x):=\prod_i\lambda_{v,i}$, $B:=\sum_vm_v\beta_v=2b\delta$:
$$e^{b(q/q_0-1)}\le\prod_{v,i}t_{v,i}^{-\beta_v}=\prod_v\Big(\frac{c_v}{\delta}\Big)^{m_v\beta_v}\,|x|^{2B}\prod_v\Delta_v(x)^{-\beta_v}.$$
Let $x$ be standard Gaussian in $E$: $x=r\omega$, $r^2\sim\chi^2_n$ independent of $\omega$, and $\Delta_v(x)=r^{2m_v}\Delta_v(\omega)$, so $\mathbb E_\omega\prod_v\Delta_v(\omega)^{-\beta_v}=\prod_v\mathbb E\Delta_v(x_v)^{-\beta_v}/\mathbb E r^{-2B}$ (the $x_v$ are independent). Bartlett decomposition: by Gram-Schmidt on the $m_v$ columns of $x_v\in M_{N\times m_v}(\mathbb K_v)$, $N=km_v$, $x_v=QR$ with $Q^*Q=\mathrm{Id}$, $R$ upper triangular with diagonal $\rho_1,\dots,\rho_{m_v}>0$, $\rho_i$ = distance of column $i$ to the right $\mathbb K_v$-span of columns $1..i-1$ (real dimension $d_v(i-1)$ a.s.); conditionally on the previous columns $\rho_i^2\sim\chi^2_{\nu_{v,i}}$, $\nu_{v,i}=d_v(N-i+1)$, so the $\rho_i^2$ are independent; and $\Delta_v(x_v)=\prod_i\lambda_i(R^*R)=\prod_i\rho_i^2$ (for $\mathbb K=\mathbb H$ via the real representation $\chi:M_m(\mathbb H)\to M_{4m}(\mathbb R)$: $\det\chi(\Gamma)=\prod\lambda_i^4$ for Hermitian $\Gamma$ and $\det\chi(R^*R)=(\det\chi(R))^2=\prod\rho_i^8$). Since $\mathbb E(\chi^2_\nu)^{-\beta}=2^{-\beta}\Gamma(\frac\nu2-\beta)/\Gamma(\frac\nu2)$ for $\beta<\frac\nu2$,
$$\mathbb E_\omega\prod_v\Delta_v(\omega)^{-\beta_v}=\frac{\Gamma(\frac n2)}{\Gamma(\frac n2-B)}\prod_v\prod_{i\le m_v}\frac{\Gamma(\frac{\nu_{v,i}}2-\beta_v)}{\Gamma(\frac{\nu_{v,i}}2)}\le\Big(\frac n2\Big)^B\prod_{v,i}\Big(\frac{\nu_{v,i}}2-\beta_v-1\Big)^{-\beta_v},$$
using $\log(t-1)<\psi(t)<\log t$ ($t>1$) for $\psi=(\log\Gamma)'$, hence $(a-\beta-1)^\beta\le\Gamma(a)/\Gamma(a-\beta)\le a^\beta$ for $a-\beta>1$. Now $\frac{\nu_{v,i}}2-\beta_v-1\ge\frac{kc_v}2(1-\kappa)$ (as $d_v(i-1)\le\delta$, $2\beta_v\le4b\delta$, $kc_v\ge k$) and $n/(kc_v)=\delta/c_v$, so the right side is $\le\prod_v(\delta/c_v)^{m_v\beta_v}(1-\kappa)^{-B}$; the factors $(c_v/\delta)^{m_v\beta_v}$ cancel. $\square$
(`check_blocks.py`: Monte Carlo for $m=2,d=1$: $\mathbb E e^{z(q/q_0-1)}=1.26,1.07,1.02$ for $k=50,200,800$, $z=2\log n$; $k\cdot\mathbb E(q/q_0-1)\approx0.99$.)

**Conclusion of §7.** With Lemma 7.1 in place of Lemmas 0.1-0.2 and Lemma 7.2 (with $b=z$) in place of Lemma 3.4 in the last step of Prop. 3.5 (the hypotheses of Lemma 3.3 hold uniformly since $q_0T\le q(\omega)T\le T$), §0.4, §1-§5 go through unchanged ($\theta=4D_k/(q_0n^2)=2+O(1/k)$; [G, Thm 16] is stated for every $\mathbb Q$-division algebra), giving $V_{k\delta}\ge(1-o(1))\frac{w\,n^2}{4e\delta}=(1-o(1))\frac{wD_k}{2e}$ for every such $(D,G)$, $k\to\infty$. **Status: un-refereed.** (Remark: the gain over [R] is again $\frac{w}{2\delta}$; we did not search Amitsur's list for pairs $(D,G)$ with matrix blocks and $\frac w{2\delta}>7.5$, and make no claim about their existence.)

## Appendix: lemma-by-lemma correspondence with Klartag's paper

| [K] | here | change |
|---|---|---|
| (13) $F_A\subseteq\mathrm{Sym}_n$ | (13$_\mathcal Q$) $F_A\subseteq\mathcal Q$ | only $G$-invariant directions are allowed |
| Lemma 2.1, Lemma A.1 | used as is | none (deterministic) |
| Lemma 2.2, Prop. 2.3, Cor. 2.4 | Prop. 1.1 | none in the proofs; in addition $A_t\in\mathcal Q$ |
| Lemma 3.1, Cor. 3.2 | Lemma 2.1, Cor. 2.2 | none (same ambient Dyson motion, same $C_0$) |
| Lemma 3.3 | Lemma 2.3 | none |
| Prop. 3.4: $N_t\ge\frac{n(n+1)}2-\frac{\lvert\partial E_t\cap L\rvert}2$ | Lemma 1.3, Prop. 2.4: $N_t\ge D_k-\frac{\lvert\partial E_t\cap L\rvert}{w}$ | free orbits of size $w$; constants kept as in [R, Prop. 2] |
| Prop. 4.1: $[M]_t\le t\lvert x\rvert^4$ | Prop. 3.1: $[M]_t\le t\,q(x)\lvert x\rvert^4$ | $\pi_t=\pi_tP$, Lemma 0.2 |
| Prop. 4.2 | Prop. 3.2 | $\Phi$ evaluated at the rescaled argument |
| Lemma 4.3 / [R, Lemma 1, Prop. 1] | Lemma 3.3, Lemma 3.4, Prop. 3.5 | time rescaled by $q(\omega)$; angular average (trivial for Hurwitz) |
| (63) Siegel | Thm 4.1 [G, Thm 16], Lemma 4.2 | random $G$-lattice $\lambda g\mathcal O^k$, $k\ge2$ |
| Prop. 5.1 / [R, Cor. 1] | Prop. 4.3 | same two Markov inequalities |
| Lemma 5.2 / [R, Thm 1] | §5 | $e^z=\frac{w q_0n^2}{4}e^{\rho}$ instead of $e^{n^2T/8}=\frac{n^2}{2}e^{\rho}$ |

Bookkeeping identity behind the constant: with $z=n^2q_0T/8$,
$\ \mathbb E\log\det A_T\lesssim 2\rho-\frac{4D_k}{q_0n^2}z+\frac{8e^{-\rho}}{w q_0n^2}e^{z}$, and $\frac{4D_k}{q_0n^2}\to2$; minimum $-2\log\frac{wq_0n^2}{4e}$.
The "gain" over [R] is $\frac{wq_0}{2}=\frac{w}{2\delta}$ ($=3$ Hurwitz, $7.5$ icosians, $\frac{m}{2\varphi(m)}$ cyclotomic with $m$ even).

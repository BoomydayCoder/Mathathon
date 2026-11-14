# Zudilin's elementary proof of Apéry's theorem: a checked exposition

Source: W. Zudilin, *An elementary proof of Apéry's theorem*, arXiv:math/0202159v1 (17 Feb 2002), 8 pages.
Retrieved and read in full. Equation numbers (1)–(16) and Lemma numbers 1–7 below are Zudilin's.
Items marked **[Z]** are in the paper; items marked **[ours]** are additions (expanded calculations,
motivation, corrections, referee remarks). Scripts named `computations/*.py` are in this folder.

**Theorem (Apéry 1978).** $\zeta(3)=\sum_{m\ge1}m^{-3}$ is irrational.

Notation: $D_n=\operatorname{lcm}(1,\dots,n)$, $D_0=1$; $H_k^{(j)}=\sum_{l=1}^k l^{-j}$ ($H_0^{(j)}=0$);
$\lambda=(\sqrt2-1)^4=17-12\sqrt2=0.02943\ldots$; and for a sequence $X_n$

$$\mathcal L[X]_n:=(n+1)^3X_{n+1}-(2n+1)(17n^2+17n+5)X_n+n^3X_{n-1}\qquad(n\ge1).$$

---

## 0. The shape of the proof

### 0.1 The irrationality criterion

> **Criterion.** Let $\xi\in\mathbb R$. If there are integers $p_n,q_n$ with $0<|q_n\xi-p_n|\to0$, then $\xi\notin\mathbb Q$.

Proof: if $\xi=p/q$ with $q>0$, then $q(q_n\xi-p_n)=q_np-p_nq$ is a non-zero integer, so
$|q_n\xi-p_n|\ge1/q$, which contradicts the limit $0$. $\square$

So one needs linear forms $I_n=u_n\zeta(3)-v_n$ with three properties:

| | property | what it means |
|---|---|---|
| (A) arithmetic | $\Delta_nu_n,\ \Delta_nv_n\in\mathbb Z$ for an explicit multiplier $\Delta_n$ | controlled denominators |
| (B) size | $0<|I_n|\le C\,n^{c}\,\rho^{\,n}$ | small **and non-zero** |
| (C) balance | $\Delta_n\rho^n\to0$ | smallness beats denominators |

Here it will be $\Delta_n=D_n^3$, $\rho=\lambda$, and (C) is $D_n^{3}\lambda^n\to0$, which holds because
$D_n^{1/n}\to e$ and $e^3\lambda=0.591\ldots<1$.

### 0.2 Zudilin's architecture

The proof uses **two different series for the same number**:

* $F_n:=-\sum_{t\ge1}R_n'(t)$ (Nesterenko's series; Apéry's linear forms). Property (A) is easy for it
  (Lemma 1). Property (B) is *not visible*: its terms are not small and change sign (see §1.4).
* $\widetilde F_n:=\sum_{t\ge1}\widetilde R_n(t)$ (Ball's series). Property (B) is easy for it, because the
  terms are positive and each is small (Lemma 4). Property (A) is only available in a weaker form
  ($D_n^4$ in place of $D_n^3$; see §6.3), and $e^4\lambda=1.607\ldots>1$ fails (C).

The bridge is $F_n=\widetilde F_n$ (Lemma 7). Zudilin proves it by showing that both sequences satisfy
the same second-order recurrence $\mathcal L[X]=0$ (Lemmas 3 and 6, each from a telescoping identity,
Lemmas 2 and 5) and agree at $n=0,1$.

So in this proof **Apéry's recurrence is only a bridge between two series**. It is not used to
estimate anything, and the recurrence for $u_n$ and $v_n$ separately is never needed.

```
Lemma 1  (F_n = u_n ζ(3) − v_n,  u_n ∈ Z, D_n^3 v_n ∈ Z) ──────────────┐
Lemma 2 (telescoping for R_n) ─► Lemma 3 (L[F]=0) ──┐                   │
                                                    ├─► Lemma 7 (F=F~) ─┼─► 0 < q D_n^3 F_n → 0
Lemma 5 (telescoping for R~_n) ─► Lemma 6 (L[F~]=0)─┘                   │      contradiction
Lemma 4  (0 < F~_n < 20 (n+1)^4 λ^n) ───────────────────────────────────┤
Prime number input:  D_n ≤ C^n with C^3 λ < 1  ─────────────────────────┘
```

---

## 1. Lemma 1: the arithmetic of $F_n$

### 1.1 Definitions **[Z]**

$$R_n(t):=\left(\frac{(t-1)(t-2)\cdots(t-n)}{t(t+1)\cdots(t+n)}\right)^{2},\qquad
F_n:=-\sum_{t=1}^\infty R_n'(t).\tag{1}$$

$R_n$ has double poles at $t=0,-1,\dots,-n$, double zeros at $t=1,\dots,n$, and
$R_n(t)=t^{-2}+O(t^{-3})$ as $t\to\infty$; so $R_n'(t)=O(t^{-3})$ and the series converges absolutely.

> **Lemma 1.** $F_n=u_n\zeta(3)-v_n$ with $u_n\in\mathbb Z$ and $D_n^3v_n\in\mathbb Z$.

### 1.2 Proof, with the omitted calculations

**Step 1 (simple partial fractions).** Let $r_n(t)=\dfrac{(t-1)\cdots(t-n)}{t(t+1)\cdots(t+n)}$. The
degree of the numerator is less than that of the denominator and the poles are simple, so
$r_n(t)=\sum_{k=0}^n\frac{a_k}{t+k}$ with $a_k=\lim_{t\to-k}(t+k)r_n(t)$. At $t=-k$:

* numerator: $\prod_{j=1}^n(-k-j)=(-1)^n\frac{(n+k)!}{k!}$;
* denominator without the factor $(t+k)$: $\prod_{j=0,j\ne k}^n(j-k)=(-1)^k\,k!\,(n-k)!$.

Hence
$$a_k=(-1)^{n-k}\frac{(n+k)!}{k!^2\,(n-k)!}=(-1)^{n-k}\binom{n+k}{n}\binom nk\in\mathbb Z.$$

**Step 2 (squaring).** For $k\ne l$, $\frac1{(t+k)(t+l)}=\frac1{l-k}\left(\frac1{t+k}-\frac1{t+l}\right)$. Therefore

$$R_n(t)=\sum_{k=0}^n\left(\frac{A_{2k}}{(t+k)^2}+\frac{A_{1k}}{t+k}\right),\qquad
A_{2k}=a_k^2=\binom{n+k}n^2\binom nk^2,\qquad A_{1k}=2\sum_{l\ne k}\frac{a_ka_l}{l-k}.\tag{2}$$

Since $1\le|l-k|\le n$, we get $A_{2k}\in\mathbb Z$ and $D_nA_{1k}\in\mathbb Z$.

*Closed form* **[ours]**: $A_{1k}=\frac{d}{dt}\bigl[(t+k)^2R_n(t)\bigr]_{t=-k}
=2A_{2k}\cdot\frac{d}{dt}\log\bigl((t+k)r_n(t)\bigr)\big|_{t=-k}$, and the logarithmic derivative is
$\sum_{j=1}^n\frac1{-k-j}-\sum_{j\ne k}\frac1{j-k}=-(H_{n+k}-H_k)-(H_{n-k}-H_k)$. So
$$A_{1k}=2A_{2k}\,(2H_k-H_{n+k}-H_{n-k}),\qquad H_m:=H_m^{(1)}.$$
(Verified for $n\le12$ in `zudilin_lemma1.py`.) This is the formula behind Apéry's
$b_n=\sum_k\binom nk^2\binom{n+k}k^2c_{n,k}$.

**Step 3 (the residue sum vanishes).** $\sum_kA_{1k}=\lim_{t\to\infty}t\,R_n(t)=0$ because
$R_n(t)=O(t^{-2})$ and $\sum_k A_{2k}t/(t+k)^2\to0$.

**Step 4 (summation).** $-R_n'(t)=\sum_k\left(\frac{2A_{2k}}{(t+k)^3}+\frac{A_{1k}}{(t+k)^2}\right)$. The
sum over $k$ is finite, and each inner series converges absolutely, so we may sum over $t$ first:
$\sum_{t\ge1}(t+k)^{-j}=\zeta(j)-H_k^{(j)}$ for $j=2,3$. Thus

$$F_n=2\sum_kA_{2k}\bigl(\zeta(3)-H_k^{(3)}\bigr)+\sum_kA_{1k}\bigl(\zeta(2)-H_k^{(2)}\bigr)
=u_n\zeta(3)-v_n,$$
where the $\zeta(2)$ term has coefficient $\sum_kA_{1k}=0$ by Step 3, and

$$u_n=2\sum_{k=0}^n\binom{n+k}n^2\binom nk^2,\qquad
v_n=2\sum_{k=0}^nA_{2k}H_k^{(3)}+\sum_{k=0}^nA_{1k}H_k^{(2)}.\tag{3}$$

**Step 5 (denominators).** For $k\le n$: $D_n^3H_k^{(3)}\in\mathbb Z$ and
$D_n^3A_{1k}H_k^{(2)}=(D_nA_{1k})(D_n^2H_k^{(2)})\in\mathbb Z$. Hence $u_n\in\mathbb Z$, $D_n^3v_n\in\mathbb Z$. $\square$

**Examples [Z].** $R_0=t^{-2}$, $R_1=\frac1{t^2}+\frac4{(t+1)^2}-\frac4t+\frac4{t+1}$, so
$$F_0=2\zeta(3),\qquad F_1=10\zeta(3)-12.\tag{4}$$
($u_1=2(1+4)=10$; $v_1=2\cdot4\cdot1+4\cdot1=12$.) $u_n/2=1,5,73,1445,\dots$ are the Apéry numbers.

### 1.3 Why this construction: purpose of each feature **[ours]**

* **Why a rational function summed over integers.** A partial fraction $\frac{A}{(t+k)^j}$ summed over
  $t\ge1$ gives $A(\zeta(j)-H_k^{(j)})$: a zeta value with coefficient $A$, minus a rational number whose
  denominator divides $D_n^j$. So any $R(t)\in\mathbb Q(t)$ with poles in $\{0,-1,\dots,-n\}$ yields a
  linear form in $1,\zeta(2),\zeta(3),\dots$ with explicit coefficients. *The whole design problem is
  to choose $R$ so that unwanted zeta values drop out, denominators are small, and the sum is tiny.*
* **Why the derivative.** $\sum_tR_n(t)$ would be a form in $1$ and $\zeta(2)$ (the simple-pole terms
  telescope to a rational because $\sum A_{1k}=0$). Differentiating raises every pole order by one:
  $(t+k)^{-2}\mapsto(t+k)^{-3}$ gives $\zeta(3)$ and $(t+k)^{-1}\mapsto(t+k)^{-2}$ gives $\zeta(2)$ with
  coefficient $\sum_kA_{1k}$, which is $0$. So the *same* fact (residue sum $=0$) that makes
  $\sum R_n$ converge kills $\zeta(2)$ in $\sum R_n'$. This is how one gets a form in $1,\zeta(3)$ only.
* **Why the square.** The square makes the poles double, so that after differentiation the top pole
  order is 3. It also makes the top coefficients $A_{2k}=a_k^2$ *integers* with no denominator, so
  the only denominators come from $H_k^{(3)}$ and from $A_{1k}H_k^{(2)}$: exactly $D_n^3$.
* **Why the zeros at $t=1,\dots,n$.** They are what makes the sum small; but for $F_n$ this is hidden
  (§1.4). In the alternative proof by Nesterenko the smallness comes from a contour integral
  representation; Zudilin avoids that by passing to Ball's series.

### 1.4 Why the size of $F_n$ is not visible from (1) **[ours]**

Because $R_n$ has double zeros at $1,\dots,n$, $R_n'(t)=0$ there and $F_n=-\sum_{t>n}R_n'(t)$. But on
$t>n$ the function $R_n$ rises from $0$ to a maximum and then decays like $t^{-2}$: with $t=\tau n$ one has
$\frac1n\log R_n(\tau n)\to g(\tau)=2\bigl(2\tau\log\tau-(\tau-1)\log(\tau-1)-(\tau+1)\log(\tau+1)\bigr)$,
with $g'(\tau)=2\log\frac{\tau^2}{\tau^2-1}>0$ and $g(\tau)\to0$ as $\tau\to\infty$. So individual terms
$R_n'(t)$ are *not* exponentially small; they are positive on the rise and negative on the decay, and
$\sum_t R_n'(t)\approx\int R_n'=0$. The exponentially small value $F_n\approx\lambda^n$ is what is left
after this cancellation. No termwise estimate can see it. This is the reason the proof needs a
second representation.

---

## 2. Lemmas 2 and 3: the recurrence for $F_n$

### 2.1 Statement **[Z]**

$$s_n(t):=4(2n+1)\bigl(-2t^2+t+(2n+1)^2\bigr),\qquad S_n(t):=s_n(t)R_n(t).\tag{5}$$

> **Lemma 2.** For $n\ge1$,
> $$(n+1)^3R_{n+1}(t)-(2n+1)(17n^2+17n+5)R_n(t)+n^3R_{n-1}(t)=S_n(t+1)-S_n(t).\tag{6}$$

### 2.2 Proof

Divide by $R_n(t)$. From the definition,

$$\frac{R_{n+1}(t)}{R_n(t)}=\Bigl(\frac{t-n-1}{t+n+1}\Bigr)^2,\qquad
\frac{R_{n-1}(t)}{R_n(t)}=\Bigl(\frac{t+n}{t-n}\Bigr)^2,\qquad
\frac{R_n(t+1)}{R_n(t)}=\Bigl(\frac{t^2}{(t-n)(t+n+1)}\Bigr)^2.$$

(For the last: $r_n(t+1)=\frac{t(t-1)\cdots(t-n+1)}{(t+1)\cdots(t+n+1)}$, and dividing by $r_n(t)$ leaves
$\frac{t}{t-n}\cdot\frac{t}{t+n+1}$.) So (6) is equivalent to the identity of rational functions

$$(n+1)^3\Bigl(\frac{t-n-1}{t+n+1}\Bigr)^2-(2n+1)(17n^2+17n+5)+n^3\Bigl(\frac{t+n}{t-n}\Bigr)^2
=s_n(t+1)\Bigl(\frac{t^2}{(t-n)(t+n+1)}\Bigr)^2-s_n(t).$$

Multiplying by $(t-n)^2(t+n+1)^2$ gives a polynomial identity in $t,n$ of degree 4 in $t$ (the $t^6$
and $t^5$ terms of the right side cancel), which is checked by expansion:
`computations/zudilin_lemma2_recurrence.py`, check (a). This is "technical verification": the
*idea* is the existence of such an $s_n$, discussed in §7. $\square$

### 2.3 Lemma 3

> **Lemma 3.** $\mathcal L[F]_n=0$ for $n\ge1$, i.e.
> $(n+1)^3F_{n+1}-(2n+1)(17n^2+17n+5)F_n+n^3F_{n-1}=0.\qquad(7)$

(The paper prints (7) with $u_{n+1},u_n,u_n$; from the proof and the later use, $F_{n+1},F_n,F_{n-1}$ is
meant. See §8.)

Proof. Differentiate (6) with respect to $t$ and sum over $t=1,\dots,T$:
$$\sum_{t=1}^{T}\bigl((n+1)^3R'_{n+1}(t)-\dots+n^3R'_{n-1}(t)\bigr)=S_n'(T+1)-S_n'(1).$$
As $t\to\infty$, $S_n(t)=s_n(t)R_n(t)=-8(2n+1)+c_1t^{-1}+\dots$ (a convergent expansion in $1/t$), so
$S_n'(t)=O(t^{-2})\to0$. Letting $T\to\infty$ and using $F_m=-\sum_tR_m'(t)$:
$$\mathcal L[F]_n=S_n'(1).$$
For $n\ge1$, $R_n$ has a double zero at $t=1$ and $s_n$ is a polynomial, so $S_n$ has a double zero at
$t=1$ and $S_n'(1)=0$. $\square$

Role of the hypotheses: $n\ge1$ is needed for $R_{n-1}$ to exist and for the zero at $t=1$. The zeros
of $R_n$ at the *start of the summation range* are what make the recurrence homogeneous.

### 2.4 The recurrence for $u_n$ and $v_n$ separately (Zudilin's footnote) **[Z, expanded]**

Not needed for the proof, but it identifies $u_n$, $v_n$ as Apéry's sequences.
Multiply (6) by $(t+k)^2$ and let $t\to-k$. The coefficient of the double pole of $R_m$ at $-k$ is
$A_{2k}^{(m)}$ (zero unless $0\le k\le m$). The double pole of $S_n(t)=s_n(t)R_n(t)$ at $t=-k$ has
coefficient $s_n(-k)A^{(n)}_{2k}$, and that of $S_n(t+1)$ at $t=-k$ is $s_n(1-k)A^{(n)}_{2,k-1}$. Hence

$$(n+1)^3A^{(n+1)}_{2k}-(2n+1)(17n^2+17n+5)A^{(n)}_{2k}+n^3A^{(n-1)}_{2k}=B_{k-1}-B_k,\qquad B_k:=s_n(-k)A^{(n)}_{2k}.$$

Sum over all $k\in\mathbb Z$: the right side telescopes to $0$ (finite support). So $\mathcal L[u]_n=0$,
and then $v_n=u_n\zeta(3)-F_n$ gives $\mathcal L[v]_n=0$. With $u_0=2,u_1=10,v_0=0,v_1=12$ these are
twice Apéry's $a_n,b_n$. (Checked in `zudilin_lemma2_recurrence.py`, checks (c), (e).)

---

## 3. Ball's series and Lemma 4: the size

### 3.1 Definitions **[Z]**

$$\widetilde R_n(t):=n!^2\,(2t+n)\,\frac{(t-1)\cdots(t-n)\cdot(t+n+1)\cdots(t+2n)}{\bigl(t(t+1)\cdots(t+n)\bigr)^4},
\qquad\widetilde F_n:=\sum_{t=1}^\infty\widetilde R_n(t).\tag{8, 9}$$

$\widetilde R_n(t)=O(t^{-2n-3})$, so the series converges absolutely.

> **Lemma 4.** $0<\widetilde F_n<20(n+1)^4\lambda^{\,n}$ for all $n\ge0$, where $\lambda=(\sqrt2-1)^4$. $\quad(10)$

### 3.2 Proof

**Positivity.** $\widetilde R_n(t)=0$ for $t=1,\dots,n$ and $\widetilde R_n(t)>0$ for integers $t>n$
(all factors positive). So $\widetilde F_n=\sum_{t\ge n+1}\widetilde R_n(t)>0$.
*This is the step that supplies the "non-zero" half of property (B).*

**Case $n=0$.** $\widetilde R_0=2t^{-3}$, $\widetilde F_0=2\zeta(3)<20$. Assume $n\ge1$ from now on
(the paper's display (11) contains $(1+1/n)^{2n}$ and tacitly assumes this).

**Step 1 (an elementary Stirling substitute).** For integers $m\ge1$,
$(1+\frac1m)^m<e<(1+\frac1m)^{m+1}$, that is
$$\frac{(m+1)^m}{m^{m-1}}<e\,m<\frac{(m+1)^{m+1}}{m^m}.$$
Multiply these for $m,m+1,\dots,m+n-1$. Both outer products telescope:
$\prod_{j=m}^{m+n-1}\frac{(j+1)^j}{j^{j-1}}=\frac{(m+n)^{m+n-1}}{m^{m-1}}$ and
$\prod_{j=m}^{m+n-1}\frac{(j+1)^{j+1}}{j^{j}}=\frac{(m+n)^{m+n}}{m^{m}}$. Hence for $m,n\ge1$

$$e^{-n}\frac{(m+n)^{m+n-1}}{m^{m-1}}<m(m+1)\cdots(m+n-1)<e^{-n}\frac{(m+n)^{m+n}}{m^m}.\tag{$*$}$$

**Step 2 (apply it to each block).** For integers $t\ge n+1$, cancelling one factor $(t+n)^4$ and
the factor $(t+2n)$:
$$\widetilde R_n(t)\cdot\frac{(t+n)^5}{(2t+n)(t+2n)}
=n!^2\cdot\frac{[(t-n)\cdots(t-1)]\cdot[(t+n)\cdots(t+2n-1)]}{[t(t+1)\cdots(t+n-1)]^4}.$$
Each bracket is a product of $n$ consecutive integers. Use the upper bound in $(*)$ for
$n!$ ($m=1$: $n!<e^{-n}(n+1)^{n+1}$), for $(t-n)\cdots(t-1)$ ($m=t-n\ge1$) and for
$(t+n)\cdots(t+2n-1)$ ($m=t+n$), and the lower bound for $t\cdots(t+n-1)$ ($m=t$). The powers of $e$
are $e^{-2n}\cdot e^{-n}\cdot e^{-n}\cdot e^{+4n}=1$, and

$$\widetilde R_n(t)\cdot\frac{(t+n)^5}{(2t+n)(t+2n)}<(n+1)^{2(n+1)}\cdot
\frac{t^{5t-4}(t+2n)^{t+2n}}{(t-n)^{t-n}(t+n)^{5(t+n)-4}}.$$

Multiply by $\dfrac{t^4}{(t+n)^4(n+1)^2}$:

$$\widetilde R_n(t)\cdot\frac{t^4(t+n)}{(2t+n)(t+2n)(n+1)^2}<(n+1)^{2n}\cdot
\frac{t^{5t}(t+2n)^{t+2n}}{(t-n)^{t-n}(t+n)^{5(t+n)}}.$$

**Step 3 (scaling).** Put $t=\tau n$, $\tau>1$. The total power of $n$ in the last fraction is
$5t+(t+2n)-(t-n)-5(t+n)=-2n$, so the right side equals

$$\Bigl(1+\frac1n\Bigr)^{2n}e^{\,nf(\tau)}<e^2\,e^{\,nf(\tau)},\qquad
f(\tau):=\log\frac{\tau^{5\tau}(\tau+2)^{\tau+2}}{(\tau-1)^{\tau-1}(\tau+1)^{5(\tau+1)}}.\tag{11}$$

**Step 4 (maximising $f$ on $\tau>1$).**
$f(\tau)=5\tau\log\tau+(\tau+2)\log(\tau+2)-(\tau-1)\log(\tau-1)-5(\tau+1)\log(\tau+1)$, and the
constants from differentiating cancel ($5+1-1-5=0$):
$$f'(\tau)=\log\frac{\tau^5(\tau+2)}{(\tau-1)(\tau+1)^5}.$$
So $f'(\tau)=0$ iff $P(\tau):=\tau^5(\tau+2)-(\tau-1)(\tau+1)^5=0$. Expanding,
$P(\tau)=-2\tau^5-5\tau^4+5\tau^2+4\tau+1$, and with $x=\tau+\tfrac12$
$$P(\tau)=-x\bigl(2x^4-5x^2-\tfrac78\bigr).$$
(The substitution $x=\tau+\frac12$ is natural: $\tau\mapsto-1-\tau$, i.e. $x\mapsto-x$, is the symmetry
$t\mapsto-n-t$ of $\widetilde R_n$; see §6.2. So $P$ is odd in $x$.) The roots are $x=0$ and
$x^2=\frac54\pm\sqrt2$; the only root with $x>\frac32$ (i.e. $\tau>1$) is
$$x_0=\sqrt{\tfrac54+\sqrt2},\qquad\tau_0=-\tfrac12+\sqrt{\tfrac54+\sqrt2}=1.1322\ldots$$
(indeed $\frac54+\sqrt2>\frac94$). *It is a maximum* **[ours; the paper does not say why]**: for
$\tau>1$, $\operatorname{sign}f'(\tau)=\operatorname{sign}P(\tau)$; $P(1)=3>0$ and $P(\tau)\to-\infty$,
with a single sign change at $\tau_0$. So $f$ increases on $(1,\tau_0]$ and decreases afterwards.

**Step 5 (the value $f(\tau_0)$).** Since $f'(\tau_0)=0$,
$$f(\tau_0)=f(\tau_0)-\tau_0f'(\tau_0)=2\log(\tau_0+2)+\log(\tau_0-1)-5\log(\tau_0+1)
=\log\frac{(x_0+\frac32)^2(x_0-\frac32)}{(x_0+\frac12)^5}.$$
*Exact evaluation by hand* **[ours; the paper states the value without proof]**. Write $x=x_0$,
$x^2=\frac54+\sqrt2$.

* Numerator: $(x+\frac32)^2(x-\frac32)=(x+\frac32)(x^2-\frac94)=(x+\frac32)(\sqrt2-1)$.
* $(x+\frac12)^2=x^2+x+\frac14=x+\frac32+\sqrt2$.
* $(x+\frac12)^4=x^2+(3+2\sqrt2)x+(\frac32+\sqrt2)^2=(3+2\sqrt2)x+\frac{11}2+4\sqrt2
  =(1+\sqrt2)^2\bigl(x+\frac12+\sqrt2\bigr)$, using $(\frac{11}2+4\sqrt2)(3-2\sqrt2)=\frac12+\sqrt2$.
* $(x+\frac12)^5=(1+\sqrt2)^2\bigl(x^2+(1+\sqrt2)x+\frac14+\frac{\sqrt2}2\bigr)
  =(1+\sqrt2)^2\bigl((1+\sqrt2)x+\frac32(1+\sqrt2)\bigr)=(1+\sqrt2)^3(x+\frac32)$.

Therefore
$$e^{f(\tau_0)}=\frac{(x+\frac32)(\sqrt2-1)}{(1+\sqrt2)^3(x+\frac32)}=(\sqrt2-1)^4=\lambda.$$
So $\sup_{\tau>1}f=\log\lambda$, and for all integers $t\ge n+1$:

$$\widetilde R_n(t)\cdot\frac{t^4(t+n)}{(2t+n)(t+2n)}<e^2(n+1)^2\lambda^n.\tag{12}$$

**Step 6 (summation).** $(2t+n)(t+2n)=2t^2+5nt+2n^2$ and $t+n>t$, so
$$\frac{(2t+n)(t+2n)}{t^4(t+n)}<\frac2{t^3}+\frac{5n}{t^4}+\frac{2n^2}{t^5},$$
$$\widetilde F_n<e^2(n+1)^2\lambda^n\bigl(2\zeta(3)+5n\zeta(4)+2n^2\zeta(5)\bigr).$$
(The paper prints $\frac2{t^5}+\frac{5n}{t^4}+\frac{2n^2}{t^3}$, exponents interchanged; see §8. The
final bound is unaffected.) Now $2\zeta(3)=2.404\ldots$, $5\zeta(4)=5.411\ldots$, $2\zeta(5)=2.073\ldots$,
all at most $c(1,2,1)$ coefficientwise with $c=\frac52\zeta(4)=2.7058\ldots$, so the bracket is at most
$c(n+1)^2$, and $e^2c=19.993\ldots<20$. Hence $\widetilde F_n<20(n+1)^4\lambda^n$. $\square$

### 3.3 What matters in Lemma 4 **[ours]**

Essential: positivity, and the exponential rate $\lambda=e^{\max f}$. Everything else (the constant
20, the power $(n+1)^4$) is irrelevant to the theorem. The rate is a Laplace-method statement: the sum
of positive terms is dominated by its largest term, at $t\approx\tau_0n$. The true asymptotic is
$\widetilde F_n\approx C n^{-3/2}\lambda^n$ (expected from the recurrence; the ratios
$\widetilde F_n/\widetilde F_{n-1}=0.0086,0.0150,\dots,0.0245$ for $n=1..8$ in
`zudilin_ball_series.py` increase slowly towards $\lambda=0.0294$), consistent with
$\lambda=17-12\sqrt2$ being the small root of $x^2-34x+1$, the characteristic polynomial of
$\mathcal L$ as $n\to\infty$ (§7.3).

---

## 4. Lemmas 5 and 6: the recurrence for $\widetilde F_n$

**[Z]** Define
$$\widetilde S_n(t):=\frac{\widetilde R_n(t)}{(2t+n)(t+2n-1)(t+2n)}\,p_n(t),\tag{13}$$
$$\begin{aligned}p_n(t)={}&-t^6-(8n-1)t^5+(4n^2+27n+5)t^4+2n(67n^2+71n+15)t^3\\
&+(358n^4+339n^3+76n^2-7n-3)t^2+(384n^5+396n^4+97n^3-29n^2-17n-2)t\\
&+n(153n^5+183n^4+50n^3-30n^2-22n-4).\end{aligned}$$

> **Lemma 5.** For $n\ge1$:
> $(n+1)^3\widetilde R_{n+1}(t)-(2n+1)(17n^2+17n+5)\widetilde R_n(t)+n^3\widetilde R_{n-1}(t)=\widetilde S_n(t+1)-\widetilde S_n(t).\quad(14)$

Proof. Divide by $\widetilde R_n(t)$; all of $\widetilde R_{n\pm1}(t)/\widetilde R_n(t)$,
$\widetilde R_n(t+1)/\widetilde R_n(t)$ are explicit rational functions of $t,n$, for instance
$$\frac{\widetilde R_{n+1}(t)}{\widetilde R_n(t)}=\frac{(n+1)^2(2t+n+1)(t-n-1)(t+2n+1)(t+2n+2)}{(2t+n)(t+n+1)^5},\qquad
\frac{\widetilde R_n(t+1)}{\widetilde R_n(t)}=\frac{(2t+n+2)\,t^5\,(t+2n+1)}{(2t+n)(t-n)(t+n+1)^5},$$
and the identity becomes a polynomial identity checked by expansion
(`computations/zudilin_ball_series.py`, check (a)). $\square$

> **Lemma 6.** $\mathcal L[\widetilde F]_n=0$ for $n\ge1$.

Proof. Sum (14) over $t=1,\dots,T$: the right side is $\widetilde S_n(T+1)-\widetilde S_n(1)$.
For $n\ge1$: $\widetilde S_n(t)=O(t^{-2n-3}\cdot t^{6-3})=O(t^{-2n})\to0$; and $\widetilde S_n(1)=0$ because
$\widetilde R_n(1)=0$ while $p_n(1)$ is finite and $(2+n)(2n)(2n+1)\ne0$. So
$\mathcal L[\widetilde F]_n=-\widetilde S_n(1)=0$. $\square$

No differentiation is needed here: Ball's series is a plain sum of $\widetilde R_n$.

---

## 5. Lemma 7 and the conclusion

> **Lemma 7.** $F_n=\widetilde F_n$ for all $n\ge0$.

Proof. Both satisfy $\mathcal L[X]_n=0$ for $n\ge1$, and the leading coefficient $(n+1)^3$ never
vanishes, so each sequence is determined by its values at $n=0,1$. $\widetilde R_0=2t^{-3}$ gives
$\widetilde F_0=2\zeta(3)=F_0$. And **[Z]**
$$\widetilde R_1(t)=-\frac2{t^4}+\frac2{(t+1)^4}+\frac5{t^3}+\frac5{(t+1)^3}-\frac5{t^2}+\frac5{(t+1)^2},$$
so $\widetilde F_1=-2\zeta(4)+2(\zeta(4)-1)+5\zeta(3)+5(\zeta(3)-1)-5\zeta(2)+5(\zeta(2)-1)=10\zeta(3)-12=F_1$. $\square$

> **Proof of Apéry's theorem.** Suppose $\zeta(3)=p/q$ with $p,q\in\mathbb Z_{>0}$. By Lemma 1,
> $N_n:=qD_n^3F_n=(D_n^3u_n)p-(D_n^3v_n)q\in\mathbb Z$. By Lemmas 7 and 4,
> $$0<N_n<20q(n+1)^4\,D_n^3\lambda^n.\tag{15}$$
> If $D_n\le C^n$ for large $n$ with a constant $C$ such that $C^3\lambda<1$, the right side tends to $0$:
> an integer strictly between $0$ and $1$. Contradiction. $\square$

### 5.1 The prime-number input (referee flag) **[ours]**

The paper says "using a trivial bound $D_n<3^n$"; then $27\lambda=0.7948\ldots<1$. This bound is
true but **not trivial**, and the paper gives no reference. What the proof needs is exactly

$$\limsup_{n\to\infty}D_n^{1/n}<(\sqrt2+1)^{4/3}=3.2388\ldots\qquad\text{i.e.}\quad
\psi(n):=\log D_n\le c\,n\ \text{ with }c<\tfrac43\log(1+\sqrt2)=1.1752\ldots$$

* The prime number theorem gives $\psi(n)\sim n$, i.e. $D_n=e^{n+o(n)}$: sufficient ($e^3\lambda=0.591$).
* The easy elementary bound $D_n\le 4^{\,n+o(n)}$ (from $\prod_{p\le n}p<4^n$, proved with
  $\binom{2m+1}{m}$, plus $O(\sqrt n)$ prime powers) is **not** sufficient: $64\lambda=1.88>1$.
* Chebyshev's 1852 method gives $\psi(n)<1.1056\,n$ for large $n$, and Hanson (1972) proved
  $D_n<3^n$ for all $n$ by elementary means. *(These two attributions are quoted from memory;
  the sources were not retrieved in this session.)* Either suffices.

So the elementary proof is complete modulo a Chebyshev-type estimate strictly stronger than
the "$4^n$" one. `computations/dn_bound_check.py` checks $D_n<3^n$ for $1\le n\le 10^5$ numerically
(evidence only, not a proof). We regard this as the single external input of the proof.

---

## 6. Remarks on the two series **[ours unless marked]**

### 6.1 Where each property lives

| | $F_n=-\sum R_n'(t)$ | $\widetilde F_n=\sum\widetilde R_n(t)$ |
|---|---|---|
| form in $1,\zeta(3)$ only | yes: $\zeta(2)$ dies since $\sum_k\operatorname{Res}=0$ | yes: $\zeta(2),\zeta(4)$ die by the symmetry of §6.2 |
| denominators | $u_n\in\mathbb Z$, $D_n^3v_n\in\mathbb Z$ (elementary) | elementary argument gives only $D_n\widetilde u_n$, $D_n^4\widetilde v_n\in\mathbb Z$ |
| positivity, size | hidden by cancellation | termwise, elementary |
| recurrence certificate | degree 2 in $t$ | degree 6 in $t$ over degree 3 |

### 6.2 Why Ball's series has no $\zeta(2)$, $\zeta(4)$: well-poised symmetry

Write $\widetilde R_n(t)=\sum_{j=1}^4\sum_{k=0}^n\frac{C_{jk}}{(t+k)^j}$. Then
$\widetilde F_n=\sum_{j=2}^4\bigl(\sum_kC_{jk}\bigr)\zeta(j)-\sum_{j,k}C_{jk}H_k^{(j)}$ (the $j=1$ terms
telescope since $\sum_kC_{1k}=0$). Now $\widetilde R_n(-n-t)=-\widetilde R_n(t)$: under $t\mapsto-n-t$
the factor $2t+n$ changes sign, $\prod_{j=1}^n(t-j)\leftrightarrow(-1)^n\prod_{j=n+1}^{2n}(t+j)$, and the
denominator is invariant. Comparing partial fractions, $C_{j,n-k}=(-1)^{j+1}C_{jk}$. For even $j$ this
gives $\sum_kC_{jk}=-\sum_kC_{jk}=0$. So the coefficients of $\zeta(2)$ and $\zeta(4)$ vanish, and
$\widetilde F_n=\widetilde u_n\zeta(3)-\widetilde v_n$ with $\widetilde u_n=\sum_kC_{3k}$.
The factor $(2t+n)$ ("very-well-poised") is there precisely to make the function *odd* under the
symmetry; this is the mechanism Ball and Rivoal use to isolate odd zeta values.

### 6.3 Why Ball's series alone is not enough

The elementary denominator bound for $\widetilde v_n$ is $D_n^4$ (pole order 4). Property (C) would
need $D_n^4\lambda^n\to0$, but $e^4\lambda=1.607>1$. By Lemma 7 and irrationality of $\zeta(3)$
a posteriori, $\widetilde u_n=u_n$, $\widetilde v_n=v_n$, so in truth $D_n^3$ suffices; Zudilin
records as an open question of Rivoal to prove this directly. (The script
`zudilin_ball_series.py`, check (e), confirms $\widetilde u_n=u_n$, $\widetilde v_n=v_n$ exactly for
$n\le8$ by partial fractions, independently of the recurrence.)

### 6.4 Lemma 7 without the recurrence **[Z]**

Zudilin notes that $F_n=\widetilde F_n$ is the case $a=3n+2$, $b=c=d=e=f=n+1$ of Bailey's
transformation (16) between a very-well-poised ${}_7F_6(1)$ and a Barnes-type integral, the latter
being Nesterenko's integral for $F_n$. So the recurrence can be replaced by a classical
hypergeometric identity. We have not verified (16) in this exposition.

---

## 7. Where the recurrence comes from **[ours]**

### 7.1 Creative telescoping as linear algebra

Lemma 2 looks like magic, but finding it is a finite linear-algebra problem. One looks for
$\alpha,\beta,\gamma$ (free of $t$) and a rational $S(t)$ with
$\alpha R_{n+1}+\beta R_n+\gamma R_{n-1}=S(t+1)-S(t)$. Because all ratios $R_{n\pm1}(t)/R_n(t)$ and
$R_n(t+1)/R_n(t)$ are rational (the summand is *hypergeometric* in both $n$ and $t$), try
$S(t)=s(t)R_n(t)$. Dividing by $R_n(t)$ and clearing denominators:

$$\alpha(t-n-1)^2(t-n)^2+\beta(t-n)^2(t+n+1)^2+\gamma(t+n)^2(t+n+1)^2=s(t+1)\,t^4-s(t)(t-n)^2(t+n+1)^2.$$

The left side has degree 4. If $s$ is a polynomial of degree $d$ with leading coefficient $\sigma$,
then, since $(t-n)^2(t+n+1)^2=t^4+2t^3+\dots$ and $s(t+1)-s(t)=\sigma d\,t^{d-1}+\dots$, the right
side has degree $\le d+3$ with $t^{d+3}$-coefficient $\sigma(d-2)$. So either $d\le1$ (which turns
out to have no non-zero solution), or $d=2$, where $t^5$ cancels automatically and the coefficients
of $t^4,\dots,t^0$ give 5 equations. Unknowns:
$\alpha,\beta,\gamma$ and the 3 coefficients of $s$: 6 unknowns, 5 homogeneous linear equations, so a
non-zero solution exists. Computation (`zudilin_lemma2_recurrence.py`, check (d), $n=1,\dots,8$)
shows: the solution space is exactly one-dimensional; allowing $\deg s$ up to 6 adds nothing;
$\deg s\le1$ gives only $0$; and a first-order recurrence ($\gamma=0$) has no solution with
$\deg s\le6$. Solving gives
$\alpha=(n+1)^3$, $\beta=-(2n+1)(17n^2+17n+5)$, $\gamma=n^3$, $s=s_n$.

*So the order 2 of Apéry's recurrence is forced by a dimension count*: with only $\alpha,\beta$
(order 1) there are too few unknowns. That count is the content of Zeilberger's algorithm here.

### 7.2 Why the coefficients look as they do

* The shape $(n+1)^3,\ n^3$: put $t=n$, then $t=-n-1$, in the cleared identity above. At $t=n$: $\alpha\cdot0+\beta\cdot0+\gamma(2n)^2(2n+1)^2=s(n+1)n^4$, and at $t=-n-1$:
  $\alpha(2n+1)^2(2n+2)^2=s(-n)(n+1)^4$. So $\gamma=\frac{n^2\,s(n+1)}{4(2n+1)^2}$ and
  $\alpha=\frac{(n+1)^2s(-n)}{4(2n+1)^2}$. With $s_n(n+1)=4(2n+1)\cdot n(2n+1)$ and
  $s_n(-n)=4(2n+1)(n+1)(2n+1)$ this gives $\gamma=n^3$, $\alpha=(n+1)^3$. The cubes are
  "(square from the double zero/pole) × (one linear factor from $s_n$)".
* The reflection $n\mapsto-1-n$ fixes $R$-ratios up to swapping $R_{n+1}\leftrightarrow R_{n-1}$
  ($\frac{t-n-1}{t+n+1}\leftrightarrow\frac{t+n}{t-n}$) and fixes $(t-n)(t+n+1)$; correspondingly
  $\alpha\leftrightarrow-\gamma$, $\beta\mapsto-\beta$, $s_n\mapsto-s_n$ (by uniqueness the reflected
  solution is a multiple of the original; the multiple is $-1$). A polynomial $\beta$ that is odd
  under $n\mapsto-1-n$ is $(2n+1)\times$(a polynomial in $n(n+1)$): here $-(2n+1)(17n(n+1)+5)$.
* The number 17: as $n\to\infty$ the recurrence becomes $X_{n+1}-34X_n+X_{n-1}=0$ with roots
  $17\pm12\sqrt2=(\sqrt2\pm1)^4$. These must be the exponential growth rates of the two kinds of
  solutions: $u_n\approx(\sqrt2+1)^{4n}$ (largest term of $\sum_k\binom nk^2\binom{n+k}k^2$, at
  $k/n=1/\sqrt2$) and $F_n\approx(\sqrt2-1)^{4n}$ (Lemma 4). So $34=(\sqrt2+1)^4+(\sqrt2-1)^4$.

A fuller derivation, including how to discover $R_n$ itself, is the subject of a separate note
planned in this folder (`recurrence-origin.md`).

### 7.3 Consistency check linking Lemma 4 to the recurrence

A second-order recurrence with characteristic roots $\lambda^{-1}>\lambda$ has (Poincaré–Perron) a
one-dimensional space of solutions decaying like $\lambda^n$; all others grow like $\lambda^{-n}$.
$F_n=u_n\zeta(3)-v_n$ is the decaying one: $\zeta(3)=\lim v_n/u_n$. Apéry's original proof gets the
size of $F_n$ from the recurrence (via the Casoratian $u_nv_{n-1}-u_{n-1}v_n=\frac{c}{n^3}$);
Zudilin gets it from Ball's series instead, and uses the recurrence only for $F=\widetilde F$.

---

## 8. Referee report on the paper **[ours]**

Checked line by line. The proof is correct. Findings:

1. **Misprint in Lemma 3, eq. (7).** Printed: $(n+1)^3u_{n+1}-(2n+1)(17n^2+17n+5)u_n+n^3u_n=0$ for
   "the quantity (1)". Intended: the same with $F_{n+1},F_n,F_{n-1}$ (it is what the proof shows and
   what Lemmas 6, 7 use). The recurrence does also hold for $u_n$ (footnote; §2.4).
2. **Slip in the last display of Lemma 4.** $\frac{(2t+n)(t+2n)}{t^4(t+n)}<\frac{2}{t^3}+\frac{5n}{t^4}+\frac{2n^2}{t^5}$,
   not $\frac2{t^5}+\frac{5n}{t^4}+\frac{2n^2}{t^3}$; correspondingly $2\zeta(3)+5n\zeta(4)+2n^2\zeta(5)$.
   Both versions are $\le\frac52\zeta(4)(n+1)^2$, and $e^2\cdot\frac52\zeta(4)=19.993<20$, so (10) stands
   (with little room: the constant 20 is nearly sharp for this chain of inequalities).
3. **$n=0$ in Lemma 4.** The chain (11) uses $(1+1/n)^{2n}$; $n=0$ must be checked separately
   ($\widetilde F_0=2\zeta(3)<20$).
4. **$\tau_0$ is a maximum** and **$f(\tau_0)=4\log(\sqrt2-1)$** are asserted without proof; proofs in §3.2
   Steps 4–5.
5. **"Trivial bound $D_n<3^n$"** is not trivial; see §5.1. This is the only place where the proof
   relies on something not proved or cited in the paper.
6. **"Verify numerically"** in the one-line proofs of Lemmas 2 and 5 means: verify the polynomial
   identity in $t,n$ by expansion. Done symbolically in our scripts.
7. Interchange of $\sum_t$ and $\sum_k$ in Lemma 1 and termwise differentiation are legitimate
   (finite sum over $k$; absolutely convergent series).
8. $S_n(t)$ does not tend to $0$ ($S_n\to-8(2n+1)$); what the proof of Lemma 3 uses is
   $S_n'(t)=O(t^{-2})$, as the paper states.

Status of our verification: every step above has been checked by hand; the polynomial identities of
Lemmas 2 and 5, the partial fractions, the initial values, and the numerical agreement
$F_n=\widetilde F_n$ are in addition checked by the scripts listed in `computations/` (see
`NOTEBOOK.md` for the run log).

---

## 9. Dictionary with Beukers' integrals (series side)

To be completed after reading the other agent's exposition of Beukers' proof (see `GOALS.md`,
"Shared"). The series-side statement to be matched is: $F_n=u_n\zeta(3)-v_n$ with $u_n=2a_n$,
$v_n=2b_n$ (Apéry's sequences), to be compared with Beukers'
$\int_0^1\!\!\int_0^1\frac{-\log xy}{1-xy}P_n(x)P_n(y)\,dx\,dy$.

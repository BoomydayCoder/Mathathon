# The energy upper bound for $\Delta_K(\zeta(3))$

Companion to `hankel-zeta3.md` (§1, §2, §5). This note proves the bound that §5 lists as "not written".

> **Theorem.** For every $K\ge2$,
> $$\log\Delta_K(\zeta(3))\le-\tfrac92\log2\cdot K^2+R_K,$$
> $$R_K=2K\log K+3K\log(1+2\pi K)+\bigl(3+\log\tfrac6\pi\bigr)K+\tfrac{6\sqrt2}{\pi}\bigl(4+\log\tfrac98+4\log K\bigr).$$
> In particular $R_K=5K\log K+O(K)$, so $\log\Delta_K(\zeta(3))\le-E^*K^2+O(K\log K)$ with $E^*=\frac92\log2$.

The argument follows Fauzan §6 (circle regularisation and the negativity of the energy of a mass-zero
signed measure), adapted to our field. **Choice of variable.** Instead of symmetrising to the $2K$ points
$\pm x_i$ (the Vandermonde of $\pm x_i$ equals $\prod_{i<j}(x_i^2-x_j^2)^4\prod_i(2x_i)^2$, which has the wrong
power and needs a correction factor), we push everything forward to $t=x^2$. Then
$\prod_{i<j}(x_i^2-x_j^2)^2=\prod_{i<j}(t_i-t_j)^2$ is an ordinary Vandermonde squared, the kernel is
$\log|t-s|$, and no factor $\prod 2x_i$ appears at all. The equilibrium data are translated in Step 4.

## 0. Inputs (used as given)

* **(I1) Andréief.** $\Delta_K(\zeta(3))=\frac1{K!}\int_{(0,\infty)^K}\prod_{i<j}(y_i^2-y_j^2)^2\prod_i\omega_K(y_i)\,dy$,
  $\omega_K(y)=w(y)/D_K(y^2)$. (Gram matrix of $1,t,\dots,t^{K-1}$ in $t=y^2$; $\det[y_k^{2i}]^2=\prod_{i<j}(y_i^2-y_j^2)^2$.)
* **(I2)** $0<w(y)\le\frac1\pi(1+2\pi y)^3e^{-2\pi y}$ for $y>0$ (`hankel-zeta3.md` §5 (iii)).
* **(I3) Equilibrium measure** (`computations/equilibrium_s3.py`). $\rho$ is the probability measure on $[0,\frac34]$
  with density $\sigma$ restricted to $[0,\frac34]$, where
  $\sigma(x)=\frac1\pi\bigl[3\operatorname{artanh}\frac sb-\operatorname{artanh}\frac sc\bigr]$, $s=\sqrt{b^2-x^2}$, $b=\frac34$, $c=\frac54$.
  Since $\sigma$ is even with mass 2, $2\int\log|x-y|\sigma(y)dy=2\int\log|x^2-y^2|\,d\rho(y)$ for $x\ge0$, so the
  Euler–Lagrange relations read: for all $x\ge0$,
  $$2\int\log|x^2-y^2|\,d\rho(y)\le V(x)-\ell,\quad\text{with equality for }x\in[0,\tfrac34],\qquad\ell=2+6\log2. \tag{EL}$$
  Also $\int V\,d\sigma=6\log2-4$, i.e. $\int V\,d\rho=3\log2-2$.
  Here $V(x)=2\pi x+L(x)$, $L(x)=\int_0^1\log(x^2+u^2)\,du$.
* **(I4) Negativity lemma** (Saff–Totik, *Logarithmic potentials with external fields*, Lemma I.1.8; also Fauzan Lemma 6.2).
  If $\nu$ is a real signed Borel measure with compact support in $\mathbb C$, $\nu(\mathbb C)=0$, and
  $\iint|\log|z-w||\,d|\nu|(z)\,d|\nu|(w)<\infty$, then $I(\nu):=\iint\log|z-w|\,d\nu(z)\,d\nu(w)\le0$.
* **(I5) Circle means.** If $\omega$ is normalised arclength on $\{|z-a|=\varepsilon\}$, then
  $\int\log|z-u|\,d\omega(z)=\log\max(|a-u|,\varepsilon)$ for every $u\in\mathbb C$ (mean-value property of the
  harmonic function $\log|z-u|$, resp. Jensen's formula).

## 1. Scaling

**Step 1.** Put $y_i=Kx_i$. Then $\prod_{i<j}(y_i^2-y_j^2)^2=K^{2K(K-1)}\prod_{i<j}(x_i^2-x_j^2)^2$, $dy=K^Kdx$, and
$D_K(K^2x^2)=K^{2K}\prod_{j=1}^K\bigl(x^2+(j/K)^2\bigr)$. Hence
$$\Delta_K(\zeta(3))=\frac{K^{-K}}{K!}\int_{(0,\infty)^K}\prod_{i<j}(x_i^2-x_j^2)^2\prod_i\frac{w(Kx_i)}{\prod_{j\le K}(x_i^2+(j/K)^2)}\,dx,$$
since $2K(K-1)+K-2K^2=-K$.

## 2. Lower bound for $D_K$ (no error term)

**Step 2.** For fixed $x\ge0$, $u\mapsto\log(x^2+u^2)$ is nondecreasing on $[0,\infty)$, so
$\log(x^2+(j/K)^2)\ge K\int_{(j-1)/K}^{j/K}\log(x^2+u^2)\,du$. Summing over $j=1,\dots,K$:
$$\sum_{j=1}^K\log\bigl(x^2+(j/K)^2\bigr)\ge K\,L(x).$$
(For $x=0$ the integral is $-2$, finite.) No $\log K$ is lost here; the $K^{2K}$ was extracted exactly in Step 1.
Numerically the slack is $\ge0$ (min $0.005$ over a test grid).

**Step 3.** By (I2) and Step 2, for $x>0$:
$$\frac{w(Kx)}{\prod_j(x^2+(j/K)^2)}\le\frac1\pi(1+2\pi Kx)^3e^{-2\pi Kx-KL(x)}=\frac1\pi(1+2\pi Kx)^3e^{-KV(x)}.$$
Therefore, with $t_i=x_i^2$,
$$\Delta_K(\zeta(3))\le\frac{K^{-K}\pi^{-K}}{K!}\int_{(0,\infty)^K}\exp\Bigl(\sum_{i\ne j}\log|t_i-t_j|-K\sum_iV(x_i)\Bigr)\prod_i(1+2\pi Kx_i)^3\,dx. \tag{3.1}$$
(The set where two $t_i$ coincide is null and may be discarded.)

## 3. The equilibrium measure in the variable $t$

**Step 4.** Let $\mu$ be the image of $\rho$ under $x\mapsto x^2$ (a probability measure on $[0,\frac9{16}]$) and
$U^\mu(t)=\int\log|t-s|\,d\mu(s)$. Then $U^\mu(x^2)=\int\log|x^2-y^2|d\rho(y)$, so (EL) becomes
$2U^\mu(x^2)\le V(x)-\ell$ for all $x\ge0$, with equality for $x\in[0,\frac34]$, i.e. $\mu$-a.e.
Integrating the equality against $\mu$: $2I(\mu)=\int V\,d\rho-\ell$. Hence
$$I(\mu)+\ell=\tfrac12\Bigl(\int Vd\rho+\ell\Bigr)=\tfrac12(3\log2-2+2+6\log2)=\tfrac92\log2=E^*. \tag{4.1}$$
(Equivalently $E^*=\int Vd\rho-I(\mu)$; $I(\mu)=-2-\frac32\log2$.)

**Step 5 (density bound).** Using $\operatorname{artanh}z=\log(1+z)-\frac12\log(1-z^2)$ and $1-s^2/b^2=x^2/b^2$:
$\operatorname{artanh}(s/b)=\log(1+s/b)+\log(b/x)\le\log(2b/x)$. Since $\operatorname{artanh}(s/c)\ge0$,
$$0\le\sigma(x)\le f(x):=\tfrac3\pi\log\tfrac3{2x}\qquad(0<x<\tfrac34)$$
(numerically $\max(\sigma-f)=-0.22$). $f$ is positive and decreasing on $(0,\frac34]$, so for any interval $J$
of length $\le\delta\le\frac34$: $\rho(J)\le\int_0^\delta f=m(\delta):=\frac3\pi\delta\bigl(1+\log\frac3{2\delta}\bigr)$.

**Step 6 (mass of small discs).** For $a\in\mathbb C$ and $0<r\le\frac9{32}$, $B(a,r)\cap[0,\infty)$ is contained in an
interval $[\alpha,\beta]$ with $\beta-\alpha\le2r$; its preimage under $x\mapsto x^2$ is $[\sqrt\alpha,\sqrt\beta]$ with
$\sqrt\beta-\sqrt\alpha\le\sqrt{\beta-\alpha}\le\sqrt{2r}\le\frac34$. Hence $\mu(B(a,r))\le m(\sqrt{2r})$.

**Step 7 (smoothing error).** For $a\in[0,\infty)$ and $0<\varepsilon\le\frac9{32}$, by (I5) and Fubini/layer-cake
($\log(\varepsilon/d)=\int_d^\varepsilon dr/r$):
$$\int\log\max(|a-s|,\varepsilon)\,d\mu(s)-U^\mu(a)=\int_0^\varepsilon\frac{\mu(B(a,r))}{r}\,dr\le\int_0^\varepsilon\frac{m(\sqrt{2r})}r\,dr=\eta(\varepsilon),$$
$$\eta(\varepsilon):=\tfrac{3\sqrt2}\pi\sqrt\varepsilon\,\bigl(4+\log\tfrac98+\log\tfrac1\varepsilon\bigr),$$
using $\int_0^\varepsilon r^{-1/2}dr=2\sqrt\varepsilon$ and $\int_0^\varepsilon r^{-1/2}\log r\,dr=2\sqrt\varepsilon\log\varepsilon-4\sqrt\varepsilon$.
($U^\mu(a)$ is finite: see Step 8.)

## 4. The discrete energy inequality

**Step 8 (finiteness).** $\mu$ has density $g(t)=\sigma(\sqrt t)/(2\sqrt t)\le\frac3{2\pi}t^{-1/2}\log\frac3{2\sqrt t}$ on
$(0,\frac9{16}]$, so $g\in L^p$ for every $p\in(1,2)$; as $\log|\cdot|\in L^{p'}_{\rm loc}$, Hölder gives
$\sup_{|a|\le R}\int|\log|a-s||\,d\mu(s)<\infty$. For a circle measure $\omega$ of radius $\varepsilon$,
$\sup_{|z|\le R}\int|\log|z-w||\,d\omega(w)<\infty$ (the singularity is a one-dimensional log on a smooth curve).
Hence all mutual and self $|\log|$-energies of $\mu,\omega_1,\dots,\omega_K$ are finite.

**Step 9 (energy inequality).** Let $t_1,\dots,t_K\in[0,\infty)$ be distinct and $0<\varepsilon\le\frac9{32}$. Then
$$\sum_{i\ne j}\log|t_i-t_j|\le2K\sum_iU^\mu(t_i)-K^2I(\mu)-K\log\varepsilon+2K^2\eta(\varepsilon). \tag{9.1}$$
*Proof.* Let $\omega_i$ be normalised arclength on $|z-t_i|=\varepsilon$, $\tau=\frac1K\sum_i\omega_i$, $\nu=\tau-\mu$.
$\nu$ is compactly supported, has total mass $0$, and by Step 8 satisfies the hypothesis of (I4); so
$0\ge I(\nu)=I(\tau)-2\int U^\mu d\tau+I(\mu)$ (all three terms absolutely convergent). By (I5),
$I(\omega_i)=\log\varepsilon$ and, applying (I5) twice,
$\iint\log|z-w|d\omega_id\omega_j=\int\log\max(|t_i-w|,\varepsilon)d\omega_j(w)\ge\int\log|t_i-w|d\omega_j(w)=\log\max(|t_i-t_j|,\varepsilon)\ge\log|t_i-t_j|$.
So $K^2I(\tau)\ge K\log\varepsilon+\sum_{i\ne j}\log|t_i-t_j|$. By Fubini, (I5) and Step 7,
$\int U^\mu d\omega_i=\int\log\max(|t_i-s|,\varepsilon)d\mu(s)\le U^\mu(t_i)+\eta(\varepsilon)$. Combine and multiply by $K^2$. $\square$

**Step 10 (reserving decay).** Put $\varphi(x)=2U^\mu(x^2)-V(x)+\ell\le0$ ($x\ge0$). Claim: $\varphi(x)+x\le c_0:=2$.
For $0\le x\le\frac34$ this is $\varphi\le0$. For $x\ge\frac34$: $\rho\subset[0,\frac34]$ gives $0\le x^2-y^2\le x^2$,
so $2U^\mu(x^2)\le4\log x$; and $L(x)\ge2\log x$. Thus $\varphi(x)+x\le g(x):=2\log x-(2\pi-1)x+\ell$, and
$g'(x)=2/x-(2\pi-1)<0$ for $x\ge\frac34$, so $g(x)\le g(\frac34)=1.6211\ldots<2$.
(Numerically the true $\sup(\varphi+x)$ is $0.75$.) Consequently, for $K\ge1$ and all $x\ge0$,
$$2U^\mu(x^2)\le V(x)-\ell-\frac xK+\frac{c_0}K: \tag{10.1}$$
if $x\le c_0$, use $\varphi\le0$; if $x>c_0$, $\varphi+\frac xK\le c_0-x(1-\frac1K)\le c_0-c_0(1-\frac1K)=\frac{c_0}K$.

**Step 11.** Insert (10.1) into (9.1) and use (4.1): for distinct $t_i=x_i^2$,
$$\sum_{i\ne j}\log|t_i-t_j|-K\sum_iV(x_i)\le-E^*K^2-\sum_ix_i+c_0K-K\log\varepsilon+2K^2\eta(\varepsilon).\tag{11.1}$$
(The $-K^2\ell$ from $K\cdot K$ copies of $-\ell$ combines with $-K^2I(\mu)$ to $-K^2E^*$.)
Sanity check: for $K=10,30$, both at quantiles of $\rho$ and at random points of $(0,1)$, LHS $<$ RHS
(e.g. $K=30$ quantiles: $-2594\le-2298$).

## 5. Integration and the explicit error term

**Step 12.** By (3.1) and (11.1),
$$\Delta_K(\zeta(3))\le\frac{(\pi K)^{-K}}{K!}e^{-E^*K^2+c_0K-K\log\varepsilon+2K^2\eta(\varepsilon)}\Bigl(\int_0^\infty(1+2\pi Kx)^3e^{-x}dx\Bigr)^K.$$
With $a=2\pi K$: $\int_0^\infty(1+ax)^3e^{-x}dx=1+3a+6a^2+6a^3\le6(1+a)^3$.

**Step 13.** Take $\varepsilon=K^{-4}$ (allowed for $K\ge2$), so $-K\log\varepsilon=4K\log K$ and
$2K^2\eta(\varepsilon)=\frac{6\sqrt2}\pi(4+\log\frac98+4\log K)$. Use $\log K!\ge K\log K-K$. Then
$$\log\Delta_K(\zeta(3))\le-E^*K^2-(K\log K-K)-K\log(\pi K)+2K+4K\log K+K\log6+3K\log(1+2\pi K)+\tfrac{6\sqrt2}\pi\bigl(4+\log\tfrac98+4\log K\bigr),$$
which is the Theorem. **$\log K$ bookkeeping:** $-K\log K$ ($1/K!$), $-K\log K$ (Jacobian $K^{-K}$),
$+4K\log K$ (self-energy of circles, $\varepsilon=K^{-4}$), $+3K\log K$ (polynomial factor $(1+2\pi Kx)^3$ of $w$);
net $5K\log K$. Nothing of order $K^2\log K$ arises: the $K^{2K^2}$ from $D_K$ cancels exactly against the
Vandermonde scaling in Step 1. (With $\varepsilon=(K\log K)^{-2}$ the coefficient drops to 3; not needed.)

## 6. Numerical check

Bound $-E^*+R_K/K^2$ versus exact $\log\Delta_K(\zeta(3))/K^2$:

| $K$ | bound | exact | consistent |
|---|---|---|---|
| 20 | $-1.802$ | $-2.866$ | yes |
| 60 | $-2.610$ | $-3.002$ | yes |
| 100 | $-2.791$ | $-3.039$ | yes |
| $10^4$ | $-3.114$ | – | – |
| $10^6$ | $-3.119$ | – | – |

($E^*=3.11916$.) The bound is valid but weak for small $K$; the exact values approach $-E^*$ from above as expected.
Script: scratchpad `mycheck/check_b.py` (not committed).

## 7. Status

**Steps I am fully confident in:** 1 (scaling), 2–3 (Riemann lower bound for $D_K$, given (I2)), 4 (translation to $t$
and identity (4.1), given (I3)), 5–7 (density bound and smoothing error), 9 (expansion, given (I4)), 10–13.

**Steps that need a second look:**
* **(I3) is imported, not re-proved here**: the inequality (EL) off the support and $\int Vd\sigma=6\log2-4$ come from
  `equilibrium_s3.py` (symbolic virial identity; off-support inequality "proved" there and spot-checked numerically,
  $\varphi\le0$ on a grid up to $x=3$). Step 10 uses only (EL), its equality on $[0,\frac34]$, and $\operatorname{supp}\rho\subset[0,\frac34]$.
* **Step 8** is stated briefly (standard $L^p$/Hölder and circle-potential estimates); it is needed only to license (I4)
  and the absolute convergence of the expansion. A referee may want it written out.
* **(I4)** is cited (Saff–Totik Lemma I.1.8), not proved here.
* Step 5 uses the closed form of $\sigma$ from `hankel-zeta3.md` §5; if $\sigma$ were only known numerically, the
  constant $\frac3\pi$ in $f$ would need its own certificate.

# Four rational parameters share one-lcm inner arithmetic

Status: exact finite-sum theorem and arithmetic integrality proof, checked independently against polynomial moments at 24 cases. A separate mathematical audit of this extension is pending. The sextic factor's integrality proof is an **arithmetic factorial-ratio lemma**, not a combinatorial solution to problem #137. No classification of all rational parameters or historical novelty is claimed.

For \(q\in\{2,3,4,6\}\), put \(\tau=1/q\) and write
\[
 \delta_k=(1-\tau)_k(1+\tau)_k=c_q^{-k}F_{q,k}.
\]
Consecutive quotients, including the value at zero, prove the following table:

| q | c_q | F_(q,k) | H_(q,r)=F_(q,r)/((r+1)(2r+1)!) |
|---:|---:|---|---|
| 2 | 16 | (2k)!(2k+1)!/(k!)² | Catalan(r) |
| 3 | 27 | (3k+1)!/k! | Raney(3,2)(r), counting two ternary trees |
| 4 | 64 | (4k+1)!/(2k)! | Catalan(2r+1) |
| 6 | 432 | (6k+1)!k!/((3k)!(2k)!) | (6r+1)!r!/((r+1)(3r)!(2r)!(2r+1)!) |

The last sequence begins \(1,35,2002,138567\). An interpretation by a specified combinatorial class is not established here. The Catalan and Raney interpretations in the first three rows are classical.

## Exact finite-sum theorem

For the positive deformation's polynomial \(\Phi_n\), value \(A_n=\Phi_n(-(n+1/2)^2)\), and associated polynomial quotient \(\mathcal N_n\), define
\[
 W^{(q)}_{nk}=\frac{C_{nk}^2(k!)^2}{F_{q,k}},\quad
 L^{(q)}_n=\mathop{\mathrm{lcm}}_{0\le k\le n}
                       \operatorname{den}(W^{(q)}_{nk}).
\]
Then
\[
 A_n=\sum_{k=0}^n c_q^kW^{(q)}_{nk},
\]
\[
 \boxed{\frac{\mathcal N_n}{c_q}
 =\sum_{k=1}^n\sum_{m=1}^k(-1)^{m-1}c_q^{k-m}W^{(q)}_{nk}
 H_{q,m-1}\,\mathrm B(n-m+1,2m).} \tag{1}
\]
All four \(H\) factors are integers. Hence, for \(d_j=\mathrm{lcm}(1,\ldots,j)\),
\[
 \boxed{L^{(q)}_nA_n\in\mathbb Z,\qquad
 L^{(q)}_nd_{2n}\mathcal N_n/c_q\in\mathbb Z.} \tag{2}
\]
This does not identify an efficient common normalizer for the entire logarithmic rational approximant: the remaining outer \(W\) denominators and finite target tail depend on the parameter.

To prove (1), the finite Newton divided-product identity gives
\[
 \mathcal N_n=\sum_{k=1}^n C_{nk}^2\frac{(k!)^2}{\delta_k}
 \sum_{m=1}^k(-1)^{m-1}\frac{\delta_{m-1}(n-m)!}{m(n+m)!}.
\]
Insert \(\delta_j=c_q^{-j}F_{q,j}\) and multiply and divide the inner coefficient by \((2m-1)!\). This is precisely (1). The beta integral integrates the integer polynomial \(x^{n-m}(1-x)^{2m-1}\), of degree at most \(2n-1\), so one \(d_{2n}\) clears it. The first three integer-factor identifications follow by canceling factorials in the classical Catalan and Raney formulas. The last row needs the next lemma.

## The sextic factorial lemma, including exceptional primes

Put
\[
 J_r=\frac{(6r)!r!}{(3r)!(2r)!^2},\qquad
 H_{6,r}=\frac{(6r+1)J_r}{(r+1)(2r+1)}.
\]
For \(0\le x<1\),
\[
 \Lambda(x)=\lfloor6x\rfloor-\lfloor3x\rfloor-2\lfloor2x\rfloor
\]
takes values \(0,1,1,0,0,1\) on the six successive intervals of length \(1/6\). Counting multiples of prime powers in factorials gives
\(v_p(J_r)=\sum_{a\ge1}\Lambda(\{r/p^a\})\ge0\), so \(J_r\) is integral.

For an odd prime power \(d\mid2r+1\), the residue \(r\bmod d=(d-1)/2\) gives \(\Lambda(\{r/d\})=1\). Thus \(J_r\) contains the entire odd number \(2r+1\).

For a prime \(p\ge7\), each power \(d=p^a\mid r+1\) has residue \(d-1\) and contributes one to \(v_p(J_r)\). The same holds for powers of 5 of size at least 25. At the missing level 5, the extra factor \(6r+1=6(r+1)-5\) supplies at least one factor 5 whenever \(5\mid r+1\).

For the two remaining primes, let \(s_p\) denote the base-p digit sum. The factorial digit-sum formula gives
\[
 v_2(J_r)=s_2(r),\qquad v_3(J_r)=s_3(2r)/2.
\]
If \(2^a\mid r+1\), the last \(a\) binary digits of \(r\) are ones, so \(s_2(r)\ge a\). If \(3^a\mid r+1\), the last \(a\) ternary digits of \(2r\) are \(a-1\) twos followed by a one. Their sum is \(2a-1\); the total digit sum is even, because \(2r\) is even and \(3\equiv1\pmod2\), so \(s_3(2r)\ge2a\). Thus \(J_r\) also supplies the required powers of 2 and 3 from \(r+1\).

Finally \(\gcd(r+1,2r+1)=1\), so all required prime powers can be combined. This proves \((r+1)(2r+1)\mid(6r+1)J_r\), including \(r=0\), where both sides equal one. The argument is purely arithmetic and does not produce a counting interpretation for \(H_{6,r}\).

## Checks and limits

Run `python codex/computations/check_rational_endpoint_inner_factors.py`. It checks 804 exact factors (four parameters, indices 0 through 200), 5,049 residue/floor cases, and 24 independent polynomial-moment/divided-quotient reconstructions (four parameters, n=1 through 6), including one-lcm denominator clearing. No symbolic generic numerator recurrence is used in those reconstructions.

The outer normalizers are treated separately at [q=2](log2-endpoint-catalan-arithmetic.md) and [q=3](log3-ternary-raney-arithmetic.md). The present common inner theorem alone proves no new irrationality statement at q=4 or q=6, and it supplies no claim that these four parameters exhaust the possible integral inner factors.

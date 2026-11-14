# The higher-Rodrigues ζ(5) phase is irrational

Status: proved algebraic and density results for the \(m=2\) form, with
independent review of that argument by the contour-note author. The general
discriminant formulas below are proved. The listed small Galois groups are
exact computational results; no all-\(m\) Galois or contour-selection claim
is made. No historical novelty is asserted.

This strengthens the elementary infinitely-many-signs argument in
higher-rodrigues-laurent-geometry.md. That argument remains useful because
it requires neither Galois theory nor equidistribution.

## 1. A cubic field excludes rational phase

Let \(K\) be the splitting field of

\[
f(x)=2x^3-3x-2.
\]

The shifted polynomial

\[
f(x+1)=2x^3+6x^2+3x-3
\]

is Eisenstein at three. Thus \(f\) is irreducible, and its discriminant
\(-216\) is a nonsquare. It follows that
\(\operatorname{Gal}(K/\mathbb Q)=S_3\), whose unique quadratic subfield
is \(\mathbb Q(\sqrt{-6})\).

Every root of unity in \(K\) generates an abelian Galois extension of
\(\mathbb Q\). Such an extension corresponds to an abelian quotient of
\(S_3\), hence has degree at most two. The only quadratic cyclotomic
fields are \(\mathbb Q(i)\) and \(\mathbb Q(\sqrt{-3})\), neither the
quadratic subfield above. Therefore \(K\) contains no roots of unity
other than \(1,-1\).

The two nonreal critical values \(\eta,\overline\eta\) are obtained by
applying the rational function

\[
R(x)=\left(\frac{2(1+x)}{2-x}\right)^4
\]

to two roots of \(f\). In particular \(\eta/\overline\eta\in K\).
It is not one, because \(\eta\) is nonreal. It is not minus one either:
then \(\eta,\overline\eta\) would be purely imaginary and opposite, so
the second elementary symmetric sum of the three values
\(\lambda,\eta,\overline\eta\) would equal
\(\eta\overline\eta>0\). But their exact polynomial

\[
256L^3-2036001L^2-17376L-256
\]

gives that sum as \(-17376/256<0\). Consequently

\[
\boxed{\eta/\overline\eta\text{ is not a root of unity};
\qquad \arg\eta/\pi\notin\mathbb Q.}                       \tag{1}
\]

This is an algebraic proof, not an inference from decimal approximations
of an angle.

## 2. Exact sign densities

The independently audited asymptotic is

\[
\mathscr L_n=\frac{n^{-5/2}}3\Re(C\eta^n)
                 +O(n^{-7/2}|\eta|^n),\qquad C\ne0.
\]

Writing \(\theta=\arg\eta\), \(\phi=\arg C\), this becomes

\[
\frac{3n^{5/2}\mathscr L_n}{|C|\,|\eta|^n}
=\cos(n\theta+\phi)+O(1/n).
\]

By (1), irrational-rotation equidistribution applies. For every fixed
\(0<\epsilon<1\) and sufficiently large \(n\),
\[
\{\cos(n\theta+\phi)>\epsilon\}
\subseteq\{\mathscr L_n>0\}
\subseteq\{\cos(n\theta+\phi)>-\epsilon\}.
\]
The limiting frequencies of the outer sets are
\(\arccos(\epsilon)/\pi\) and
\(\arccos(-\epsilon)/\pi\), respectively. Letting \(\epsilon\) tend
to zero squeezes the upper and lower natural densities to \(1/2\).
The same argument applies to negative terms. Zero terms eventually lie
in \(\{|\cos(n\theta+\phi)|\le\epsilon\}\), whose limiting frequency
tends to zero with \(\epsilon\). Thus

\[
\boxed{\operatorname{dens}\{n:\mathscr L_n>0\}
=\operatorname{dens}\{n:\mathscr L_n<0\}=\tfrac12,\qquad
\operatorname{dens}\{n:\mathscr L_n=0\}=0.}                \tag{2}
\]

This argument uses only fixed intervals followed by
\(\epsilon\downarrow0\). It does not require an unproved shrinking-target
estimate, a Diophantine bound for the angle, or eventual nonvanishing.

## 3. A general Galois criterion for critical-value angles

The following criterion explains which part of the proof can extend.

**Lemma.** Let \(K/\mathbb Q\) be a finite Galois extension acting
transitively on roots \(x_1,\dots,x_d\), \(d\ge3\). Suppose its commutator
subgroup acts doubly transitively on those roots. Let \(R\in\mathbb Q(x)\)
be defined and nonzero at all roots, and suppose the values \(R(x_i)\)
are not all equal. Then no quotient \(R(x_i)/R(x_j)\), \(i\ne j\), is
a root of unity.

**Proof.** The commutator subgroup fixes every root of unity belonging to
\(K\), because the Galois action on a cyclotomic subfield is abelian.
If \(R(x_i)/R(x_j)=\zeta\), double transitivity therefore gives the
same equation \(R(x_a)/R(x_b)=\zeta\) for every ordered pair of distinct
indices. Taking three distinct indices first forces two of their values
to be equal and then forces \(\zeta=1\). All values would then be equal,
a contradiction. \(\square\)

In particular, if a critical value \(\alpha=R(x_i)\) is nonreal,
its conjugate is a distinct value in the same orbit, and the lemma gives
\(\arg\alpha/\pi\notin\mathbb Q\).
For \(S_d\), \(d\ge4\), the commutator \(A_d\) is doubly transitive,
so the criterion applies. It also applies to some groups that are not
full symmetric groups. For \(S_3\), the commutator is not doubly
transitive; Section 1 used the particular quadratic subfield instead.

This is a statement about algebraic phases. Densities for a new sequence
would additionally require a proved asymptotic selecting one dominant
conjugate pair. A list of critical values, or their Galois groups, does
not supply that analytic theorem.

## 4. Exact discriminants for both general critical polynomials

Write

\[
F_{m,\epsilon}(x)=(m-x)(1+x)^m-\epsilon x^{m+1},
\qquad \epsilon=\pm1.
\]

For the plus sign its degree is \(m+1\), with leading coefficient \(-2\).
For the minus sign its degree is \(m-1\), with leading coefficient
\(m(m+1)/2\) when \(m\ge2\). In either case there is the exact identity

\[
xF_{m,\epsilon}'(x)
=(m+1)F_{m,\epsilon}(x)-m(m+1)(1+x)^{m-1}.                \tag{3}
\]

No root is \(0\) or \(-1\); thus at any root \(\xi\),
\[
F_{m,\epsilon}'(\xi)
=-m(m+1)(1+\xi)^{m-1}/\xi\ne0.
\]
All these critical roots are simple.

For the plus polynomial, the products over its \(m+1\) roots are
\[
\prod\xi=(-1)^m m/2,\qquad \prod(1+\xi)=1/2.
\]
For the minus polynomial, over its \(m-1\) roots, they are
\[
\prod\xi=(-1)^{m-1}\frac2{m+1},\qquad
\prod(1+\xi)=\frac2{m(m+1)}.
\]
Insert these products and (3) into
\[
\operatorname{Disc}(F)
=(-1)^{d(d-1)/2}a^{d-2}\prod_{F(\xi)=0}F'(\xi),
\]
where \(d\) and \(a\) are the degree and leading coefficient. This gives
\[
\boxed{\operatorname{Disc}(F_{m,+})
=(-1)^{m(m+3)/2}\,2m^m(m+1)^{m+1},}                     \tag{4}
\]
\[
\boxed{\operatorname{Disc}(F_{m,-})
=(-1)^{(m-1)(m-2)/2}\,2m^{m-3}(m+1)^{m-2}
\quad(m\ge2).}                                         \tag{5}
\]
For \(m=2\), the second formula is \(1\), agreeing with the
discriminant convention for a linear polynomial.

The square class in (4) always has odd \(2\)-adic valuation; thus it
cannot give \(\mathbb Q(i)\) or \(\mathbb Q(\sqrt{-3})\). If the plus
polynomial has full symmetric Galois group, the particular
roots-of-unity argument of Section 1 therefore also applies.
But the full symmetric group is an additional hypothesis, not a
consequence of the discriminant.

## 5. Exact small-group checks and a warning against extrapolation

The standard exact Galois-group algorithm in SymPy gives:

| branch | \(m\) | degree | group order | commutator order | commutator doubly transitive |
|---|---:|---:|---:|---:|---|
| plus | 2 | 3 | 6 | 3 | no |
| plus | 3 | 4 | 24 | 12 | yes |
| plus | 4 | 5 | 120 | 60 | yes |
| plus | 5 | 6 | 720 | 360 | yes |
| minus | 3 | 2 | 2 | 1 | no |
| minus | 4 | 3 | 6 | 3 | no |
| minus | 5 | 4 | 24 | 12 | yes |
| minus | 6 | 5 | 120 | 60 | yes |
| minus | 7 | 6 | 120 | 60 | yes |

The last row is important: the degree-six minus polynomial at \(m=7\)
does **not** have Galois group \(S_6\). Its discriminant is a negative
square, so its splitting field contains \(i\). The assertion that all
these splitting fields have only the roots of unity \(1,-1\) would
therefore be false. Nevertheless the commutator criterion above applies
to this computed group.

For every row where that criterion applies, every nonreal critical
value from that branch has irrational argument relative to \(\pi\).
No new saddle selection, decay theorem, or sign-density result for the
higher-\(m\) forms is inferred from these finite algebraic computations.

Reproduce the shifted cubic, discriminants, derivative identities, and
finite Galois-group/commutator checks with
check_higher_rodrigues_phase_arithmetic.py.

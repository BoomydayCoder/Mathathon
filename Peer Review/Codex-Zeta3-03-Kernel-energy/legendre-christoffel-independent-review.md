# Independent review: finite Legendre mode deletion

Verdict: accepted. Root reviewed legendre-christoffel-cubic-tail.md in full.
Adjacent Legendre endpoint signs make Q_k a polynomial of degree k-1.
Its orthogonality follows by canceling the weight; its norm uses the leading
coefficient ratio -kappa_k/kappa_(k-1)=2(2k-1)/k and the Legendre norm.
The factor2/k is therefore exact.

The Volterra expression for G_n differs from K_0 by a polynomial of degree<n.
Fubini against a polynomial q of degree<n reduces to the Legendre pairing
with x^(-1) integral_0^x q, which has the same degree. Thus G_n is orthogonal
to that entire finite space. The difference of consecutive G's is orthogonal
to degree<n-1 and has leading coefficient kappa_n/n. Comparing with Q_n
forces the coefficient -1/n. The n=1 case uses the zero correction G_0-K_0
and remains valid. This proves the finite projection identity without using
any infinite basis or derivative ladder.

Both K_0 and G_n have integrable square logarithmic growth at zero and finite
limits at one. Hence Pythagoras and every inner product used are legitimate.
The product substitution in I_(n,0) gives <K_0,G_n>_w, which equals its squared
norm by the just-proved orthogonality. Tonelli evaluates ||K_0||²=2zeta(3).
This proves the exact tail and its strict positivity, since a polynomial
cannot cancel the logarithmic singularity. The statement is a normalization
and motivation lemma, not a complete irrationality proof by itself.

The author's exact Fraction checker covers40 norms/deletions,780 polynomial
orthogonalities and820 logarithmic orthogonalities. The root's separate
Markov checker also independently checks the same40 weighted mode norms.

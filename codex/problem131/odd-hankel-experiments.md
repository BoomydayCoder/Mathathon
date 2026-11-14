# Exact small Hankel experiments for ζ(3), ζ(5), ζ(7)

Status: finite exact rational-polynomial calculations plus numerical evaluation, not an asymptotic theorem. Independently reviewed the algebra and analytic boundary arguments of `odd-moment-generalization.md`; no error found. Hermite's formula remains its explicitly stated classical input. The positive Gram conclusion actually permits any nonnegative integer r, since D_N(y²)>0, although the stated even-r restriction is harmless.

## Experiment

For each s in {3,5,7}, K=1,...,6, N=0 and additionally N=1 when K>=2, set r=6 and h=K-N. Compute exactly

    G_s(X)[i,j]=mu_X(D_N(t)^6 t^(i+j)/D_K(t)), 0<=i,j<h.

Polynomial division is done first. Every polynomial quotient term contributes its Bernoulli moment; every simple pole contributes its full rational functional, including constant terms. This deliberately avoids the dangerous shortcut of estimating only favorable residues.

Clear the determinant's coefficient denominators and divide the resulting integer coefficients by their positive gcd. The resulting primitive integer polynomial P_s,K,N is the smallest positive rational scaling that makes this particular determinant integral. Positivity of its value follows from the Gram argument, independently of floating-point evaluation.

Exact coefficients, normalization factors, numerical log-values and sizes are in `codex/computations/odd_hankel_results.json`. Reproduce with `python codex/computations/explore_odd_hankel.py` (SymPy and mpmath). The experiment printed below uses 150 decimal digits; a separate evaluation at 300 digits agreed in every logarithm to 85 decimal digits. These are high-precision checks, not certified interval bounds.

## An exact leading-coefficient theorem

For odd s>=3, r>=1, K>N>=0 and h=K-N,

    [X^h] det G_s(X)
      =(-1)^(h(h-1)/2) product_{j=N+1}^K j^(s-1) D_N(-j²)^(r-1).

In particular its degree is exactly h, for all these s. Proof: after cancellation write the common denominator Dtail=D_K/D_N and numerator D_N^(r-1)t^(i+j). The coefficient of X in the moment matrix is V diag(c_j) V^T, with V[i,j]=(-j²)^i and c_j=j^(s-1)D_N(-j²)^(r-1)/Dtail'(-j²). The product of Dtail' over its h roots equals (-1)^(h(h-1)/2) times the squared Vandermonde determinant. Taking the determinant cancels those factors and proves the formula. The code verifies it for every experiment.

This extends the degree calculation as easily as positivity. Neither guarantees useful arithmetic normalization.

## Results

Entries below are -log(P_s,K,0(ζ(s)))/degree, rounded; larger positive values mean smaller values per degree.

| K | s=3 | s=5 | s=7 |
|---|---:|---:|---:|
| 1 | 1.5992 | -0.1378 | -0.0247 |
| 2 | 2.8236 | 1.4406 | 0.3655 |
| 3 | 4.1281 | 2.1107 | 0.0857 |
| 4 | 5.1175 | 2.2235 | 0.0183 |
| 5 | 5.9008 | 2.5104 | -1.2531 |
| 6 | 7.2980 | 3.1720 | -1.6830 |

For N=1 at K=6 the corresponding values are 5.3580, 1.3890, -4.0841. The largest coefficients have 40,49,61 decimal digits respectively. Thus a determinant can be positive and have quite structured coefficients while its primitive integer value is greater than 1. The s=7 samples demonstrate this directly.

**What this supports:** the analytic moment model and full-degree property survive changing s, but these particular small normalized candidates have substantially worse values at s=7. It is worth investigating where the normalization cost increases.

**What this does not support:** no impossibility theorem for s=7, no growth law from six K values, and no verification of Fauzan's asymptotic choice N=3K/40, h=37K/40. In fact N=0 or fixed N=1 has a different limiting parameter ratio. Even good-looking ζ(5) values do not establish irrationality without a proof that the displayed ratio tends to infinity along an infinite sequence.

## Arithmetic obstacle and next concrete test

The generic monomial moment contains B_(2e+2) times a polynomial in e of degree s-2 divided by (s-1)!. As s rises, the prime denominators from Bernoulli numbers persist while other factorial and residue factors change. The pole functional also contains j^(s-1) H_j^(s); these changes must be analyzed as a complete rational functional, not pole by pole in isolation.

The next tractable experiment is to factor the rational normalization prime by prime and compare it to a valid uniform analytic upper bound for the determinant. Then vary N/K and r together. Small K is useful for finding missing factors and testing exact local lemmas; it cannot validate the needed n²-scale margin. A faithful test of Fauzan's first proportional choice starts at K=40,N=3,h=37 and is a much larger symbolic determinant problem than the computations reported here.

# Dependency map: how the ideas lead to irrationality

All three ζ(3) proofs, and the ζ(5) proof, instantiate one scheme:

```
 (A) arithmetic:   Δ_n · (object_n) has integer coefficients
 (B) analysis:     0 < |object_n(ζ)| ≤ ε_n
 (C) balance:      Δ_n · ε_n · (cost of plugging ζ = p/q) → 0
 ───────────────────────────────────────────────────────────
           ζ = p/q  ⇒  a positive integer tends to 0
```

## 1. Zudilin's proof (exposition.md)

```
partial fractions of R_n = r_n^2 ──► Lemma 1:  F_n = u_n ζ(3) − v_n,  u_n ∈ Z, D_n^3 v_n ∈ Z     (A)
        │
        ├─ polynomial identity (⋆) [recurrence-origin.md §3] ──► Lemma 2 ──► Lemma 3: L[F] = 0 ─┐
        │                                                                                      ├─► Lemma 7: F = F~
Ball's R~_n ── degree-6 certificate ──► Lemma 5 ──► Lemma 6: L[F~] = 0 ─────────────────────────┘        │
        │                                                                                               ▼
        └─ positivity of terms + Stirling-type bound + max of f(τ) ──► Lemma 4: 0 < F~_n < 20(n+1)^4 λ^n   (B)

lcm bound:  limsup D_n^{1/n} < (√2+1)^{4/3}   [PNT, or Hanson, or the period-6 Chebyshev bound]          (C)
```

## 2. The shorter route (short-proof.md)

```
Lemma 1 ───────────────────────────────────────────────────────────────────────────► (A)
identity (⋆) ──► L[F] = 0  (differentiate, sum over t)
            └──► L[u] = 0  (double-pole coefficients, sum over k)
                      │
                      ▼
        Casoratian  u_n F_{n−1} − u_{n−1} F_n = 24/n^3
                      │        + trivial bound |F_n| ≤ 2  (so F_n/u_n → 0)
                      ▼
        F_n = u_n Σ_{m>n} 24/(m^3 u_{m−1} u_m)  > 0,   F_n < 12/(n^2 u_{n+1})
                      │        + growth u_{n+1}/u_n → 34 − 1/30 (from L[u] = 0 again)
                      ▼
                     (B) with ε_n ≈ ρ^{−n}, any ρ < 33.96
period-6 Chebyshev bound  e^{3c} = 33.50 < 33.96 ───────────────────────────────────► (C)
```

Single non-trivial computation: the degree-4 polynomial identity (⋆), solvable by hand.

## 3. Beukers' proof (via the dictionary, exposition.md §9)

```
Legendre P_n ──Mellin──► r_n(t);   ∫∫ (−log xy)/(1−xy) P_n P_n = −Σ R_n'(t) = F_n          (A) = Lemma 1
orthogonality + integration by parts ──► positive triple integral, max of kernel = λ       (B)
lcm bound                                                                                  (C)
```

## 4. Fauzan's ζ(5) proof (fauzan-status.md)

```
functional μ_X: rational moments (Bernoulli), pole values affine in X
   │
   ├─ at X = ζ(5): integration against w > 0 ──► Gram matrix ──► det > 0                    (B, non-vanishing)
   ├─ Andréief + weighted logarithmic energy ──► det ≤ exp(U K^2 + …),  U = −1.3670         (B, size)
   └─ p-adic valuation of det as a whole (3 ranges of primes) + PNT ──► cost exp(A K^2),
                                                              A_200 = 1.3496               (A)
polynomial of degree h = 37n in ζ(5); plugging p/q costs q^{37n} = exp(O(K))               (C): A + U < 0
```

## 5. What each proof needs from outside

| proof | external inputs |
|---|---|
| Zudilin as printed | "$D_n<3^n$" (unproved there) |
| short route | none beyond undergraduate analysis, if the period-6 lcm bound is included |
| Beukers | an lcm bound of the same strength |
| Fauzan | prime number theorem; computer-assisted energy estimate; (as printed) two gaps filled in the formalizations |

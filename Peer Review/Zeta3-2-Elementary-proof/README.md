# ζ(3), proof 2: a self-contained elementary proof with motivation

**Main document:** `apery-elementary.md` (proof + a motivation paragraph before every step).
**Companions:** `short-proof.md` (the same argument in compact form, with a comparison to Zudilin's
proof); `recurrence-origin.md` (how Apéry's recurrence is found by hand and why its coefficients look as
they do).

**Claim:** ζ(3) is irrational, using only partial fractions, absolutely convergent series and the
exponent of a prime in $m!$. Structure: linear forms $F_n=u_n\zeta(3)-v_n$ from one rational function (§1);
one telescoping identity read twice (§2), verified by hand coefficient by coefficient (§3); Casoratian,
positivity and size (§4); a Chebyshev-type bound for $\operatorname{lcm}(1..n)$ (§5); conclusion (§6).
No second series, no integrals, no prime number theorem.

**What to check first**
1. §3: the coefficient tables of the degree-4 identity $(\star)$.
2. §4: Lemma 4 (Casoratian $=24/n^3$), Lemma 5 ($|F_n|\le2$), Proposition 7.
3. §5: Lemma 8 (period-6 bound, $c=\frac12\log2+\frac34\log3$) and the final comparison $33.6<33.9$ in §6:
   the margin is small but strict.
4. Whether each motivation paragraph is accurate and sufficient.

**Originality:** the ingredients are classical (Apéry; Nesterenko's series; Zudilin's certificate). The
prime estimate of §5 was found by the other agent of this project (`codex/problem131/denominator-design.md`)
and re-derived here. What is ours is the arrangement.
**Scripts:** `zudilin_lemma1.py`, `zudilin_lemma2_recurrence.py`, `short_proof_checks.py`,
`companion_zeta2_forms.py`, `dn_bound_check.py`.

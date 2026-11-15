# Sources used in the paper and the slides

This file lists every outside work that the paper (`paper/paper.tex`) and the slides
(`presentation/slides.tex`) use or rely on. For each one it says what we took from it, where it
appears, and how the citation was checked.

**How the citations were checked (15 November 2026).** Bibliographic details were confirmed against
publisher, Numdam, EuDML, arXiv or library records found by web search. The full texts of arXiv
papers and of Zudilin's home page could not be opened from the session that compiled this list
(those hosts were blocked), so statements about their *content* come from search snippets or from
the team's own reading notes in `Peer Review/`. Items marked **(content not re-read)** should be
checked against the original before the paper is submitted.

## 1. Works whose results or constructions we use

| Work | What we use | Where |
|---|---|---|
| R. Apéry, *Irrationalité de ζ2 et ζ3*, Astérisque 61 (1979), 11–13. [Numdam](https://numdam.org/item/AST_1979__61__11_0/) | The theorem; the recurrence $(n+1)^3u_{n+1}-(2n+1)(17n^2+17n+5)u_n+n^3u_{n-1}=0$ and the numbers $a_n$. | Paper Part I intro, Parts II–III; slides throughout |
| F. Beukers, *A note on the irrationality of ζ(2) and ζ(3)*, Bull. London Math. Soc. 11 (1979), 268–272. [doi:10.1112/blms/11.3.268](https://doi.org/10.1112/blms/11.3.268) | The integral proof of Part I: kernel $-\log(xy)/(1-xy)$ (equivalently his triple integral), Legendre polynomials, positivity and the rate $(\sqrt2-1)^4$. **(content not re-read: the team never retrieved the paper itself; the exact form of his bound comes from the team's notes)** | Paper Part I; slides Part I |
| F. Beukers, *Legendre polynomials in irrationality proofs*, Bull. Austral. Math. Soc. 22 (1980), 431–438. [doi:10.1017/S0004972700006742](https://doi.org/10.1017/S0004972700006742) | The Legendre-polynomial / integration-by-parts point of view (Step 2). | Paper Part I |
| Yu. V. Nesterenko, *A few remarks on ζ(3)*, Mat. Zametki 59:6 (1996), 865–880; Math. Notes 59:6 (1996), 625–636. [Math-Net.Ru](https://mathnet.ru/eng/mzm1785) | The linear forms written as the sum $-\sum_t R_n'(t)$ (called Nesterenko's series in the team's Zudilin notes). **(attribution of this exact form from the team's notes)** | Paper Section 8 |
| W. Zudilin, *An elementary proof of Apéry's theorem*, arXiv:[math/0202159](https://arxiv.org/abs/math/0202159) (2002); published as *Apéry's theorem. Thirty years after*, Int. J. Math. Comput. Sci. 4 (2009), no. 1, 9–19. | The rational function $r_n(t)$ and $R_n=r_n^2$; the certificate $4(2n+1)(-2t^2+t+(2n+1)^2)R_n$ (our $W$); Ball's series, its degree-6 certificate, the $n=1$ check; the remark "using a trivial bound $D_n<3^n$". Formulas checked symbolically by us for small $n$. | Paper Sections 8–9, Part III; slides Part II, III |
| K. Ball and T. Rivoal, *Irrationalité d'une infinité de valeurs de la fonction zêta aux entiers impairs*, Invent. Math. 146 (2001), 193–207. [doi:10.1007/s002220100168](https://doi.org/10.1007/s002220100168) | The well-poised series type behind Ball's series. Zudilin credits the series and the observation that it approximates ζ(3) to K. Ball (unpublished). | Paper Section 9; slides "Ball's series" |
| T. Rivoal, *La fonction zêta de Riemann prend une infinité de valeurs irrationnelles aux entiers impairs*, C. R. Acad. Sci. Paris Sér. I 331 (2000), 267–270. [arXiv:math/0008051](https://arxiv.org/abs/math/0008051) | Background to Ball's series; per Zudilin, Rivoal found (with Zeilberger's program) that Ball's sequence coincides with Apéry's. | Paper Section 9 |
| D. Zeilberger, *The method of creative telescoping*, J. Symbolic Comput. 11 (1991), 195–204. [doi:10.1016/S0747-7171(08)80044-2](https://doi.org/10.1016/S0747-7171(08)80044-2) | Creative telescoping: how Zudilin's certificates are found by computer; our "guess" is a hand version. | Paper Section 9 |
| A. van der Poorten, *A proof that Euler missed … Apéry's proof of the irrationality of ζ(3). An informal report*, Math. Intelligencer 1 (1979), 195–203. [doi:10.1007/BF03028234](https://doi.org/10.1007/BF03028234) | Contemporary account of Apéry's proof; origin of the term "creative telescoping". | Paper Part I intro, Section 9 |
| P. L. Tchebichef (Chebyshev), *Mémoire sur les nombres premiers*, J. Math. Pures Appl. (1) 17 (1852), 366–390. [Numdam](https://numdam.org/item/JMPA_1852_1_17__366_0/) | Our identity $N!=\prod_k D_{\lfloor N/k\rfloor}$ is Chebyshev's identity $\log N!=\sum_k\psi(N/k)$ after taking logs. **(the identity is standard and due to Chebyshev; we did not confirm which of his two 1852 memoirs states it first)** | Paper Section 6; slides Step 5 and "Denominators" |
| A.-M. Legendre, *Essai sur la théorie des nombres*, 2nd ed., 1808. | Legendre's formula for the exponent of $p$ in $N!$ (proof of the identity above). | Paper Section 6 |
| D. Hanson, *On the product of the primes*, Canad. Math. Bull. 15 (1972), 33–37. [doi:10.4153/CMB-1972-007-7](https://doi.org/10.4153/CMB-1972-007-7) | $\mathrm{lcm}(1,\dots,n)<3^n$ for all $n$, as cited by later papers: the non-trivial fact behind Zudilin's "trivial bound". We do not need it. | Paper Section 9 |
| J. Hadamard, *Sur la distribution des zéros de la fonction ζ(s) et ses conséquences arithmétiques*, Bull. SMF 24 (1896), 199–220 [doi:10.24033/bsmf.545](https://doi.org/10.24033/bsmf.545); C.-J. de la Vallée Poussin, *Recherches analytiques sur la théorie des nombres premiers*, Ann. Soc. Sci. Bruxelles 20 (1896), 183–256. | The prime number theorem, which Beukers used and which our N! lemma replaces. (Page range of de la Vallée Poussin taken from secondary listings.) | Paper remark in Section 6 |

## 2. Prior work on the maximum-principle idea (cited so we claim no novelty)

| Work | Relation to Part III |
|---|---|
| A. B. Mingarelli, *Principal solutions of recurrence relations and irrationality questions in number theory*, arXiv:[math/0608577](https://arxiv.org/abs/math/0608577); Mathematics 11(2) (2023), 262, [doi:10.3390/math11020262](https://doi.org/10.3390/math11020262). | Treats Apéry's recurrence as a disconjugate three-term recurrence and the linear forms as its recessive (principal) solution; describes this as a restatement of Apéry's proof. Closest prior art to our maximum principle. |
| A. Mahboubi and T. Sibut-Pinote, *A formal proof of the irrationality of ζ(3)*, Log. Methods Comput. Sci. 17(1) (2021), 16:1–16:25, arXiv:[1912.06611](https://arxiv.org/abs/1912.06611). | Coq formalization; per the team's notes it proves $33^n=O(a_n)$ and uses $27<33$ (Lemma 5.6, **content not re-read**). |
| F. Chyzak, A. Mahboubi, T. Sibut-Pinote, E. Tassi, *A computer-algebra-based formal proof of the irrationality of ζ(3)*, ITP 2014, LNCS 8558, 160–176. [doi:10.1007/978-3-319-08970-6_11](https://doi.org/10.1007/978-3-319-08970-6_11) | Earlier conference version; recurrence certificates checked in Coq. |
| K. Rajkumar, *A simplification of Apéry's proof of the irrationality of ζ(3)*, arXiv:[1212.5881](https://arxiv.org/abs/1212.5881) (2012). | Gets the decay rate from the recurrence via Poincaré–Perron and links the construction to a Ramanujan continued fraction. |
| P. G. Doyle and J. L. Snell, *Random Walks and Electric Networks*, Carus Math. Monographs 22, MAA, 1984; arXiv:[math/0001057](https://arxiv.org/abs/math/0001057). | Standard reference for the voltage/current/Kirchhoff language used for the current law. |

## 3. Surveys and expositions (background, cited once)

- H. Cohen, *Démonstration de l'irrationalité de ζ(3) (d'après R. Apéry)*, Sém. Théorie des Nombres de Grenoble 6 (1977–1978), exposé 6, 1–9. [Numdam](https://www.numdam.org/item/STNG_1977-1978__6__A6_0/)
- S. Fischler, *Irrationalité de valeurs de zêta (d'après Apéry, Rivoal, …)*, Sém. Bourbaki exp. 910, Astérisque 294 (2004), 27–62. [arXiv:math/0303066](https://arxiv.org/abs/math/0303066)

Not cited, but useful further reading: D. Huylebrouck, *Similarities in irrationality proofs for π, ln 2, ζ(2) and ζ(3)*, Amer. Math. Monthly 108 (2001), 222–231; S. D. Miller, *An easier way to show ζ(3) ∉ ℚ* (1998, [note](https://ow3.math.rutgers.edu/~sdmiller/simplerzeta3.pdf)); G. Szegő, *Orthogonal Polynomials*, AMS Colloq. Publ. 23 (1939) and [DLMF §18](https://dlmf.nist.gov/18.3) for shifted Legendre polynomials; the Lean 4 ([arXiv:2503.07625](https://arxiv.org/abs/2503.07625)) and Isabelle ([AFP Zeta_3_Irrational](https://isa-afp.org/browser_info/current/AFP/Zeta_3_Irrational/document.pdf)) formalizations of Beukers' proof.

## 4. The team's own sources

- *Why ζ(3) is irrational: a proof built from its requirements* (Codex notes),
  `Peer Review/Zeta3-4-Motivated-proof-Codex/zeta3-motivated-proof.pdf` — Part I follows its Sections 1–8.
- *The irrationality of ζ(3): a maximum-principle proof* (reading edition, guess-and-check
  recurrence) — Part III follows it; earlier drafts are in `Peer Review/Digest/`.
- Zudilin exposition and referee notes, `Peer Review/Zeta3-1-Zudilin-exposition/exposition.md`.

All of the team's documents have been checked by AI only, not by a human referee.

## 5. Open items before submission

1. Read Beukers (1979) and confirm the form of his bound and his use of the prime number theorem.
2. Confirm in Zudilin's paper the attribution of Ball's series to K. Ball and the exact quote "trivial bound".
3. Confirm Lemma 5.6 of Mahboubi–Sibut-Pinote and Mingarelli's statement, since Part III's novelty remark rests on them.
4. Confirm which Chebyshev memoir to cite for $\log N!=\sum\psi(N/k)$.

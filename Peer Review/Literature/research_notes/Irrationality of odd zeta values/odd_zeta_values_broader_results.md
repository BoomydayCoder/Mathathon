# Irrationality and linear independence results for collections of odd zeta values ζ(3), ζ(5), ζ(7), … (2000 – 15 November 2026)

All retrievals were made on 15 November 2026. How each item was verified is marked as follows.

- **[arXiv]**: the arXiv record (and, where a statement is quoted, the arXiv PDF) was retrieved in this session.
- **[Zenodo]**: the Zenodo record was retrieved in this session.
- **[cited]**: the journal page could not be opened in this session (the hosts were blocked, see the last section). The bibliographic data come from the reference lists of retrieved arXiv papers or from web-search result snippets, and are labelled as such.

Nothing below is described from memory except the one-line entries explicitly marked "unverified (from memory)" in the last section.

## Status on 15 November 2026: current records, and the separation from "ζ(5) is irrational"

### Takeaway
None of the theorems of the types "infinitely many", "at least one of ζ(5),…,ζ(s)", "at least k of them" or "dimension ≥ …" proves that ζ(5), or any specific odd zeta value other than ζ(3), is irrational. The records I could verify are: dimension ≥ 0.21·√(s/log s) (Fischler); number of irrationals ≥ 2.78·√(s·log log s/log s) (Fischler–Lai–Sprang–Zudilin, arXiv preprint of 12 Nov 2026) or ≥ 1.284579·√(s/log s) (Lai, cited as published); one of ζ(5), ζ(7), ζ(9), ζ(11) (Zudilin 2001); two of ζ(5),…,ζ(35) (Lai–Zhou). Separately, a Zenodo preprint by A. Fauzan dated 17 September 2026 claims that ζ(5) itself is irrational, by a different (Hankel-determinant) method, and the newest papers in this bibliography cite it as a proof.

### Cited Findings

**Record table (each line verified on 15 Nov 2026 from the page linked)**

| Type of statement | Best result found | Status | Source |
|---|---|---|---|
| dim_Q Span(1, ζ(3), …, ζ(s)), s odd, large | ≥ 0.21·√(s/log s) | arXiv v3 of 27 Nov 2025; cited as J. Lond. Math. Soc. (2) 113 (2026), no. 4, e70535 | [arXiv:2109.10136](https://arxiv.org/abs/2109.10136); journal data as cited in [arXiv:2611.16156](https://arxiv.org/pdf/2611.16156) |
| Number of irrationals among ζ(3), …, ζ(s), s large | ≥ (√(4ζ(2)ζ(3)/ζ(6)) − ε)·√(s·log log s/log s), constant 2.788258… | arXiv preprint, 12 Nov 2026, not refereed | [arXiv:2611.16156](https://arxiv.org/abs/2611.16156) |
| Same, best result cited as published | ≥ 1.284579·√(s/log s) | arXiv 2025; cited as J. Théor. Nombres Bordeaux 38 (2026), no. 2, 469–487 | [arXiv:2501.05321](https://arxiv.org/abs/2501.05321) |
| At least one irrational in a short initial block | one of ζ(5), ζ(7), ζ(9), ζ(11) | published 2001 / 2004 | [arXiv:math/0206176](https://arxiv.org/abs/math/0206176) |
| At least two irrationals in a short initial block | two of ζ(5), ζ(7), …, ζ(35) | arXiv 2021; cited as Publ. Math. Debrecen 101 (2022), 353–372 | [arXiv:2103.00904](https://arxiv.org/abs/2103.00904) |
| 1, ζ(3), ζ(j) linearly independent for some odd j | j ≤ 139 (after 169, 145, 165) | Math. Ann. 347 (2010), cited | [arXiv:2103.00904 (introduction)](https://arxiv.org/pdf/2103.00904) |

- The Fischler–Lai–Sprang–Zudilin paper of 12 Nov 2026 still presents Zudilin's four-number theorem as the reference result of the "at least one" type: "This was refined by the fourth author [17]: at least one of ζ(5), ζ(7), ζ(9), ζ(11) is irrational. But the irrationality of ζ(5) itself was still very open." — [arXiv:2611.16156 (PDF, introduction)](https://arxiv.org/pdf/2611.16156)
- The same paper states that the dimension bound "was significantly improved by the first author [8]; this is however another story", and that it only improves the count of irrational values. — [arXiv:2611.16156 (PDF)](https://arxiv.org/pdf/2611.16156)
- Statements of the open status before September 2026, in the sources' own words:
  - "Malgré ces développements récents, il n'existe aucun entier s ≥ 5 impair pour lequel on sache si ζ(s) est rationnel ou non." (Fischler, Bourbaki exposé, 2002) — [arXiv:math/0303066](https://arxiv.org/pdf/math/0303066)
  - "the question whether ζ(5) is irrational remains open" (Fischler–Sprang–Zudilin, 2018) — [arXiv:1803.08905](https://arxiv.org/pdf/1803.08905)
  - "More than four decades after Apéry's breakthrough that ζ(3) ∉ Q, we still do not know any other single ζ(2k + 1) to be irrational." (Lai–Zhou, 2021) — [arXiv:2103.00904](https://arxiv.org/pdf/2103.00904)
  - "there is still no odd s ≥ 5 for which ζ(s) is known to be irrational" (Fischler, v3 of Nov 2025) — [arXiv:2109.10136](https://arxiv.org/pdf/2109.10136)
  - "So far, the irrationality of ζ(5) remains open" (Lai, Jan 2025) — [arXiv:2501.05321](https://arxiv.org/pdf/2501.05321)
- The separate claim about ζ(5) itself: the Zenodo record "ζ(5) is irrational" by Aabir Fauzan (Aalto University) has publication date 17 Sept 2026, resource type Preprint, version v1, DOI 10.5281/zenodo.22826419. Its description says the proof "constructs integer polynomials Qₙ of degree 37n with 0 < Qₙ(ζ(5)) < exp(−139n²/5)" that "are rationally normalized Hankel determinants", and gives |ζ(5) − a/b| > b^(−260) for large b. — [Zenodo record 22826419](https://zenodo.org/records/22826419)
- The same Zenodo record carries the note: "This version records the initial presentation of the main construction and argument, while the full manuscript is in preparation", and says later revisions "will be maintained on arXiv once the paper is announced there." — [Zenodo record 22826419](https://zenodo.org/records/22826419)
- An arXiv API author search for "Fauzan" restricted to 2026 returned five physics and computer-science papers and no paper on ζ(5), so the preprint was not on arXiv on 15 Nov 2026. — [arXiv API query](https://export.arxiv.org/api/query?search_query=au:Fauzan+AND+submittedDate:%5B202601010000+TO+202611160000%5D&max_results=20)
- Specialists cite the Fauzan preprint as a proof: "In a recent breakthrough, Fauzan proved the irrationality of ζ(5)" (Lai, 5 Oct 2026) — [arXiv:2610.06838](https://arxiv.org/abs/2610.06838); "the proof that ζ(5) is irrational posted by Fauzan on September 17, 2026" (Fischler, Lai, Sprang, Zudilin, 12 Nov 2026) — [arXiv:2611.16156](https://arxiv.org/abs/2611.16156)

### Inferences
- Lai–Zhou's "at least two of ζ(5),…,ζ(35)" implies unconditionally that at least one of ζ(7), ζ(9), …, ζ(35) is irrational. This is the sharpest verified statement about values beyond ζ(5), and it does not depend on the Fauzan preprint.
- If the Fauzan preprint is correct, "one of ζ(5), ζ(7), ζ(9), ζ(11)" is subsumed, but none of the counting or dimension theorems is. The report should therefore present two separate lists: (a) theorems about collections, and (b) the single claimed theorem about ζ(5).
- The only refereed-journal statuses I could attach to 2026 items are those given in other authors' reference lists. I found no evidence that the Fauzan preprint has been refereed.

### Gaps
- I did not read or check the Fauzan PDF. Its validity is outside this assignment; only its existence, date, stated result and how others cite it are recorded.
- A web-search summary reported that Wikipedia lists the Fauzan preprint as a "claimed proof" with a claimed Lean formalization, and that MathWorld mentions a formal verification by Romik assuming the prime number theorem (GitHub danromik/zeta5-irrationality). I could not open either page (hosts blocked), so these are unverified search snippets.
- A web-search result also listed a page titled "ζ(5) is irrational" on arxivlens.com. I could not determine which paper it mirrors.

## Rivoal (C. R. 2000), Ball–Rivoal (Invent. Math. 2001), Nesterenko's criterion and its refinements

### Takeaway
Rivoal and Ball–Rivoal proved dim_Q Span(1, ζ(3), …, ζ(s)) ≥ ((1−ε)/(1+log 2))·log s using a well-poised hypergeometric series, whose symmetry removes the even zeta values, together with Nesterenko's 1985 linear independence criterion. The result gives no information about which odd zeta values are irrational. Later work refined the criterion (Fischler–Zudilin, Fischler) or replaced it by Siegel's criterion plus a Shidlovsky-type zero estimate.

### Cited Findings

**T. Rivoal, "La fonction zêta de Riemann prend une infinité de valeurs irrationnelles aux entiers impairs"** (7 Aug 2000) [arXiv]
- Links: [arXiv:math/0008051](https://arxiv.org/abs/math/0008051), [PDF](https://arxiv.org/pdf/math/0008051). Type: research announcement, 4 pages. Published as C. R. Acad. Sci. Paris Sér. I Math. 331 (2000), no. 4, 267–270, DOI 10.1016/S0764-4442(00)01624-4 (DOI from the arXiv record; pages as cited in [arXiv:2611.16156](https://arxiv.org/pdf/2611.16156)).
- Result (Théorème 1): for every ε > 0 there is N(ε) such that for n > N(ε), dim_Q(Q + Qζ(3) + ⋯ + Qζ(2n+1)) ≥ ((1−ε)/(1+log 2))·log n. — [arXiv:math/0008051 (PDF)](https://arxiv.org/pdf/math/0008051)
- Method: the note says the proof "s'inspire du travail de Nikishin" on type I Padé approximants of polylogarithms. — [arXiv:math/0008051 (PDF)](https://arxiv.org/pdf/math/0008051)
- Limitation: the theorem is asymptotic and gives a dimension, not the irrationality of any named value.

**K. Ball, T. Rivoal, "Irrationalité d'une infinité de valeurs de la fonction zêta aux entiers impairs"** (2001) [cited]
- Reference: Invent. Math. 146 (2001), no. 1, 193–207, as cited in [arXiv:2611.16156](https://arxiv.org/pdf/2611.16156) and [arXiv:2109.10136](https://arxiv.org/pdf/2109.10136). No arXiv version was found; the Springer page could not be opened.
- Statements, as reported in Fischler's Bourbaki exposé — [arXiv:math/0303066 (PDF, §2.1)](https://arxiv.org/pdf/math/0303066):
  - Théorème 2.1: for odd ℓ ≥ 3 let δ_ℓ be the dimension of the Q-span of 1, ζ(3), ζ(5), …, ζ(ℓ). For every ε > 0 there is ℓ₀ with δ_ℓ ≥ ((1−ε)/(1+log 2))·log ℓ for all odd ℓ ≥ ℓ₀.
  - Remarque 2.2: with 1/3 in place of (1−ε)/(1+log 2) one can take ℓ₀ = 3.
  - Théorème 2.3: there is an odd ℓ ≤ 169 such that 1, ζ(3), ζ(ℓ) are linearly independent over Q.
- Method, as described by Fischler–Sprang–Zudilin: the series n!^(s−2r)·Σ_{t≥1} Π_{j=0}^{(2r+1)n}(t − rn + j) / Π_{j=0}^{n}(t + j)^(s+1) "happens to be a Q-linear combination of 1 and odd zeta values when s is odd and n is even", and it is combined with Nesterenko's criterion. — [arXiv:1803.08905 (PDF)](https://arxiv.org/pdf/1803.08905)
- Fischler's exposé names the two essential ingredients: "l'absence de ζ(2), ζ(4), …, ζ(ℓ − 1) d'une part, et la minoration en log(ℓ) de la dimension d'autre part." — [arXiv:math/0303066 (PDF)](https://arxiv.org/pdf/math/0303066)
- Why the log bound resisted improvement: "To improve on this bound using the same strategy, one has to find linear combinations that are considerably smaller, with not too large coefficients, — it comes out to be a rather difficult task." — [arXiv:1803.08905 (PDF)](https://arxiv.org/pdf/1803.08905)
- Origin of the series: Zudilin reports that K. Ball's well-poised series "produces linear forms in 1 and ζ(3) only and can be evaluated by elementary means; however, its 'obvious' arithmetic does not allow one to prove the irrationality of ζ(3)", and that Rivoal generalised it in the spirit of Nikishin. — [arXiv:math/0206176 (PDF)](https://arxiv.org/pdf/math/0206176)

**Yu. V. Nesterenko, "On the linear independence of numbers"** (1985) [cited]
- Reference: Vestnik Moskov. Univ. Ser. I Mat. Mekh. 40 (1985), no. 1, 46–49; English translation Moscow Univ. Math. Bull., 69–74 — as cited in [arXiv:2109.10136 (reference 23)](https://arxiv.org/pdf/2109.10136). Not retrieved.
- Role: the criterion turns a sequence of small integer linear forms, with controlled coefficient growth, into a lower bound on the dimension. Fischler's abstract notes that it needs a lower bound on the linear forms, which Siegel's criterion avoids. — [arXiv:1609.09770](https://arxiv.org/abs/1609.09770)

**S. Fischler, W. Zudilin, "A refinement of Nesterenko's linear independence criterion with applications to zeta values"** (2010) [cited]
- Reference: Math. Ann. 347 (2010), no. 4, 739–763, as cited in [arXiv:2610.06838 (reference 21)](https://arxiv.org/pdf/2610.06838). No arXiv version was found by title search; the Springer page could not be opened.
- Result as reported by others: "Fischler and Zudilin [FZ10] lowered 169 to 139 by refining Nesterenko's linear independence criterion" — [arXiv:2103.00904 (PDF)](https://arxiv.org/pdf/2103.00904). Rivoal–Zudilin phrase it as two distinct odd m, n ∈ {3, 5, …, 139} with 1, ζ(m), ζ(n) linearly independent — [arXiv:1803.03160 (PDF)](https://arxiv.org/pdf/1803.03160).
- A web-search summary of the Springer abstract says the refinement exploits large common divisors of some coefficients of the linear forms; this is a snippet, not a retrieved page.

**S. Fischler, "Nesterenko's criterion when the small linear forms oscillate"** (12 Jan 2012) [arXiv]
- Links: [arXiv:1201.2651](https://arxiv.org/abs/1201.2651). Status: "To appear in Archiv der Math" (arXiv comment).
- Result: a version of the criterion for linear forms with oscillating behaviour (as produced by the saddle-point method), giving a dimension bound and a measure of simultaneous approximation. Application: "an explicit measure of simultaneous approximation to ζ(5), ζ(7), ζ(9), and ζ(11), using Zudilin's proof that at least one of these numbers is irrational." — [arXiv:1201.2651](https://arxiv.org/abs/1201.2651)

**S. Fischler, "Nesterenko's linear independence criterion for vectors"** (10 Feb 2012) [arXiv]
- Links: [arXiv:1202.2279](https://arxiv.org/abs/1202.2279). Published: Monatsh. Math. 177 (2015), 397–419 (arXiv journal reference).
- Result: a lower bound for the rank of p vectors in R^k from linear forms small at k points; the case k = 1 is Nesterenko's criterion. Proof by Minkowski's convex body theorem. — [arXiv:1202.2279](https://arxiv.org/abs/1202.2279)

**S. Fischler, "Distribution of irrational zeta values"** (7 Oct 2013) [arXiv]
- Links: [arXiv:1310.1685](https://arxiv.org/abs/1310.1685). Published: Bull. Soc. Math. France 145 (2017), no. 3, 381–409 (as cited in [arXiv:2610.06838](https://arxiv.org/pdf/2610.06838)).
- Result: for odd a large in terms of ε there are [(1−ε)·log a/(1+log 2)] odd integers s between 3 and a, at distance at least a^ε from one another, at which ζ takes Q-linearly independent values. Method: hypergeometric-type series, saddle point, criterion for vectors. — [arXiv:1310.1685](https://arxiv.org/abs/1310.1685)

**S. Fischler, "Shidlovsky's multiplicity estimate and Irrationality of zeta values"** (30 Sep 2016) [arXiv]
- Links: [arXiv:1609.09770](https://arxiv.org/abs/1609.09770). Published: J. Austral. Math. Soc. 105 (2018), no. 2, 145–172 (as cited in [arXiv:2109.10136](https://arxiv.org/pdf/2109.10136)).
- Result: a general Shidlovsky lemma for Padé problems at several points, giving "a new proof of the Ball-Rivoal theorem … Shidlovsky's lemma is used to replace Nesterenko's linear independence criterion with Siegel's, so that no lower bound is needed on the linear forms in zeta values." — [arXiv:1609.09770](https://arxiv.org/abs/1609.09770)
- Connection: this is the technical precursor of Fischler's 2021 Siegel-lemma theorem (see the 2018–2026 section).

**Related criterion papers** [arXiv]
- Fischler, Hussain, Kristensen, Levesley, "A converse to linear independence criteria, valid almost everywhere": proves optimality of several linear independence criteria. — [arXiv:1302.1952](https://arxiv.org/abs/1302.1952)
- Fischler, Rivoal, "Irrationality exponent and rational approximations with prescribed growth" (Proc. Amer. Math. Soc. 138 (2010)): discusses multivariate questions "connected to Fischler-Zudilin's new proof of Nesterenko's linear independence criterion". — [arXiv:0910.4448](https://arxiv.org/abs/0910.4448)
- Cresson, Fischler, Rivoal, "Phénomènes de symétrie dans des formes linéaires en polyzêtas": generalises, in arbitrary depth, "the symmetry phenomenon used by Ball-Rivoal". — [arXiv:math/0609744](https://arxiv.org/abs/math/0609744)

### Inferences
- The Ball–Rivoal method has two separable parts: a parity (well-poised symmetry) argument that removes even zeta values, and a size-versus-denominator comparison that fixes the growth rate. Every later result in this bibliography modifies one of the two.
- Linear independence criteria are needed only for dimension statements. The "at least one" and counting theorems avoid them, which is why they could be made elementary in 2018.

### Gaps
- The Invent. Math. paper, the Math. Ann. paper and Nesterenko's 1985 note were not retrieved; their statements are taken from surveys and later papers.
- I did not find lecture notes by Fischler devoted specifically to linear independence criteria (see the surveys section for the one related item found by search).

## "At least one of ζ(5), …, ζ(s)" and other small-range results, 2001–2012 (Rivoal, Zudilin, Krattenthaler–Rivoal)

### Takeaway
Rivoal obtained one irrational among ζ(5),…,ζ(21) by differentiating a very-well-poised series twice to remove ζ(3). Zudilin reduced the range to ζ(5),…,ζ(11) by a more general construction and a p-adic saving factor Φ_n, and Krattenthaler–Rivoal proved the "denominator conjecture", which gives ζ(5),…,ζ(19) without Φ_n. Zudilin's four-number result holds with a very small numerical margin, and the sources say that only a further, conjectural, drop in denominators could reach three numbers.

### Cited Findings

**T. Rivoal, "Irrationalité d'au moins un des neuf nombres ζ(5), ζ(7), …, ζ(21)"** (24 Apr 2001) [arXiv]
- Links: [arXiv:math/0104221](https://arxiv.org/abs/math/0104221), [PDF](https://arxiv.org/pdf/math/0104221). Published: Acta Arith. 103 (2002), no. 2, 157–167, DOI 10.4064/aa103-2-5 (DOI from the arXiv record; pages as cited in [arXiv:2611.16156](https://arxiv.org/pdf/2611.16156)).
- Result: at least one of the nine numbers ζ(5), ζ(7), …, ζ(21) is irrational. — [arXiv:math/0104221](https://arxiv.org/abs/math/0104221)
- Method: a series with the extra factor (t + n/2), differentiated twice, is "a Q-linear combination of 1 and odd zeta values starting from ζ(5)" (used with s = 20) — [arXiv:1803.08905 (PDF)](https://arxiv.org/pdf/1803.08905). The linear forms are estimated by the saddle-point method following Nesterenko — [arXiv:math/0104221 (PDF)](https://arxiv.org/pdf/math/0104221). In the Bourbaki account, with a = 20 the integer forms d_{φ(n)}^22·S_{φ(n)}(1) in 1, ζ(5), …, ζ(21) tend to 0 and are nonzero — [arXiv:math/0303066 (PDF, §3.2)](https://arxiv.org/pdf/math/0303066).
- No linear independence criterion is needed: "Proofs of these results do not require use of linear independence criteria". — [arXiv:1803.08905 (PDF)](https://arxiv.org/pdf/1803.08905)
- Documented limitation: "Si on savait démontrer la conjecture … on pourrait appliquer la même méthode avec a = 18, et démontrer ainsi que l'un au moins des nombres ζ(5), ζ(7), …, ζ(19), est irrationnel" (Remarque 3.3). — [arXiv:math/0303066 (PDF)](https://arxiv.org/pdf/math/0303066)

**W. Zudilin, "Irrationality of values of zeta-function"** (26 Apr 2001) [arXiv]
- Links: [arXiv:math/0104249](https://arxiv.org/abs/math/0104249). Type: 8+8 page conference paper (English and Russian). The arXiv record carries the DOI of the full paper "Irrationality of values of the Riemann zeta function", Izv. Math. 66 (2002), no. 3, 489–542 (DOI 10.1070/IM2002v066n03ABEH000387; journal data as cited in [arXiv:2610.06838](https://arxiv.org/pdf/2610.06838)).
- Results — [arXiv:math/0104249 (PDF)](https://arxiv.org/pdf/math/0104249):
  - Theorem 1: each of {ζ(5), …, ζ(21)}, {ζ(7), …, ζ(37)}, {ζ(9), …, ζ(53)} contains an irrational number.
  - Theorem 2: for each odd b ≥ 1, the collection ζ(b+2), ζ(b+4), …, ζ(8b−1) contains an irrational number.
  - Theorem 3: there are odd a₁ ≤ 145 and a₂ ≤ 1971 such that 1, ζ(3), ζ(a₁), ζ(a₂) are linearly independent over Q.
  - Theorem 4: for each odd a ≥ 3, δ(a) > 0.395·log a > (2/3)·log a/(1+log 2).
- Method: a generalisation of Rivoal's construction with parameters (a, b, c), saddle-point asymptotics, and improved denominator estimates. A footnote says Rivoal obtained the first collection of Theorem 1 independently. — [arXiv:math/0104249 (PDF)](https://arxiv.org/pdf/math/0104249)
- Theorem 2 is quoted later in the equivalent form (b = 2n + 1): "for any integer n ≥ 0, there exists at least one irrational number amongst … ζ(2n + 3), ζ(2n + 5), …, ζ(16n + 7)", and is used with n = 5 in the Rivoal–Zudilin proof. — [arXiv:1803.03160 (PDF)](https://arxiv.org/pdf/1803.03160)

**W. Zudilin, "One of the numbers ζ(5), ζ(7), ζ(9), ζ(11) is irrational"** (2001) [cited]
- Reference: Uspekhi Mat. Nauk 56 (2001), no. 4, 149–150; translation Russian Math. Surveys 56 (2001), no. 4, 774–776 — as cited in [arXiv:2610.06838 (reference 60)](https://arxiv.org/pdf/2610.06838). The mathnet.ru page could not be opened; a search snippet gives DOI 10.1070/RM2001v056n04ABEH000427.
- Result: the statement in the title. The detailed proof is in the next entry.

**W. Zudilin, "Arithmetic of linear forms involving odd zeta values"** (18 Jun 2002) [arXiv]
- Links: [arXiv:math/0206176](https://arxiv.org/abs/math/0206176), [PDF](https://arxiv.org/pdf/math/0206176). Published: J. Théor. Nombres Bordeaux 16 (2004), no. 1, 251–291 (arXiv journal reference). 42 pages.
- Result: "A general hypergeometric construction of linear forms in (odd) zeta values … allows to recover the records of Rhin and Viola for the irrationality measures of ζ(2) and ζ(3), as well as to explain Rivoal's 'infinitely-many' result … and to prove that at least one of the four numbers ζ(5), ζ(7), ζ(9), and ζ(11) is irrational." — [arXiv:math/0206176](https://arxiv.org/abs/math/0206176)
- Method and the size of the margin (Bourbaki §3.3): a p-adic analysis of binomial coefficients gives an integer Φ_n that divides the coefficients; the asymptotics are limsup |S̃_n(1)|^(1/n) ≤ e^(−227.58…), limsup |Φ_n^(−1)|^(1/n) ≤ e^(−176.75…), and denominators growing like e^(403n); "Comme 403 < 227,58 + 176,75 on obtient la conclusion cherchée." — [arXiv:math/0303066 (PDF)](https://arxiv.org/pdf/math/0303066)
- Useful sections: §7 (arithmetic of the rational-function "bricks"), §8 (well-poised origin of Rivoal's construction and the four-number theorem), §9 (possible group structures, conditional on an explicit arithmetic condition). — [arXiv:math/0206176 (PDF, introduction)](https://arxiv.org/pdf/math/0206176)
- Documented limitation: "Zudilin conjecture ([Zu4], §9) que des compensations ont lieu quand z = 1, ce qui permettrait de trouver un dénominateur plus petit pour les P_j(1). Peut-être pourrait-on alors démontrer que parmi ζ(5), ζ(7) et ζ(9) l'un au moins est irrationnel" (Remarque 3.4). — [arXiv:math/0303066 (PDF)](https://arxiv.org/pdf/math/0303066)

**W. Zudilin, "Well-poised hypergeometric service for diophantine problems of zeta values"** (2003) [cited]
- Reference: J. Théor. Nombres Bordeaux 15 (2003), no. 2, 593–626 (Numdam and EuDML records found by web search; not opened). The search summary lists: a permutation group for linear forms in 1 and ζ(4), a second-order Apéry-like recursion for ζ(4), low-order recursions for linear forms in odd zeta values, and a permutation group for Euler-type multiple integrals.
- Retrieved companion: "Multiple-integral representations of very-well-poised hypergeometric series" — "The multiple integrals are of Euler-type and can be viewed as a natural generalization of Beukers' integrals for ζ(2) and ζ(3)." — [arXiv:math/0206177](https://arxiv.org/abs/math/0206177)

**C. Krattenthaler, T. Rivoal, "Hypergéométrie et fonction zêta de Riemann"** (7 Nov 2003, v4 20 Dec 2004) [arXiv]
- Links: [arXiv:math/0311114](https://arxiv.org/abs/math/0311114), [PDF](https://arxiv.org/pdf/math/0311114). Published: Mem. Amer. Math. Soc. 186 (2007), no. 875 (arXiv journal reference). 73 pages.
- Result: proof of Rivoal's "denominator conjecture" on common denominators of the coefficients of the linear forms. Corollaries: "at least one of the eight numbers ζ(5), ζ(7), …, ζ(19) is irrational, and there exists an odd integer j between 5 and 165 such that 1, ζ(3) and ζ(j) are linearly independent over Q." Also proofs of a conjecture of Vasilyev and of Zudilin's conjecture on rational approximations of ζ(4). — [arXiv:math/0311114](https://arxiv.org/abs/math/0311114)
- Method: "a hypergeometric identity between a single sum and a multiple sum due to Andrews". — [arXiv:math/0311114](https://arxiv.org/abs/math/0311114)
- What the conjecture says (Bourbaki, Conjecture 2.6 and Remarque 2.13): in the integrality statement one can replace d_n^a by d_n^(a−1), because "pour z = 1 on a des compensations particulières qui font chuter le dénominateur", while examples show such a gain is impossible for general z. — [arXiv:math/0303066 (PDF)](https://arxiv.org/pdf/math/0303066)
- Documented limitation: "We hope that it will be possible to apply our construction to the more general linear forms constructed by Zudilin [56], with the ultimate goal of strengthening his result that one of the numbers ζ(5), ζ(7), ζ(9), ζ(11) is irrational." — [arXiv:math/0311114](https://arxiv.org/abs/math/0311114)
- The arXiv comment records that v4 was "completely rewritten" so that "our theorems now hold unconditionally". — [arXiv:math/0311114](https://arxiv.org/abs/math/0311114)

**Progress of the "1, ζ(3), ζ(j)" statement**
- 169 (Ball–Rivoal), then 145 (Zudilin, "grâce à un raffinement du lemme 2.12") — [arXiv:math/0303066 (PDF)](https://arxiv.org/pdf/math/0303066); 165 by the denominator conjecture alone — [arXiv:math/0311114](https://arxiv.org/abs/math/0311114); 139 (Fischler–Zudilin 2010) — [arXiv:2103.00904 (PDF)](https://arxiv.org/pdf/2103.00904).

### Inferences
- Zudilin's four-number proof closes with 403 < 404.33, a margin of about 0.3 % of the exponent. This explains why no purely numerical re-optimisation of the same construction has reached three numbers.
- The denominator conjecture (now a theorem for Rivoal's forms) and the Φ_n factor are two independent sources of arithmetic saving. The sources indicate that combining both for Zudilin's general forms would be the natural route to "one of ζ(5), ζ(7), ζ(9)", and that this combination was not available in the verified literature.
- All results in this section produce linear forms in several zeta values at once. They cannot isolate ζ(5), because the constructions that remove ζ(3) leave at least four odd zeta values with nonzero coefficients.

### Gaps
- A preprint claims the three-number statement: S. Suman, "At least one of the three numbers ζ(5), ζ(7), ζ(9) is irrational" (15 Sep 2026, category math.GM, 6 pages) — [arXiv:2609.22316](https://arxiv.org/abs/2609.22316). I did not check the proof. Reasons for caution: the same author's two earlier ζ(5) submissions carry the author's own comments "Found incorrect" — [arXiv:2407.07121](https://arxiv.org/abs/2407.07121) — and "Not a good article as of now. Needs improvement" — [arXiv:2312.00298](https://arxiv.org/abs/2312.00298); and the specialists' paper posted two months later still gives the four-number theorem as the state of the art and does not cite it — [arXiv:2611.16156 (PDF)](https://arxiv.org/pdf/2611.16156). It should not be reported as a record.
- I could not retrieve the Uspekhi note or the 2003 "service" paper themselves.

## The 2018–2026 wave: elementary proofs, elimination, and the growth-rate improvements

### Takeaway
In 2018 Zudilin's "twist by half" gave an elementary proof (prime number theorem and Stirling only) that one of ζ(5),…,ζ(25) is irrational, which is weaker in range than his 2001 theorem. Sprang's general twists and the Fischler–Sprang–Zudilin elimination technique then raised the number of irrational values among ζ(3),…,ζ(s) from about log s to 2^((1−ε)·log s/log log s), Lai–Yu to c·√(s/log s), and the November 2026 preprint to c·√(s·log log s/log s). Fischler's Siegel-lemma method is a genuinely different, non-explicit construction and is the only one that improved the dimension bound.

### Cited Findings

**W. Zudilin, "One of the Odd Zeta Values from ζ(5) to ζ(25) Is Irrational. By Elementary Means"** (30 Jan 2018) [arXiv]
- Links: [arXiv:1801.09895](https://arxiv.org/abs/1801.09895), [PDF](https://arxiv.org/pdf/1801.09895). Published: SIGMA 14 (2018), 028, 8 pages, DOI 10.3842/SIGMA.2018.028.
- Result: at least one of ζ(5), ζ(7), …, ζ(25) is irrational, by "a proof whose technical ingredients are limited to the prime number theorem and Stirling's approximation formula". With d_n < 3^n in place of the prime number theorem, the same argument works with s = 33. — [arXiv:1801.09895 (PDF)](https://arxiv.org/pdf/1801.09895)
- Method: two well-poised series r_n = Σ R_n(ν) and r̂_n = Σ R_n(ν − 1/2) whose ζ(3)-coefficients are proportional with factor 7, so that 7r_n − r̂_n lies in Q + Qζ(5) + ⋯ + Qζ(s); for s = 25, lim (7r_n − r̂_n)^(1/n) = exp(−25.292363…). — [arXiv:1801.09895 (PDF)](https://arxiv.org/pdf/1801.09895)
- Relation to earlier work: a simplification, not a new record. Sprang calls it "a slightly weaker form of his important result". — [arXiv:1802.09410](https://arxiv.org/abs/1802.09410)
- Documented limitation: the new elimination idea "can be further used with the arithmetic method in [10] to significantly reduce the size of s. Since this does not let s be down to s = 9, hence leaving the achievement 'at least one of the four numbers ζ(5), ζ(7), ζ(9), ζ(11) is irrational' unchanged, we do not discuss this generalisation in greater details." — [arXiv:1801.09895 (PDF)](https://arxiv.org/pdf/1801.09895)
- Companion: Zudilin, "Some hypergeometric integrals for linear forms in zeta values" (Bull. Austral. Math. Soc. 98 (2018), 372–375) proves integral representations of the forms of arXiv:1801.09895 and arXiv:1803.08905. — [arXiv:1804.04129](https://arxiv.org/abs/1804.04129)

**J. Sprang, "Infinitely many odd zeta values are irrational. By elementary means"** (26 Feb 2018) [arXiv]
- Links: [arXiv:1802.09410](https://arxiv.org/abs/1802.09410). Status: arXiv preprint, 11 pages; still cited as "Preprint (2018)" in November 2026 — [arXiv:2611.16156 (PDF)](https://arxiv.org/pdf/2611.16156).
- Result and method: an elementary proof that infinitely many odd zeta values are irrational, generalising Zudilin's "twists by half" to "higher twists". — [arXiv:1802.09410](https://arxiv.org/abs/1802.09410)
- Limitation: it reproves the qualitative "infinitely many" statement, not the dimension bound of Ball–Rivoal.

**S. Fischler, J. Sprang, W. Zudilin, "Many odd zeta values are irrational"** (23 Mar 2018) [arXiv]
- Links: [arXiv:1803.08905](https://arxiv.org/abs/1803.08905), [PDF](https://arxiv.org/pdf/1803.08905). Published: Compositio Math. 155 (2019), 938–952, DOI 10.1112/S0010437X1900722X. 14 pages.
- Result: at least 2^((1−ε)·log s/log log s) of the values ζ(3), ζ(5), …, ζ(s) are irrational for s large in terms of ε. "This lower bound is asymptotically larger than any power of log s; it improves on the bound ((1−ε)/(1+log 2))·log s that follows from the Ball–Rivoal theorem." — [arXiv:1803.08905](https://arxiv.org/abs/1803.08905)
- Method: "construction of several linear forms in odd zeta values with related coefficients", so that unwanted zeta values can be eliminated. — [arXiv:1803.08905](https://arxiv.org/abs/1803.08905)
- Announcement: Fischler, Sprang, Zudilin, "Many values of the Riemann zeta function at odd integers are irrational", C. R. Math. Acad. Sci. Paris 356 (2018), no. 7, 707–711 [cited] — as cited in [arXiv:2611.16156 (reference 9)](https://arxiv.org/pdf/2611.16156). The authors of the C. R. note are Fischler, Sprang and Zudilin, not Lai.
- Limitation: "This strategy provides only a lower bound on the number of irrational odd zeta values, but nothing like (1.1) or Theorem 1 about linear independence." — [arXiv:2109.10136 (PDF)](https://arxiv.org/pdf/2109.10136)

**C. Krattenthaler, W. Zudilin, "Hypergeometry inspired by irrationality questions"** (24 Feb 2018) [arXiv]
- Links: [arXiv:1802.08856](https://arxiv.org/abs/1802.08856). Published: Kyushu J. Math. 73:1 (2019), 189–203.
- Result (Theorem 1): for any real λ, each of two explicit collections of numbers of the form ζ(2m+1) − λ·c_m·π^(2m+1) (with rational c_m given by Bernoulli numbers; m = 1,…,19 for the first and m = 1,…,21 for the second) contains at least one irrational number. The authors compare it with a result implying the same for m = 1,…,169. — [arXiv:1802.08856 (PDF)](https://arxiv.org/pdf/1802.08856)
- Relevance: an application of the twist idea to mixed forms in odd zeta values and odd powers of π; it says nothing about ζ(5) alone.

**T. Rivoal, W. Zudilin, "A note on odd zeta values"** (8 Mar 2018) [arXiv]
- Links: [arXiv:1803.03160](https://arxiv.org/abs/1803.03160), [PDF](https://arxiv.org/pdf/1803.03160). Published: Séminaire Lotharingien de Combinatoire 81 (2020), Art. B81b, 13 pp.
- Result: Theorem 1: at least two irrational numbers among the 33 values ζ(5), ζ(7), …, ζ(69). Equivalent Theorem 2: for each m with 1 ≤ m ≤ 33 there is an irrational number among the 32 values ζ(2n+3), n ∈ {1,…,33}∖{m}. — [arXiv:1803.03160 (PDF)](https://arxiv.org/pdf/1803.03160)
- Method: Rivoal's twice-differentiated series combined with the half-twist; "the additional factor (t − n + 1/2)_{3n} is the principal innovation"; the saddle-point method is used. — [arXiv:1803.03160 (PDF)](https://arxiv.org/pdf/1803.03160)
- Documented limitations: "our method can't decide if 1 and the two irrational odd zeta values in Theorem 1 are linearly independent over Q." The paper also notes the alternative localisation: two irrational values ζ(m), ζ(n) with m ∈ {5, 7, 9, 11} and n ∈ {13, 15, …, 87}. — [arXiv:1803.03160 (PDF)](https://arxiv.org/pdf/1803.03160)

**L. Lai, P. Yu, "A note on the number of irrational odd zeta values"** (19 Nov 2019) [arXiv]
- Links: [arXiv:1911.08458](https://arxiv.org/abs/1911.08458). Published: Compositio Math. 156 (2020), 1699–1717, DOI 10.1112/S0010437X20007307.
- Result: for odd s ≥ s₀(ε) there are at least (c₀ − ε)·√(s/log s) irrational numbers among ζ(3), ζ(5), …, ζ(s), with c₀ = 1.192507… in closed form. — [arXiv:1911.08458](https://arxiv.org/abs/1911.08458)
- Method: the Fischler–Sprang–Zudilin scheme with "an optimal design for the zeros of the auxiliary rational functions, which relates to the inverse of Euler totient function". — [arXiv:1911.08458](https://arxiv.org/abs/1911.08458)
- Relation: a refinement of the construction inside the same elimination framework; it changes the growth rate, not only the constant. The arXiv comment records that v2 "improved the constant 1/10 to about 1.19". — [arXiv:1911.08458](https://arxiv.org/abs/1911.08458)

**L. Lai, L. Zhou, "At least two of ζ(5), ζ(7), …, ζ(35) are irrational"** (1 Mar 2021) [arXiv]
- Links: [arXiv:2103.00904](https://arxiv.org/abs/2103.00904), [PDF](https://arxiv.org/pdf/2103.00904). Published: Publ. Math. Debrecen 101 (2022), no. 3–4, 353–372 (as cited by Lai in [arXiv:2610.06838](https://arxiv.org/pdf/2610.06838)).
- Result: the statement in the title, "which improves a recent result of Rivoal and Zudilin by lowering 69 to 35"; also at least one of β(2), β(4), …, β(10) is irrational. — [arXiv:2103.00904](https://arxiv.org/abs/2103.00904)
- Method: "To rule out the unwanted ζ(3), there are two different methods. The first method is 'taking twice derivatives' … the second one is 'inserting rational roots' … In [RZ20], the authors combined these two methods; in this paper, we use purely the second method in an elaborated way." The construction combines Zudilin 2004 and Lai–Yu, with shifts θ ∈ {1, 1/2, 1/3, 2/3} and Φ_n factors. — [arXiv:2103.00904 (PDF)](https://arxiv.org/pdf/2103.00904)
- Relation: a refinement of parameters and of the elimination device, not a new principle. The authors remark that Rivoal and Zudilin "did not pursue the full strength of their method for simplicity". — [arXiv:2103.00904 (PDF)](https://arxiv.org/pdf/2103.00904)

**S. Fischler, "Linear independence of odd zeta values using Siegel's lemma"** (21 Sep 2021, v3 27 Nov 2025) [arXiv]
- Links: [arXiv:2109.10136](https://arxiv.org/abs/2109.10136), [PDF](https://arxiv.org/pdf/2109.10136). Published: J. London Math. Soc. (2) 113 (2026), no. 4, e70535, 42 pp. (as cited in [arXiv:2611.16156](https://arxiv.org/pdf/2611.16156); the Wiley page returned an access challenge and a web search did not confirm the journal data independently). 44 pages.
- Result: Theorem 1: for sufficiently large odd s, dim_Q Span(1, ζ(3), ζ(5), …, ζ(s)) ≥ 0.21·√(s/log s); "the first asymptotic improvement on the lower bound, logarithmic in s, obtained by Ball-Rivoal in 2001". Theorem 2: for algebraic z with |z| ≤ 1, z ∉ {0, 1}, dim over Q(z) of Span(1, Li₁(z), …, Li_s(z)) ≥ (0.26/[Q(z):Q])·√(s/log s). — [arXiv:2109.10136 (PDF)](https://arxiv.org/pdf/2109.10136)
- Method: "Siegel's lemma to construct non-explicit linear forms … instead of using explicit well-poised hypergeometric series. A new refinement of Siegel's linear independence criterion is applied, together with a multiplicity estimate (namely a generalization of Shidlovsky's lemma)." This is a genuinely different construction from all hypergeometric ones. — [arXiv:2109.10136](https://arxiv.org/abs/2109.10136)
- Documented correction: the v3 comment says the new criterion "is used to correct an error pointed out by a referee and Ludovic Mistiaen. It turns out that the coefficient matrix of the linear forms we construct does not have always maximal rank: it satifies only a weaker assumption". — [arXiv:2109.10136](https://arxiv.org/abs/2109.10136)
- Limitation: asymptotic only; "0.21 is the rounded value of a real number that we did not try to compute exactly"; nothing for small s. — [arXiv:2109.10136 (PDF)](https://arxiv.org/pdf/2109.10136)

**L. Lai, "Small improvements to the Ball–Rivoal theorem and its p-adic variant"** (19 Jul 2024, v3 25 Sep 2026) [arXiv]
- Links: [arXiv:2407.14236](https://arxiv.org/abs/2407.14236). Status: 50 pages, "to appear in FCNT"; cited as Front. Comb. Number Theory 4 (2026), 1–44 in [arXiv:2611.16156](https://arxiv.org/pdf/2611.16156).
- Result: the dimension of the Q-span of ζ(3), ζ(5), …, ζ(s−1) is at least (1.119·log s)/(1+log 2) for large even s; the same bound for p-adic zeta values ζ_p(3), …, ζ_p(s−1), refining Sprang 2020. — [arXiv:2407.14236](https://arxiv.org/abs/2407.14236)
- Method: Zudilin's arithmetic factor Φ_n inserted into a slight generalisation of the Ball–Rivoal construction; "The asymptotic impact of the arithmetic factor Φ_n was not anticipated prior to this work." — [arXiv:2407.14236](https://arxiv.org/abs/2407.14236)
- Relation: a refinement of the constant only; the abstract says the archimedean result "is subsumed by a much stronger recent development by Fischler (2026)". — [arXiv:2407.14236](https://arxiv.org/abs/2407.14236)
- Version discrepancy: Lai's January 2025 paper quotes the constant of this work as 1.108/(1+log 2), whereas the v3 abstract gives 1.119. — [arXiv:2501.05321 (PDF)](https://arxiv.org/pdf/2501.05321); [arXiv:2407.14236](https://arxiv.org/abs/2407.14236)

**L. Lai, "A note on the number of irrational odd zeta values, II"** (9 Jan 2025) [arXiv]
- Links: [arXiv:2501.05321](https://arxiv.org/abs/2501.05321). Published: J. Théor. Nombres Bordeaux 38 (2026), no. 2, 469–487 (as cited in [arXiv:2611.16156](https://arxiv.org/pdf/2611.16156)). 17 pages.
- Result (Theorem 1.1): #{odd i ∈ [3, s] : ζ(i) ∉ Q} ≥ 1.284579·√(s/log s) for sufficiently large s. — [arXiv:2501.05321 (PDF)](https://arxiv.org/pdf/2501.05321)
- Method: "combines the elimination technique of Fischler-Sprang-Zudilin (2019) with the Φ_n factor method of Zudilin (2001)". A refinement of the constant from 1.192507 to 1.284579. — [arXiv:2501.05321](https://arxiv.org/abs/2501.05321)

**S. Fischler, L. Lai, J. Sprang, W. Zudilin, "Many more odd zeta values are irrational"** (12 Nov 2026) [arXiv]
- Links: [arXiv:2611.16156](https://arxiv.org/abs/2611.16156), [PDF](https://arxiv.org/pdf/2611.16156). Status: arXiv preprint v1, 17 pages, math.NT, MSC 11J72; 11M06, 33C20. See the next section for the verification details and what the paper says about ζ(5).
- Result (Theorem 1.1): for any ε > 0 and all sufficiently large odd s, N_s ≥ (√(4ζ(2)ζ(3)/ζ(6)) − ε)·√(s·log log s/log s), where N_s counts the irrational values among ζ(3), …, ζ(s) and √(4ζ(2)ζ(3)/ζ(6)) = 3√(70ζ(3))/π² = 2.788258…. — [arXiv:2611.16156 (PDF)](https://arxiv.org/pdf/2611.16156)

### Inferences
- The 2018–2026 results fall into three groups. New constructions: the half-twist and general twists (Zudilin, Sprang), the elimination technique (Fischler–Sprang–Zudilin), and Siegel's lemma (Fischler). Refinements of the construction within a fixed scheme: Lai–Yu (zero design), Rivoal–Zudilin and Lai–Zhou (two irrationals in a short range). Refinements of constants: Lai 2024 and Lai 2025 (Φ_n factors).
- The elimination technique trades linear independence for a count, so the counting record (√(s·log log s/log s)) and the dimension record (√(s/log s) with constant 0.21) are results of different kinds and should not be merged.
- None of these methods isolates a single value: elimination removes a chosen set of "unwanted" zeta values at the cost of enlarging s, so the output is always "at least one among the remaining ones".

### Gaps
- I found no result between 2022 and November 2026 improving "two of ζ(5),…,ζ(35)", and no verified result of the type "at least three of ζ(5),…,ζ(s)" for explicit small s.
- The journal data for Fischler (JLMS 2026), Lai (JTNB 2026) and Lai (FCNT 2026) rest on other authors' reference lists.

## arXiv:2611.16156 (verified) and the post-September-2026 landscape, kept separate from the collection theorems

### Takeaway
arXiv:2611.16156 exists and is relevant: it is "Many more odd zeta values are irrational" by Fischler, Lai, Sprang and Zudilin, submitted 12 November 2026, and it gives the current counting record. It does not prove anything about ζ(5); it borrows a Hankel-determinant device from Fauzan's separate preprint claiming the irrationality of ζ(5), and it declares that its initial proof was produced by an AI system and then verified and rewritten by the authors.

### Cited Findings

**Verification of arXiv:2611.16156** [arXiv]
- Title "Many more odd zeta values are irrational"; authors Stéphane Fischler, Li Lai, Johannes Sprang, Wadim Zudilin; submitted 12 Nov 2026 01:32:10 UTC; 17 pages; no journal reference. — [arXiv:2611.16156](https://arxiv.org/abs/2611.16156)
- Abstract: at least 2.78·√(s·log log s/log s) irrational values among ζ(3), ζ(5), …, ζ(s) for large odd s; "This improves on the best known estimate by a constant factor times √(log log s). The proof combines the elimination technique and an arithmetic lemma of Fischler–Sprang–Zudilin, the hypergeometric constructions of Lai–Yu, and the Hankel determinant method that appeared in the proof that ζ(5) is irrational posted by Fauzan on September 17, 2026." — [arXiv:2611.16156](https://arxiv.org/abs/2611.16156)
- Strategy: instead of one polynomial P_n with L_n(P_n) ∈ Z∖{0} tending to 0, the authors use the Hankel determinant Δ_{n,m} = det(L_n(t^(i+j)))_{0≤i,j<m}. The contradiction combines a denominator bound (§3), an upper bound for |Δ_{n,m}| (§4) and a non-vanishing result for at least one m in a range (§5). — [arXiv:2611.16156 (PDF)](https://arxiv.org/pdf/2611.16156)
- Where the gain comes from: "While using Hankel determinants in [7] gives one access to better arithmetic than by traditional techniques, the improvement in Theorem 1.1 over the previous bound (1.3) comes entirely from the analytic estimates of the determinants." The design of the functional "deviates significantly from its prototype from [7]". — [arXiv:2611.16156 (PDF)](https://arxiv.org/pdf/2611.16156)
- Declared AI use: "The initial proof was completely produced by ChatGPT (GPT-6 Astra, OpenAI). The authors have verified and simplified the argument, rewritten the exposition, and take full responsibility for the mathematical content of this paper." — [arXiv:2611.16156 (PDF)](https://arxiv.org/pdf/2611.16156)
- Stated outlook: the theorem "has the potential to be generalised to a bound on dimensions of the Q-linear spaces generated by odd zeta values … at least this is something that AI tools suggest"; the authors add that understanding the structures "seems to us a more exciting task than targeting particular numbers that still await their irrationality status". — [arXiv:2611.16156 (PDF)](https://arxiv.org/pdf/2611.16156)
- Antecedents named for the determinant method: Zudilin, "A determinantal approach to irrationality" (Constr. Approx. 45 (2017), 301–310) — [arXiv:1507.05697](https://arxiv.org/abs/1507.05697); Brown, "Mellin transforms, transfinite diameter and rational approximations of integrals" (22 Apr 2026, "Submitted to Experimental Maths") — [arXiv:2604.20741](https://arxiv.org/abs/2604.20741); Calegari–Dimitrov–Tang — [arXiv:2408.15403](https://arxiv.org/abs/2408.15403). Also Haynes–Zudilin, "Hankel Determinants of Zeta Values" (SIGMA 11 (2015), 101) — [arXiv:1510.01901](https://arxiv.org/abs/1510.01901).

**Claims about ζ(5) itself (listed separately; none is an "at least one of" theorem)**
- A. Fauzan, "ζ(5) is irrational", Zenodo preprint v1, 17 Sep 2026, DOI 10.5281/zenodo.22826419 [Zenodo]. Stated result and status are given in the first section. — [Zenodo record 22826419](https://zenodo.org/records/22826419)
- How Lai describes the key idea: known linear forms a_{i,j} + b_{i,j}·ζ(5) "do not suffice to prove the irrationality of ζ(5). Fauzan arranged a suitable sequence of such linear forms into a Hankel determinant … assuming that ζ(5) ∈ Q, Fauzan was able to show that v_p(Ĩ_n) is much larger than the entrywise lower bound". — [arXiv:2610.06838 (PDF)](https://arxiv.org/pdf/2610.06838)
- Lai's summary of the earlier state: "The closest results to the irrationality of ζ(5) were 'at least one of ζ(5), ζ(7), ζ(9), ζ(11) is irrational' by Zudilin [60], and '0 < |ζ(5) − p_n/q_n| < q_n^(−0.86) for an explicit rational sequence p_n/q_n' by Brown–Zudilin [9]." — [arXiv:2610.06838 (PDF)](https://arxiv.org/pdf/2610.06838)
- D. Niedbala Giraudin, "A new upper bound for the irrationality exponent of ζ(5)" (1 Nov 2026): exponent at most 171 instead of 260; "The result is conditional on Fauzan's preprint." — [arXiv:2611.01708](https://arxiv.org/abs/2611.01708)
- Brown–Zudilin, "On cellular rational approximations to ζ(5)" (7 Oct 2022, v3 29 Jan 2026): "the traditional strategy for approaching the irrationality of ζ(5) is to try to construct linear forms in 1 and ζ(5) with rational coefficients and good arithmetic properties. At the time of writing, it has not succeeded despite many years of effort." Their cellular integrals on M_{0,8} give 0 < |ζ(5) − p/q| < 1/q^0.86. — [arXiv:2210.03391](https://arxiv.org/pdf/2210.03391)

**Follow-up results using the same determinant method (as listed by Lai on 15 Oct 2026)**
- L. Lai, "The linear independence of 1, ζ(3), and ζ(4)" (5 Oct 2026, v2 15 Oct 2026), 26 pages, arXiv preprint: 1, ζ(3), ζ(4) are linearly independent over Q, hence ζ(3)/π⁴ is irrational. It combines the determinant method with simultaneous approximations of Zudilin and Dauguet–Zudilin; the proof compares Σ_p v_p(J_n)·log p ≥ −(30.83 + o(1))·n² with |J_n| ≤ e^(−31.61·n²). — [arXiv:2610.06838](https://arxiv.org/abs/2610.06838)
- Lai declares that the proof was produced by "ChatGPT (GPT-6 Astra, OpenAI)" from a prompt sent on 1 October 2026, then verified and rewritten by the author. — [arXiv:2610.06838 (PDF)](https://arxiv.org/pdf/2610.06838)
- Lai lists further follow-ups (Kilani on ζ₂(7) and ζ₃(5); Lean 4 formalisations by Liu, Q. Tang and OpenAI) and warns: "In the AI age, many recent results have not appeared on arXiv, which makes it difficult to keep track of all of them." — [arXiv:2610.06838 (PDF)](https://arxiv.org/pdf/2610.06838)

### Inferences
- The November 2026 counting record depends on the Fauzan preprint only for an idea (Hankel determinants of linear functionals), not for its theorem: the abstract and strategy describe a self-contained argument. Its correctness therefore does not logically require ζ(5) to be irrational.
- Both October–November 2026 papers are days-to-weeks-old preprints with AI-generated initial proofs. For a peer-review context they should be labelled "unrefereed preprint, authors take responsibility", not "established".
- The margin in Lai's 1, ζ(3), ζ(4) argument (31.61 against 30.83, about 2.5 %) is again small, as in Zudilin's 2001 theorem.

### Gaps
- I did not verify any proof in this section.
- The items by Liu, Q. Tang and OpenAI are known only through Lai's reference list (GitHub repositories and an "OpenAI Math Release preprint, 24 September 2026"); none was retrieved.
- Whether the Fauzan preprint has been independently refereed or formally verified could not be established from retrievable pages.

## Related constants and analogues, where they clarify what works and what fails for ζ(5)

### Takeaway
The same hypergeometric and elimination machinery transfers to Dirichlet beta values, polylogarithms, q-zeta values and p-adic zeta values, and gives the same kind of output: many irrational values in a family, or one in a short block. Individual values were reached only by different tools: modular or Padé methods for ζ_p(2), ζ_p(3) with p = 2, 3; arithmetic holonomy bounds for L(2, χ₋₃) and ζ₂(5); and, in function fields, t-motivic methods that determine all algebraic relations among Carlitz zeta values.

### Cited Findings

**Catalan's constant and Dirichlet beta values**
- T. Rivoal, W. Zudilin, "Diophantine properties of numbers related to Catalan's constant", Math. Ann. 326:4 (2003), 705–721 [cited]: "at least one of the seven numbers β(2), β(4), …, β(14) is irrational, and … there are infinitely many irrational numbers among the even beta values" — reference and statement as given in [arXiv:1804.09922 (PDF)](https://arxiv.org/pdf/1804.09922).
- W. Zudilin, "Arithmetic of Catalan's constant and its relatives" (26 Apr 2018; Abh. Math. Semin. Univ. Hambg. 89:1 (2019), 45–53): at least one of β(2), β(4), …, β(12) is irrational; Catalan's constant's "irrationality remains an open problem". — [arXiv:1804.09922](https://arxiv.org/abs/1804.09922)
- Lai–Zhou: at least one of β(2), β(4), …, β(10) is irrational. — [arXiv:2103.00904](https://arxiv.org/abs/2103.00904)
- S. Fischler, "Irrationality of values of L-functions of Dirichlet characters" (4 Apr 2019; DOI 10.1112/jlms.12290): generalises the Fischler–Sprang–Zudilin count to Dirichlet L-functions and Hurwitz zeta; "The strategy is different and less elementary: the construction is related to a Padé approximation problem, and a generalization of Shidlovsky's lemma is used to apply Siegel's linear independence criterion." — [arXiv:1904.02402](https://arxiv.org/abs/1904.02402)
- L. Mistiaen, "Linear independence of values of Dirichlet L functions" (2 Dec 2025, in French): dimension of order √(s/log s) for characters mod N with 4 ∤ N. — [arXiv:2512.02612](https://arxiv.org/abs/2512.02612)
- V. Kumar, R. Tosi, "Irrationality of values of L-functions for real Dirichlet characters" (30 Oct 2026): at least c·√(s/log s) irrational values among the first s non-critical values; "The lack of positivity arguments is overcome by applying an irrationality criterion by Lai and Sprang originating from the p-adic setting." — [arXiv:2610.41458](https://arxiv.org/abs/2610.41458)
- L. Lai, J. Li, "A partial result towards the Chowla–Milnor conjecture" (v3 22 Oct 2026, "to appear in J. Lond. Math. Soc"): a log q dimension bound for Hurwitz zeta differences, by a proof "similar to those of Ball–Rivoal (2001) and Zudilin (2002)". — [arXiv:2505.12687](https://arxiv.org/abs/2505.12687)
- Auxiliary: Krattenthaler–Rivoal, "On a linear form for Catalan's constant" — [arXiv:0810.1927](https://arxiv.org/abs/0810.1927); Zudilin, "An Apéry-like difference equation for Catalan's constant" — [arXiv:math/0201024](https://arxiv.org/abs/math/0201024).
- 2026 claims about Catalan's constant itself, listed for caution: Z.-W. Sun, "Catalan's constant is irrational" (3 Sep 2026) — [arXiv:2609.04176](https://arxiv.org/abs/2609.04176); D. Wachs's note concludes that "the published derivation does not justify the claimed estimate, so the proof of the main theorem is incomplete. Nothing is concluded about the irrationality of G." — [arXiv:2609.22339](https://arxiv.org/abs/2609.22339). Lai separately cites a Lean 4 formalised proof by OpenAI of the irrationality of Catalan's constant (not retrieved). — [arXiv:2610.06838 (PDF)](https://arxiv.org/pdf/2610.06838)

**Polylogarithms and G-functions**
- For nonzero algebraic z with |z| < 1, dim over Q(z) of Span(1, Li₁(z), …, Li_s(z)) ≥ ((1−ε)/((1+log 2)·[Q(z):Q]))·log s, "due to Rivoal [26] for z ∈ R, to Marcovecchio [21] in the general case" (Rivoal, J. Théor. Nombres Bordeaux 15 (2003), 551–559; Marcovecchio, Ann. Sc. Norm. Sup. Pisa V (2006), 1–11) [cited]. Fischler improves this to √(s/log s). — [arXiv:2109.10136 (PDF)](https://arxiv.org/pdf/2109.10136)
- Fischler–Rivoal, "Linear independence of values of G-functions" (31 Jan 2017): a log S lower bound for spaces generated by values of F_n^[s], encompassing the polylogarithm case. — [arXiv:1701.09051](https://arxiv.org/abs/1701.09051). Part II, outside the disk of convergence (Ann. Math. Québec 45 (2021), 53–93) — [arXiv:1811.08758](https://arxiv.org/abs/1811.08758).
- L. Mistiaen, polylogarithms with periodic coefficients outside the disk of convergence (7 Oct 2026). — [arXiv:2610.09774](https://arxiv.org/abs/2610.09774)

**q-analogues**
- C. Krattenthaler, T. Rivoal, W. Zudilin, "Séries hypergéométriques basiques, q-analogues des valeurs de la fonction zêta et séries d'Eisenstein" (J. Inst. Math. Jussieu 5 (2006), 53–79): arithmetic of ζ_q(s) = Σ_k q^k Σ_{d|k} d^(s−1) for 1/q an integer. — [arXiv:math/0311033](https://arxiv.org/abs/math/0311033)
- F. Jouhet, E. Mosaki, "Irrationalité aux entiers impairs positifs d'un q-analogue de la fonction zêta de Riemann" (11 Dec 2007, 33 pages; no journal reference in the arXiv record): a new lower bound for the dimension of the span of 1, ζ_q(3), ζ_q(5), …, ζ_q(A−1), improving Krattenthaler–Rivoal–Zudilin; "for 1/q ∈ Z∖{−1; 1}, at least one of the numbers ζ_q(3), ζ_q(5), ζ_q(7), ζ_q(9) is irrational." — [arXiv:0712.1762](https://arxiv.org/abs/0712.1762)
- Jouhet–Mosaki, q-analogues of Dirichlet's beta function (26 Nov 2008). — [arXiv:0811.4287](https://arxiv.org/abs/0811.4287)
- Zudilin, "Diophantine problems for q-zeta values" (Math. Notes 72 (2002), 858–862). — [arXiv:math/0206179](https://arxiv.org/abs/math/0206179)

**p-adic zeta values**
- F. Calegari, "Irrationality of certain p-adic periods for small p" (16 Aug 2004): the p-adic analogue of ζ(3) is irrational for p = 2 and 3, and so is the 2-adic analogue of Catalan's constant, via overconvergent p-adic modular forms. — [arXiv:math/0408214](https://arxiv.org/abs/math/0408214)
- F. Beukers, "Irrationality of some p-adic L-values" (12 Mar 2006): irrationality of ζ_p(k) for p = 2, 3 and k = 2, 3 "using classical continued fractions discovered by Stieltjes". — [arXiv:math/0603277](https://arxiv.org/abs/math/0603277)
- J. Sprang, "A Linear independence result for p-adic L-values" (Duke Math. J. 169 (2020), 3439–3476): a p-adic Ball–Rivoal theorem, dim_K(K + Σ_{i=2}^{s+1} L_p(i, χω^(1−i))·K) ≥ (1−ε)·log s/(2[K:Q](1+log 2)). — [arXiv:1809.07714](https://arxiv.org/abs/1809.07714)
- L. Lai, "On the irrationality of certain 2-adic zeta values" (3 Apr 2023): for every s ≥ 0 there is an odd j ∈ [s+3, 3s+5] with ζ₂(j) irrational; in particular one of ζ₂(7), ζ₂(9), ζ₂(11), ζ₂(13). "The most difficult step is proving that certain Volkenborn integrals are nonzero." — [arXiv:2304.00816](https://arxiv.org/abs/2304.00816)
- L. Lai, J. Sprang, "Many p-adic odd zeta values are irrational" (17 Jun 2023; "to appear in Michigan Math. J"): at least (c_p − ε)·√(s/log s) irrational values among ζ_p(3), …, ζ_p(s). "The main difficulty consists in proving the non-vanishing of the resulting linear forms. We overcome this problem by using a new irrationality criterion." — [arXiv:2306.10393](https://arxiv.org/abs/2306.10393)
- L. Lai, C. Lupu, J. Sprang, "On the irrationality of certain p-adic zeta values" (29 May 2025): "the p-adic analogue of Zudilin's theorem"; for any prime p ≥ 5 there is an odd i ∈ [3, p + p/log p + 5] with ζ_p(i) irrational. — [arXiv:2505.23088](https://arxiv.org/abs/2505.23088)
- F. Calegari, V. Dimitrov, Y. Tang, "The linear independence of 1, ζ(2), and L(2, χ₋₃)" (27 Aug 2024, 218 pages, "the submitted version"): irrationality of L(2, χ₋₃) by "a new kind of arithmetic holonomy bound" applied to a construction of Zagier. — [arXiv:2408.15403](https://arxiv.org/abs/2408.15403)
- L. Lai, J. Sprang, W. Zudilin, "A note on the irrationality of ζ₂(5)" (IMRN 2026, Issue 16, rnag180): Apéry-style approximations give a new proof of the irrationality of ζ₂(5), "the result established recently by Calegari, Dimitrov and Tang using a different method", with μ(ζ₂(5)) ≤ 20.342…. — [arXiv:2505.05005](https://arxiv.org/abs/2505.05005)
- M. Kawashima, A. Poëls, a linear independence criterion for p-adic Hurwitz zeta values at odd integers (11 Oct 2026), generalising Beukers. — [arXiv:2610.15525](https://arxiv.org/abs/2610.15525)
- I. Kilani, Zenodo preprints claiming ζ₂(7) irrational (30 Sep 2026) and ζ₃(5) irrational (7 Oct 2026), both by Hankel determinants of Volkenborn-integral moments [Zenodo]. — [Zenodo 23058762](https://zenodo.org/records/23058762); [Zenodo 23203596](https://zenodo.org/records/23203596)

**Function fields**
- C.-Y. Chang, M. Papanikolas, J. Yu, "Frobenius difference equations and algebraic independence of zeta values in positive equal characteristic" (Algebra & Number Theory 5 (2011), 111–129): among Carlitz zeta values ζ_r(n) for all r and n, all algebraic relations are those within each individual family. — [arXiv:0804.0038](https://arxiv.org/abs/0804.0038)
- Chang–Papanikolas–Yu, "Geometric Gamma values and zeta values in positive characteristic" (IMRN 2010): all algebraic relations among these special values come from the standard functional equations, the Euler–Carlitz relations and the Frobenius p-th power relations. — [arXiv:0905.2876](https://arxiv.org/abs/0905.2876)

### Inferences
- In every setting, well-poised hypergeometric constructions give collection-type theorems of the same shape as for odd zeta values (one of four q-zeta values, one of five beta values, one in a p-adic window). The transfer of the limitation is as systematic as the transfer of the method.
- Where a single value was reached before 2026 (ζ₂(3), ζ₃(3), L(2, χ₋₃), ζ₂(5)), the proofs used extra structure: modularity, classical continued fractions, or holonomy bounds. This supports the reading that ζ(5) required a new arithmetic input and not a better choice of hypergeometric parameters.
- The function-field theorems show the expected final answer (no relations beyond the obvious ones), but their tools (Frobenius difference equations, t-motives) have no known counterpart for ζ(2n+1) over Q.

### Gaps
- Rivoal–Zudilin (Math. Ann. 2003), Rivoal (polylogarithms, 2003) and Marcovecchio (2006) were not retrieved.
- I did not retrieve the older function-field literature (for example the transcendence of all Carlitz zeta values); only the two Chang–Papanikolas–Yu records above were verified.
- A journal reference for Calegari's 2004 preprint was not confirmed (the arXiv record says "Preprint").

## Surveys, lecture material and recorded talks

### Takeaway
The best verified entry points are Fischler's Séminaire Bourbaki exposé (in French, with full statements and proof sketches up to 2002), the Krattenthaler–Rivoal memoir for the arithmetic of the linear forms, Zudilin's D.Sc. monograph (in Russian), and the introductions of the 2018–2026 papers for the modern state. I found talk announcements but no recording of a talk on odd zeta values by Zudilin, Fischler, Sprang or Lai.

### Cited Findings

**Verified texts** [arXiv]
- S. Fischler, "Irrationalité de valeurs de zêta (d'après Apéry, Rivoal, …)" (5 Mar 2003): Séminaire Bourbaki 2002–2003, exposé no. 910 (Nov. 2002); Astérisque 294 (2004), 27–62; 45 pages, in French. Part 1: all known proofs of Apéry's theorem and their connections. Part 2: proof of the Ball–Rivoal theorem. Part 3: quantitative statements (ζ(5),…,ζ(21); ζ(5),…,ζ(11)). — [arXiv:math/0303066](https://arxiv.org/abs/math/0303066)
- W. Zudilin, "Apéry's theorem and problems for the values of Riemann's zeta function and their q-analogues" (25 Dec 2013): 118 pages, in Russian; described as the author's habilitation (D.Sc.) thesis; Chapter 1 is about the four-number theorem. — [arXiv:1312.6919](https://arxiv.org/abs/1312.6919)
- C. Krattenthaler, T. Rivoal, Mem. Amer. Math. Soc. 186 (2007): 73 pages, in French; the reference for denominators and for very-well-poised identities. — [arXiv:math/0311114](https://arxiv.org/abs/math/0311114)
- F. Brown, "Irrationality proofs for zeta values, moduli spaces and dinner parties" (19 Dec 2014; Mosc. J. Comb. Number Theory 6 (2016), as cited by Lai): a geometric framework on M_{0,n} that "contains a generalisation of the linear forms used by Ball and Rivoal". — [arXiv:1412.6508](https://arxiv.org/abs/1412.6508)
- W. Zudilin, "Hypergeometric heritage of W.N. Bailey" (2016; ICCM Notices 7:2 (2019)): reviews Bailey's hypergeometric work with applications in number theory. — [arXiv:1611.08806](https://arxiv.org/abs/1611.08806)
- Useful prerequisites named in the sources: the prime number theorem in the form d_n = lcm(1,…,n) = e^(n+o(n)) and Stirling's formula (elementary proofs) — [arXiv:1801.09895 (PDF)](https://arxiv.org/pdf/1801.09895); the saddle-point method (Rivoal 2002, Zudilin 2001–2004, Rivoal–Zudilin 2018) — [arXiv:1803.03160 (PDF)](https://arxiv.org/pdf/1803.03160); Siegel's lemma and Shidlovsky's lemma (Fischler) — [arXiv:2109.10136](https://arxiv.org/abs/2109.10136).

**Found by web search only; pages not opened** [cited]
- T. Rivoal, "Séries hypergéométriques et irrationalité des valeurs de la fonction zêta de Riemann", J. Théor. Nombres Bordeaux 15 (2003), no. 1, 351–365 (Numdam item JTNB_2003__15_1_351_0; EuDML doc 249077).
- T. Rivoal, doctoral thesis "Propriétés diophantiennes de la fonction zêta de Riemann aux entiers impairs", Université de Caen, 2001 (HAL tel-00004519, per a search summary of Wikipedia). I found no habilitation memoir by Rivoal.
- Notes by S. Fischler of a lecture by T. Rivoal at Institut Henri Poincaré on the Ball–Rivoal proof, including Nesterenko's criterion (file inde_fischlerrivoal.pdf on Fischler's Orsay page).
- M. Waldschmidt, "Recent Diophantine results on zeta values: a survey" (RIMS 2009, PDF on his IMJ-PRG page), and 2009 lecture slides covering the criteria of Siegel, Nesterenko and Fischler.
- Talk announcements without recordings: W. Zudilin, "Irrationality of odd zeta values and other treasures", MPIM Bonn, 2018 (two search summaries disagree on the month, May or September); W. Zudilin, UNSW number theory seminar, 2018; W. Zudilin, "Many odd zeta values are irrational", Radboud Geometry Seminar, 3 Oct 2018; Li Lai, "On the irrationality of certain 2-adic zeta values", PKU/BICMR number theory seminar, 5 Apr 2023 (researchseminars.org/talk/PekiNT/2).
- Recorded talks on the neighbouring Calegari–Dimitrov–Tang work: Y. Tang, "The arithmetic of power series and applications to irrationality", conference in memory of Jan Nekovář, 9–13 Oct 2023 (carmin.tv); V. Dimitrov, "On an arithmetic criterion for holonomicity", BIRS (UBC cIRcle, 58 minutes); F. Calegari, "Searching for Sequences: Irrationality Beyond Apéry", DIMACS, Sept 2024 (video on Vimeo per the listing).
- Journalism: Quanta Magazine, "Rational or Not? This Basic Math Question Took Decades to Answer", 8 Jan 2025 (on Calegari–Dimitrov–Tang; quotes Zudilin, per a search summary).

### Inferences
- For a reader who needs exact statements with proofs in one place, the Bourbaki exposé is the only verified survey; it stops in 2002 and must be supplemented by the introductions of arXiv:1803.08905, arXiv:2109.10136 and arXiv:2611.16156.
- There is no verified survey covering 2018–2026. The introduction of arXiv:2610.06838 (October 2026) is the most recent short overview, but it is a research preprint, not a survey.

### Gaps
- No YouTube recording of a talk by Zudilin, Fischler, Sprang or Lai on odd zeta values was found in three searches. Recordings may exist on institutional pages that could not be opened.
- Zudilin's survey "Arithmetic hypergeometric series" (Russian Math. Surveys 66 (2011), no. 2, 369–420) and "Apéry's theorem. Thirty years after" (Int. J. Math. Comput. Sci. 4 (2009), 9–19) are known here only from Lai's reference list — [arXiv:2610.06838 (PDF)](https://arxiv.org/pdf/2610.06838).
- Fischler's and Zudilin's home pages (publication lists, slides, course notes) could not be opened.

## Items not retrieved, and access limits of this search

### Takeaway
Only arxiv.org (pages, PDFs, API) and zenodo.org were directly reachable in this session, plus web-search result snippets. Publisher and institutional pages were blocked, so every journal reference for a paper without an arXiv journal field is second-hand, and several items requested in the brief could only be confirmed bibliographically.

### Cited Findings
- Hosts that failed to resolve or were refused in this session: link.springer.com, www.numdam.org, www.mathnet.ru, www.imo.universite-paris-saclay.fr, www.math.ru.nl, en.wikipedia.org, www.mpim-bonn.mpg.de, jtnb.centre-mersenne.org, comptes-rendus.academie-sciences.fr, journals.impan.pl; londmathsoc.onlinelibrary.wiley.com returned HTTP 403 with a challenge page. These would have to be allowed in the session's network settings for direct verification. (Session observation; no external source.)
- Confirmed bibliographically only (existence, title, authors and journal confirmed through reference lists of retrieved papers or search snippets; contents not read):
  - Ball–Rivoal, Invent. Math. 146 (2001), 193–207 — cited in [arXiv:2611.16156](https://arxiv.org/pdf/2611.16156).
  - Zudilin, Russian Math. Surveys 56 (2001), no. 4, 774–776 — cited in [arXiv:2610.06838](https://arxiv.org/pdf/2610.06838).
  - Zudilin, Izv. Math. 66 (2002), no. 3, 489–542 — cited in [arXiv:2610.06838](https://arxiv.org/pdf/2610.06838); only the short conference version [arXiv:math/0104249](https://arxiv.org/abs/math/0104249) was read.
  - Zudilin, "Well-poised hypergeometric service for diophantine problems of zeta values", J. Théor. Nombres Bordeaux 15 (2003), 593–626 — search snippet only.
  - Fischler–Zudilin, Math. Ann. 347 (2010), no. 4, 739–763 — cited in [arXiv:2610.06838](https://arxiv.org/pdf/2610.06838).
  - Nesterenko, Vestnik Moskov. Univ. 1985 — cited in [arXiv:2109.10136](https://arxiv.org/pdf/2109.10136).
  - Fischler–Sprang–Zudilin, C. R. Math. Acad. Sci. Paris 356 (2018), no. 7, 707–711 — cited in [arXiv:2611.16156](https://arxiv.org/pdf/2611.16156).
  - Rivoal–Zudilin, Math. Ann. 326:4 (2003), 705–721 — cited in [arXiv:1804.09922](https://arxiv.org/pdf/1804.09922).
  - Fischler–Rivoal, "Approximants de Padé et séries hypergéométriques équilibrées", J. Math. Pures Appl. 82 (2003), no. 10, 1369–1394 — cited in [arXiv:2109.10136](https://arxiv.org/pdf/2109.10136).
  - Zudilin, "A third-order Apéry-like recursion for ζ(5)", Math. Notes 72 (2002), no. 5, 733–737 — cited in [arXiv:2610.06838](https://arxiv.org/pdf/2610.06838).
- Retrieved but judged unreliable, and excluded from the record table:
  - S. Suman, [arXiv:2609.22316](https://arxiv.org/abs/2609.22316) (three-number claim; see the 2001–2012 section).
  - M. Kabongo, "La Zeta de Riemann est irrationnelle aux impairs positifs" (4 Oct 2018), which claims ζ(j) irrational for every j ≥ 3 from a two-dimensional linear system; no journal reference, and no later paper in this bibliography cites it. — [arXiv:1810.02361](https://arxiv.org/abs/1810.02361)

### Inferences
- The core chain of results is nevertheless well anchored: every theorem in the record table is stated in at least one retrieved arXiv abstract or PDF, and the pre-2018 statements are cross-checked against the Bourbaki exposé.
- Journal volume and page data for 2026 publications should be re-checked on the publishers' pages before being printed in a final report.

### Gaps
- Not retrieved; one line each, unverified (from memory):
  - W. Zudilin, "One of the eight numbers ζ(5), ζ(7), …, ζ(17), ζ(19) is irrational", Math. Notes 70 (2001) — unverified (from memory); the eight-number statement itself is verified through Krattenthaler–Rivoal.
  - T. Rivoal, habilitation memoir (Grenoble) — unverified (from memory); not found by search.
  - J. Yu, "Transcendence and special zeta values in characteristic p", Ann. of Math. (1991) — unverified (from memory).
  - Journal version of Calegari's p-adic periods paper (Int. Math. Res. Not., 2005) — unverified (from memory).
- I could not determine whether any 2023–2026 work other than those listed improves the small-range statements; the arXiv searches used (author searches for Lai, Sprang, Fischler, Zudilin, Rivoal, and keyword searches on titles and abstracts for July–November 2026) returned nothing further.

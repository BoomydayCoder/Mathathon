# R1_lit_perfect -- Voronoi theory, random perfect / extreme lattices, small-dimension data

Agent wall clock: 07:00 -> ~07:12 UTC (2026-11-25). Normalisation used below: V(L) = 2^n * (packing density of L).
Status tags: [PAPER] = stated in the paper, read in full text here; [ABS] = abstract only; [MEM] = from memory, not verified;
[MINE] = my own re-computation / inference from the paper's numbers (not a claim of the paper).

## 0. Access log (IDs in the task statement were partly wrong)

| paper | correct arXiv id | local copy (BASE/papers/...) | access |
|---|---|---|---|
| Andreanov, Scardicchio, "Random perfect lattices and the sphere packing problem", PRE 86, 041117 (2012) | **1202.5673** (not 1202.3275) | andreanov_scardicchio_1202.5673/ | full text (pdftotext), tables + conclusions read |
| Kallus, "Statistical mechanics of the lattice sphere packing problem", PRE 87, 063307 (2013) | **1305.1310** (not 1305.1301) | kallus_statmech_1305.1310/ | full text skimmed, Table I read |
| Kallus, Marcotte, Torquato, "Jammed lattice sphere packings", PRE 88, 062151 (2013) | **1310.5971** (not 1304.6928) | kallus_marcotte_torquato_1310.5971/ | full text, density section + Table I read |
| Andreanov, Scardicchio, Torquato, "Extreme lattices: symmetries and decorrelations" (2016) | 1309.1301 | andreanov_scardicchio_torquato_1309.1301/ | full text skimmed (Table I, Sec. III, conclusions) |
| Marcotte, Torquato, "Efficient LP algorithm to generate the densest lattice sphere packings" (2013) | 1304.5003 | marcotte_torquato_1304.5003/ | grep only |
| Dutour Sikiric, van Woerden, "The lattice packing problem in dimension 9 by Voronoi's algorithm" (v2 Feb 2026) | 2508.20719 | dutour_vanwoerden_dim9_2508.20719/ | full text, Tables 1,3,7,8 read |
| Bacher, "On the number of perfect lattices" | 1704.02234 | bacher_number_perfect_1704.02234/ | abstract/theorem statement only |
| van Woerden, "An upper bound on the number of perfect quadratic forms" (Adv. Math. 2020) | 1901.04807 | vanwoerden_upper_bound_1901.04807/ | abstract only |
| Vanoni, Guo, Torquato, "Analysis of conjectural improvements to Minkowski's lower bound..." (Sep 2026) | 2608.00920 | conjectural_minkowski_2608.00920/ | abstract only |
| Menabe, "A combinatorial approach to Voronoi theory" (Aug 2026) | 2608.06897 | combinatorial_voronoi_2608.06897/ | abstract only (new proof of Voronoi finiteness; nothing on density) |
| Cohn, Elkies, "New upper bounds on sphere packings I" | math/0110009 | cohn_elkies_math0110009/ | Table 3 only (used to cross-check the CSV) |

NOT reachable from this session (DNS / proxy allow-list): www.math.rwth-aachen.de (Nebe-Sloane catalogue) and cohn.mit.edu.
Hence the CSV is built from memory + cross-checks (see Section 6). A mis-download (Kallus, "Pessimal packing shapes", 1305.0289)
sits in papers/kallus_pessimal_shapes_1305.0289_WRONG_ID/ and is irrelevant.

## 1. Andreanov-Scardicchio 2012 (random walk on the Voronoi graph) [PAPER]

* Method: randomised Voronoi algorithm. At a perfect form Q, pick a random extreme ray of the cone of facets of the
  Ryshkov polyhedron at Q (= a random Voronoi neighbour) and hop there ("infinite temperature"); a Metropolis variant with
  energy e = -(1/d) ln(phi) and inverse temperature beta biases towards dense forms. Start: A_d (d<=12) or D_d (larger d).
* Dimensions d = 8..19. Sample sizes 10^6 (d=8-12), 10^5 (d=13-16), 2*10^5 (d=17,18), 1.5*10^5 (d=19). Isometry test
  only for d < 13. 5*10^4 lattices in d=19 took 5-7 days. They say d <= 40 looks technically feasible.
* Table II, mean energy <e> of the walk, sigma_e, and best known (natural log):

| d | <e> | sigma_e | e(best known) | [MINE] V_typ = 2^d exp(-d<e>) | V_typ/d^2 | V_typ/d^3 | V_best |
|---|---|---|---|---|---|---|---|
| 8 | 0.180572 | 0.0215 | 0.171465 | 60.4 | 0.943 | 0.118 | 64.9 |
| 9 | 0.23352 | 0.0182 | 0.21396 | 62.6 | 0.773 | 0.086 | 74.6 |
| 10 | 0.26828 | 0.0142 | 0.23857 | 70.0 | 0.700 | 0.070 | 94.2 |
| 11 | 0.29347 | 0.0123 | 0.25511 | 81.2 | 0.671 | 0.061 | 123.8 |
| 12 | 0.31505 | 0.0080 | 0.25041 | 93.4 | 0.649 | 0.054 | 202.9 |
| 13 | 0.33106 | 0.0033 | 0.27178 | 110.7 | 0.655 | 0.050 | 239.3 |
| 14 | 0.34277 | 0.0024 | 0.27386 | 135.0 | 0.689 | 0.049 | 354.3 |
| 15 | 0.35405 | 0.0027 | 0.27218 | 161.8 | 0.719 | 0.048 | 552.5 |
| 16 | 0.36507 | 0.0020 | 0.26370 | 190.4 | 0.744 | 0.046 | 964.1 |
| 17 | 0.37205 | 0.0028 | 0.27833 | 234.8 | 0.812 | 0.048 | 1155.0 |
| 18 | 0.38322 | 0.0024 | 0.28489 | 264.7 | 0.817 | 0.045 | 1554.1 |
| 19 | 0.39000 | 0.0039 | 0.28903 | 317.3 | 0.879 | 0.046 | 2160.8 |

* Their two fits of <e>(d), declared "equally good", resolvable "only for d >> 40":
  (37) <e> = (0.58 +- 0.04) - (0.9 +- 1.0) log(d)/d - (0.8 +- 0.6)/d  ==> phi ~ 2^{-(0.84 +- 0.06) d}  (0.58/ln2 = 0.84);
  (38) <e> = (0.066 +- 0.04) log d + (0.27 +- 0.04) - (1.4 +- 0.2)/d  ==> phi ~ d^{-(0.06 +- 0.04) d}.
  So the exact statement is: "two hypotheses between which we cannot distinguish" -- NOT a claim that typical perfect
  lattices beat Minkowski exponentially.
* [MINE] How solid: weak. 12 points, 3-parameter fits, 100% error on the log(d)/d coefficient. In V-normalisation the same
  data satisfy 0.65 <= V_typ/d^2 <= 0.94 for ALL d = 8..19, i.e. they are equally well described by the Minkowski exponent
  with a polynomial correction, phi_typ ~ c d^nu 2^{-d}. Log-log fit of V_typ = C d^nu: nu = 2.03 (d=8..19), 2.68 (d=12..19),
  2.79 (d=14..19). The apparent "0.84 < 1" exponent is exactly what a d^2..d^3 prefactor looks like at d <= 19.
* Kissing numbers: typical perfect lattice has Z of order d^2; <Z>/(d(d+1)) lies between 1.0 and ~1.5 for d=8..19 (Fig. 16),
  with zero inter-quartile range at d = 18, 19 (i.e. essentially all sampled forms have exactly Z = d(d+1)).
* Fraction of sampled perfect forms that are eutactic (= extreme), Table V: 0.997 (8), 0.830 (9), 0.738 (10), 0.479 (11),
  0.134 (12), 5.1e-3 (13), 3.0e-4 (14), 1.3e-4 (15), 8e-5 (16), 6e-5 (17), 1.25e-5 (18), 2e-5 (19).
* Random walk finds the densest lattice ~80% of the time in d=8, <1% in d=10; plain walk recovers best known for d<=12,
  temperature bias needed for 13..18, many restarts for d=19 ("glassy landscape"). Degree distribution of the Voronoi
  graph conjectured scale-free, exponent 2.6 +- 0.1.
* [MINE] CAVEAT on "typical": the walk's stationary measure is degree-biased, not uniform on perfect forms. Proof that the
  bias is huge: exact d=9 fraction of extreme forms is 7 338 582 / 2 237 251 040 = 0.33% (Section 5) whereas the walk
  reports 0.83 (AS) / 0.14 (AST Table I). So AS "typical" densities are biased UP relative to a uniformly random vertex.

## 2. Andreanov-Scardicchio-Torquato 2016 (extreme lattices) [PAPER]

* Same stochastic Voronoi algorithm (weak Metropolis bias, walks of 10^6-10^7 steps); representative samples of extreme
  lattices only for d = 8..13 (2*10^4 extreme lattices from 6*10^6 steps in d=13); for d = 14..19 only 200 - few thousand
  extreme lattices, from biased walks, explicitly "not representative".
* Findings: "extreme lattices are typically denser than perfect lattices" (Fig. 4); bulk dominated by typical lattices with
  similar density and kissing number; best and worst packers in long tails; A_d is consistently the LEAST dense extreme
  lattice found (Coxeter's conjecture, not disproved). Typical extreme lattices have far fewer symmetries and decorrelate
  faster (pair correlation closer to Poisson) than the densest ones; d=13 is an exception where the densest is "typical-like".
* No asymptotic density formula is claimed; the density scaling of typical extreme lattices is posed as an open problem
  ("Do typical extreme lattices provide (further) improvement of the lower bound?").

## 3. Kallus 2013 and Kallus-Marcotte-Torquato 2013 [PAPER]

* Kallus: Monte Carlo on the space of lattices in an isobaric ensemble (reduced pressure p plays the role of inverse
  temperature), quasi-static compression; de novo recovery of the densest known lattices for d = 9..20 with success rates
  1, 1, 0.55, 1, 0.70, 0.60, 0.90, 0.95, 0.8, 0.6, 0.3, 0.14 (d = 9..20; 20 runs each, 50 in d=20). Interpretation: evidence
  of a first-order "crystallisation" transition as d -> infinity; no glass transition seen up to d=20 but left open.
  No scaling law for typical densities. Table I confirms 1/delta for the best lattices in d = 9..20 (used for the CSV).
* KMT: ensemble = final states of the Marcotte-Torquato sequential-LP compression from random initial lattices, keeping only
  ISOSTATIC extreme lattices (exactly d(d+1) contacts = minimum possible for a perfect form). Runs / isostatic kept:
  d=13: 10000/365; 14: 10000/1625; 15: 10000/5196; 16: 10000/6761; 17: 10000/9235; 18: 10000/9590; 19: 20000/19200;
  20: 20000/19085; 21: 10000/9473; 22: 10000/9406; 23: 10000/9281; 24: 10000/9205. "We suspect that this fraction should
  theoretically tend to one" (losses attributed to numerical error). All kept lattices in d = 15..24 have distinct determinants.
* **Density law (numerical, d = 15..24): 2^d <phi> = 0.0567 d^{3.01}, nu = 3.01 +- 0.01**; relative std 3.3% (d=15) -> 2.9%
  (d=24). In our normalisation: V_typ(jammed isostatic lattice) ~ 0.0567 n^3, e.g. 197 (n=15), 400 (n=19), 809 (n=24).
  They contrast this with nu in [1,2] for classical jammed/RSA hard spheres and note that the naive heuristic phi ~ z 2^{-d}
  would give only d^2 2^{-d}.
* Other statistics: near-contact pair correlation g(r) ~ (r-1)^{-gamma}, gamma = 0.307-0.320; number of near contacts
  <Z(1+xi)> - d(d+1) ~ A xi^{1-gamma} with A = 1.43e-3 d^{3.30}; contact-force distribution P(f) ~ f^theta, theta = 0.35-0.38.
* [MINE] KMT/AS ratio at equal d: 1.22 (d=15), 1.26 (d=19): extreme isostatic lattices are ~25% denser than walk-typical
  perfect ones, same polynomial order. Extrapolating the KMT law is NOT justified by the paper (10 dimensions, one decade
  of nothing); the exponent 3.01 +- 0.01 is a fit error, not a systematic one.

## 4. What this means for the n^2 question [MINE, heuristic]

* Klartag's process freezes exactly at a perfect form (his eq. (44): the x (x) x for contact points span Sym_n, so
  |bd E cap L| >= n(n+1)); his paper never uses the words perfect / Voronoi / Ryshkov (grep), and his volume bound is
  for a fixed time T with positive probability, not for the frozen form.
* Perfection alone gives nothing: A_n is perfect and extreme with phi ~ n^{-n/2} (e(A_d) ~ (1/2) ln d, AS eqs. 31-32),
  i.e. V(A_n) -> 0 super-exponentially. Any bound must use typicality / the path, not the endpoint property.
* Empirically (n <= 24) generic local optima have V ~ 0.7-0.9 n^2 (walk-typical perfect) or ~0.057 n^3 (SLP-jammed
  isostatic extreme). Both are the same order as, or one power of n above, Klartag's c n^2; nothing in the data suggests
  generic perfect/extreme forms do better than a fixed power n^2..n^3, and the gap to the best known lattices widens
  quickly with n (V_best/V_typ: 1.07 at n=8, 2.2 at n=12, 6.8 at n=19; Leech 32373 vs KMT-typical 809 at n=24, factor 40).
* Counting heuristic for the exponent: an ellipsoid of normalised volume V has ~ V n eps lattice points in a shell of
  relative width eps (Poisson heuristic). n(n+1) contacts with V ~ n^2 corresponds to eps ~ 1/n; V ~ n^3 to eps ~ 1/n^2.
  So "nu = 3" would mean a local optimiser squeezes the contact shell by an extra factor n compared with Klartag's
  accounting ("a constant number of contacts per unit of swept volume").

## 5. Number of perfect forms; dimension 9 statistics

* Counts of non-similar perfect forms [PAPER: DSvW Table 1]: d=2: 1 (Lagrange); 3: 1 (Gauss); 4: 2, 5: 3 (Korkine-Zolotarev);
  6: 7 (Barnes 1957); 7: 33 (Jaquet-Chiffelle 1993); 8: 10 916 (Dutour Sikiric-Schuermann-Vallentin 2005/07);
  9: >= 500 000 (DSV), >= 23 000 000 (van Woerden 2018), **2 237 251 040 (complete, Dutour Sikiric-van Woerden 2025/26)**.
  Corollaries there: Lambda_9 is the unique densest 9-dim lattice, gamma_9 = 2; 7 338 582 extreme forms; possible kissing
  numbers of 9-dim lattices are exactly 2*{1..91, 99, 120..129, 136}. ~3 million core hours.
* Asymptotics: Bacher [ABS] e^{d^{1-eps}} < p_d < e^{d^{3+eps}} eventually; van Woerden [ABS] improves the upper bound to
  e^{O(d^2 log d)} (volumetric argument, Voronoi's first reduction); DSvW quote the pair as exp(Omega(d)) < p_d <
  exp(O(d^2 log d)) and remark the true growth looks closer to the upper bound. (Soule's older e^{O(d^4 log d)}: secondary
  source only.) I found NO result on the typical/expected determinant or Hermite invariant of perfect forms.
* Dimension 9, perfect forms by half kissing number s = |Min|/2 [PAPER Table 3/8; percentages MINE]:
  s=45 (minimum, "isostatic"): 1 353 947 672 = 60.5%; 46: 21.1%; 47: 12.0%; 48: 3.78%; 49: 1.67%; 50: 0.60%; then roughly
  exponential decay; maximum 136 (Lambda_9). "Just as in dimension 8 the minimal kissing number is attained by more than
  half of the perfect forms."
* Extreme forms in d=9 by s: 45: 986 573 (13.4%); 46: 1 494 396 (20.4%); 47: 1 597 977 (21.8%, the mode); 48: 1 264 705
  (17.2%); 49: 861 126 (11.7%). Probability that a perfect form with s pairs is extreme: 0.073% (s=45), 0.32% (46), 0.60%
  (47), 1.5% (48), 2.3% (49), 3.9% (50), 7.7% (52), 10.9% (54). So in d=9 extreme forms are NOT typically isostatic,
  whereas KMT find > 92-96% isostatic among SLP end states for d >= 17 (different ensemble, higher dimension).
* 96.7% of the 9-dim perfect forms have trivial automorphism group {+-1}; > 99.9% have |Aut| <= 4.

## 6. Rigorous links between Voronoi-graph random walks and density lower bounds

None found other than Klartag 2025. The physics papers above are purely numerical. Vanoni-Guo-Torquato 2026 [ABS] is
conditional on the Torquato-Stillinger realizability conjecture (pair-correlation LP; gives phi >~ d^beta 2^{-d} for every
fixed beta > 1 and exponential rates 2^{-0.7786 d}, 2^{-0.6226 d} conditionally) and has no lattice / Voronoi content.

## 7. best_known_lattices.csv (53 rows: n = 1..48, 56, 64, 72, 80, 128)

Columns: n, center_density, log2_center_density, name, kissing_number, source, packing_density, V_eq_2^n_density,
V_over_n, V_over_n2, V_over_n3. Generated by gen_csv.py (same directory). Reliability:
* n = 1..24: exact closed forms (laminated / K_n lattices), from memory, cross-checked against Cohn-Elkies Table 3 and
  (n = 9..20) Kallus Table I. Reliable. Proven optimal among lattices: n <= 8, n = 9, n = 24.
* n = 25..36: values equal Cohn-Elkies Table 3 ("best packing known", 2001); that these are lattices and their names
  (B27, B28, B29, Q30..Q34, B35, Ks36) are from memory. Kissing numbers only where I am confident.
* n = 37..47: **lower bounds only** -- laminated lattices via delta(Lambda_{24+k}) = 2^{k/2} delta(Lambda_k) (rule recalled
  from SPLAG Ch. 6; anchors Lambda_32 = 1, Lambda_40 = 16, Lambda_48 = 2^12 verified by the |u|u+v| construction). The
  Nebe-Sloane table very likely lists denser lattices in several of these dimensions (e.g. Mordell-Weil lattices near 44-47);
  do not treat these rows as records.
* n = 48, 56, 72: exact from (minimum, determinant) of extremal even unimodular lattices. n = 64: Ne64, delta = 3^16 [MEM].
  n = 80, 128: log2 delta = 40.14 and 97.40 for MW80, MW128 [MEM, 2-decimal precision]; rigorous fallbacks 2^40 (Bachoc-Nebe
  L80/M80) and 2^64 (BW128) are in the source column.
* Pattern check for the project: V/n^2 for the best known lattice is 0.6-1.0 for n <= 11, 3.8 (n=16), 56 (n=24), 2832 (n=48),
  1.3e5 (n=72): the records are exponentially, not polynomially, above the n^2 scale in every accessible dimension.

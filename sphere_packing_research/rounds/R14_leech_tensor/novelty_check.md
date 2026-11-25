# 14회차 문헌 점검 (novelty check) — 하위 에이전트 보고 원문

(2026-11-25 09:17–09:22 UTC 실행. 아래는 문헌 점검 에이전트의 보고를 그대로 옮긴 것이다. 에이전트가 직접 연 자료와 기억/검색 스니펫에 의존한 부분이 항목별로 구분되어 있다. 검색 스니펫에만 근거한 미확인 정보는 본 보고서 본문에 사용하지 않았다.)

## Summary

Literature check (run 2026-11-25, about 12 minutes). Bottom line: (1) No paper found states a Klartag-type c*n^2 bound with an explicit constant above 1/(2e), and none uses Hurwitz, icosian or other division-algebra module lattices with Klartag's method. The only explicit constant is Rezvyakova's 1/(2e) - o(1) for plain lattices; both n^2 log log n papers (Abuya-Gargava-Zhao, Gargava-Groebner-Simkin) are purely cyclotomic and give only an unspecified universal constant. (2) No paper found connects the Leech lattice, SL_2(p) Weil-representation / Gross lattices, or Mordell-Weil lattices of y^2 = x^p - x with Klartag's method or Venkatesh-type averaging. The closest precedent for a non-free action is Gargava-Groebner-Simkin, where a cyclic group acts on a multi-cyclotomic module with orbits of varying size. (3) Not answered from sources: the minimum of Gross's rank 2(p-1) lattices could not be checked because www.math.rwth-aachen.de (Nebe's catalogue) and www.ams.org (Gross's JAMS paper) do not resolve from this session and would need to be allowed in the network settings. The norm-6 finding in rank 72 is consistent with the literature only by indirect inference. (4) No Rogers-type higher-moment formula for modules over orders in non-commutative division algebras was found. The division-algebra case has only the first moment (Gargava 2021, Theorem 16; attributed upstream to Weil). The higher-moment papers checked (Kim, Hughes, Gargava-Serban-Viazovska, and the 2025 follow-up with Viglino) are all number-field only.

Coverage limits: the arXiv search was an abstract-phrase search on 'lattice packing(s)' / 'lattice sphere packings' (first 50 of 106 hits by date, through 24 Nov 2026), so a paper worded differently could be missed. The last third of the Gargava-Groebner-Simkin HTML and Gargava's EPFL thesis were not opened. Most pages were read through a summarising fetch tool, so quoted theorem statements are its rendering; five PDFs were read as extracted text and are marked as such. I did not check the correctness of the proposed constants 3/(2e) and 15/(4e), only whether they appear in the literature.

File note: despite the no-files instruction I downloaded five arXiv PDFs into a scratchpad subfolder for text extraction and deleted that folder before finishing. The fetch tool also auto-saved one PDF under /root/.claude/projects/-home-user-Mathathon/57e53de1-9d0a-5c79-adaf-33a30f7db1ad/tool-results/, which I left in place.

## Findings

### (1a) Does Abuya-Gargava-Zhao (arXiv:2606.05105) state explicit constants or use Hurwitz / icosian / division-algebra lattices?

**Answer.** No. Theorem 1.1: there is a universal c > 0 and an infinite sequence of dimensions N with Delta_N^lat >= c N^2 log log N 2^{-N}. Theorem 2.2: for large m there is a mu_m-stable lattice of covolume vol(B_N) and an ellipsoid E with no nonzero lattice points and vol(E) >= c*m*N*r*vol(B_N), where r = ceil((log phi(m))^4) and N = phi(m)*r; this is c*(m/phi(m))*N^2, which is where the log log N comes from. c is never made explicit and no constant involving e appears. The setting is module lattices g*O_K^r over K = Q(zeta_m), with g in SL_r(K_infinity). Quaternions, icosians and non-commutative orders are not mentioned; 'Hurwitz' appears only in the bibliography (Vance's title). Division algebras appear only as a citation to Gargava's Siegel formula (his Theorem 16). The paper states the mu_m action is 'free away from the origin'. Its only moment tool is a first-moment number-field Siegel formula (Theorem 2.1), attributed to Venkatesh (rank 2), Gargava (division algebras) and Kim (first-moment case). Context from a search snippet only: the first author's name is a pseudonym, and Aubrun's Bourbaki note reportedly describes the preprint as an AI-generated argument.

**Source.** https://arxiv.org/abs/2606.05105 and https://arxiv.org/html/2606.05105 (v2, 4 June 2026)

**Seen or memory.** Seen: abstract page and full HTML, both through the summarising fetch tool. The remark about AI generation and the pseudonym context is from a search snippet, not opened.

### (1b) Does Gargava-Groebner-Simkin (arXiv:2610.12619) state explicit constants or use non-commutative structures?

**Answer.** No explicit constant. Theorem 1.1: there are an absolute c > 0 and d_0 such that Delta_d^L >= c d^2 log log d 2^{-d} for every d >= d_0. The algebra is commutative: A = Z[X]/(F) with F a product of cyclotomic polynomials, averaged over G = product of SL_{r_n}(K_n tensor R) modulo an arithmetic subgroup. Theorem 3.1 gives an admissible gluing family with M/C >> log log d; constants such as delta_0 = 1/100 are local to lemmas. No quaternion, Hurwitz, icosian, division-algebra, Leech, SL_2(p) or Weil-representation mention in the part read. The moment tool is a Siegel-type first-moment formula 'for vectors with prescribed support' (Theorem 2.3); Rogers is not cited in the part read.

**Source.** https://arxiv.org/abs/2610.12619 and https://arxiv.org/html/2610.12619 (8 Oct 2026)

**Seen or memory.** Seen: abstract page and the first 100,000 of about 151,000 characters of the HTML, through the summarising fetch tool. The appendices and reference list were not read.

### (1c) What does Rezvyakova (arXiv:2608.22590) prove, and does it mention symmetric versions or larger constants?

**Answer.** Theorem 1: Delta_n >= (1/(2e) - o(1)) n^2 2^{-n} as n -> infinity, for the optimal lattice packing density. The proof produces a lattice-point-free ellipsoid with det A_t <= 4e^2 n^{-4}, giving covolume at most kappa_n * 2e * n^{-2} (1+o(1)); in the 'V_n >= c n^2' normalisation this is c = 1/(2e). It is an 8-page refinement of Klartag's argument for plain lattices only. No symmetric, cyclotomic, quaternionic, Hurwitz or icosian version, no mention of 3/(2e) or 15/(4e), and no remark on optimality of 1/(2e).

**Source.** https://arxiv.org/abs/2608.22590 and https://arxiv.org/html/2608.22590 (23 Aug 2026)

**Seen or memory.** Seen: abstract page and HTML, through the summarising fetch tool.

### (1d) Does any later or other 2025-2026 paper state V_n >= (3/(2e)) n^2 (Hurwitz), (15/(4e)) n^2 (icosians), or any explicit constant above 1/(2e)?

**Answer.** None found. An arXiv abstract search sorted by announcement date (first 50 of 106 hits, reaching 24 Nov 2026) shows only five Klartag-method papers: Klartag 2504.05042 (c n^2 2^{-n}, c unspecified; published in Inventiones 2026), Aubrun 2606.13313 (Bourbaki expose; abstract gives no explicit constant), Abuya-Gargava-Zhao 2606.05105, Rezvyakova 2608.22590, and Gargava-Groebner-Simkin 2610.12619. No 2025-2026 entry involves quaternions, Hurwitz, icosians or division algebras together with Klartag's method. Other recent hits are unrelated (for example Zubrilina 2610.35830 on upper bounds in small dimensions, Gargava 2607.11695 on counterexamples in dimensions 8 and 9). The largest explicit constant found is 1/(2e). The pre-Klartag division-algebra results are linear-order: Gargava 2107.04844 Theorem 21 gives density at least #G_0 * 2^{-d} in dimension d = 2 dim_Q D, and Theorem 30 gives C d (log log d)^{7/24} along a sequence, with lattices symmetric under a non-commutative finite group.

**Source.** arXiv abstract search for 'lattice packing' / 'lattice packings' / 'lattice sphere packings', newest first; https://arxiv.org/abs/2504.05042; https://arxiv.org/abs/2606.13313; https://arxiv.org/pdf/2107.04844

**Seen or memory.** Seen: the arXiv search listing (page 1 only) and the Aubrun abstract through the fetch tool; Gargava 2107.04844 Theorems 21 and 30 read in extracted PDF text. Limitation: a paper whose abstract avoids the searched phrases, or one on page 2 or later (all older), would be missed. Aubrun's full text was not opened.

### (2) Is there any paper using non-free group actions, the Leech lattice, or SL_2(p) Weil-representation lattices (Gross's globally irreducible lattices, Mordell-Weil lattices of y^2 = x^p - x) with Klartag's method or Venkatesh-type symmetric averaging?

**Answer.** None found. Targeted searches combining Klartag or Venkatesh with 'Leech lattice', 'non-free', 'Weil representation' and 'globally irreducible' returned only the cyclotomic papers. Abuya-Gargava-Zhao explicitly uses a free action. Gargava 2107.04844 uses finite subgroups of division algebras from Amitsur's classification (binary tetrahedral, octahedral and icosahedral groups, and the families G_{m,r} and T* x G_{m,r}), which act freely because they sit inside D*. The closest precedent for a non-free action is Gargava-Groebner-Simkin: a cyclic group of order M acts on a multi-cyclotomic module with orbits of varying size (Lemma 2.1 gives #Orb(v) = M_S), with the remark that 'vectors in a full cyclic orbit impose only one constraint on an invariant quadratic form', and the Siegel formula is stratified by support. That is still an abelian, cyclotomic setting. Recent Leech-related packing papers found (cross-sections and antipode constructions from Lambda_24 and P_48p) are fixed-dimension constructions unrelated to these methods.

**Source.** Web searches; https://arxiv.org/html/2606.05105; https://arxiv.org/html/2610.12619; https://arxiv.org/pdf/2107.04844; arXiv listing entries 2607.20359 and 2505.02394

**Seen or memory.** Seen: the two HTML papers through the fetch tool and the Amitsur table in Gargava's extracted PDF text. Absence is based on web search results and the arXiv listing, so it is not exhaustive. In the extracted Amitsur table the dimension entries for the binary octahedral and icosahedral rows read 16 and 20, which I could not reconcile with the 8-dimensional quaternion algebras over Q(sqrt 2) and Q(sqrt 5); likely a table-extraction artefact, unchecked.

### (3) What is known about the minimum of the even unimodular rank 2(p-1) lattices with SL_2(p) symmetry (Gross 1990; Nebe 1998; Tiep; Dummigan), especially rank 72 (p = 37) and rank 120 (p = 61)?

**Answer.** Not answered from primary sources. Nebe's lattice catalogue (www.math.rwth-aachen.de) and the AMS page for Gross's paper (www.ams.org) do not resolve from this session, and no search result stated the minimum of the rank-72 or rank-120 Gross lattice. Confirmed bibliographic facts: Gross, 'Group representations and lattices', J. Amer. Math. Soc. 3 (1990) 929-960; Nebe, 'Some cyclo-quaternionic lattices', J. Algebra 199 (1998) 472-498, DOI 10.1006/jabr.1997.7163, with keywords quaternion algebras, Mordell-Weil lattices, unimodular lattices, extremal modular lattices. Confirmed general bound: an even unimodular lattice of rank n has minimum at most 2*floor(n/24) + 2, so 8 in rank 72 and 12 in rank 120. Confirmed: existence of an extremal rank-72 lattice was a long-standing open problem until Nebe's 2010 construction, whose automorphism group contains (PSL_2(7) x SL_2(25)):2, not SL_2(37). Inference, not a sourced statement: since the SL_2(37) lattice has been available since 1990, it was evidently not known to have minimum 8, which is consistent with your norm-6 vectors. From memory only: the family requires p = 1 mod 4 and gives E_8 at p = 5, the Leech lattice at p = 13 and a rank-32 lattice of minimum 4 at p = 17; Nebe's 1998 abstract says the rank 2(p-1) faithful SL_2(p) lattices are described and studied by computer for small p, yielding a new extremal rank-48 lattice; no extremal even unimodular lattice of rank 120 is known. I could not confirm whether the minimum 6 for p = 37 is already recorded in Nebe 1998 or the catalogue. Dummigan's and Tiep's papers were not located.

**Source.** Search snippets from the MaRDI portal record for Nebe 1998, the AMS record for Gross 1990, Nebe's 'An even unimodular 72-dimensional lattice of minimum 8' and Nebe-Parker 2014; blocked: https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/unimodular.html and https://www.ams.org/jams/1990-03-04/S0894-0347-1990-1071117-8/

**Seen or memory.** Search snippets plus memory; no primary source opened. Both host failures were DNS resolution failures consistent with the session's network allowlist. The p = 5, 13, 17 values, the Nebe abstract wording and the rank-120 statement are memory and unverified.

### (4a) Kim, 'Adelic Rogers integral formula' (arXiv:2205.03138): exact generality?

**Answer.** Number fields only. Theorem 1.2: F a number field, X_n = GL(n,F)\G_n in the adelic setting, 1 <= k < n, and f a Borel integrable function on (A_F^n)^k; this is the full Rogers formula for all k < n over the adeles of a number field, plus second-moment formulas in 'a few important cases'. There is no division-algebra or quaternion setting: the only occurrence of 'division' in the paper is a bibliography entry (Tamagawa, 'On the zeta-functions of a division algebra'). Abuya-Gargava-Zhao note that Kim's first-moment result is stated under a class-number-one assumption.

**Source.** https://arxiv.org/abs/2205.03138 and https://arxiv.org/pdf/2205.03138

**Seen or memory.** Seen: abstract page through the fetch tool, and the PDF as extracted text (keyword search plus the theorem statement lines). The class-number-one remark is as reported in the AGZ HTML.

### (4b) Gargava-Serban-Viazovska (arXiv:2308.15275) and Gargava-Serban-Viazovska-Viglino (arXiv:2510.12893): exact generality?

**Answer.** Number fields only. 2308.15275 treats lattices that are modules over the ring of integers O_K of a number field K: moments of the number of omega_K-tuples of lattice points in a ball of volume V converge to Poisson moments with mean V/omega_K, as the rank grows or as the degree grows under a uniform Weil-height lower bound. Division rings are mentioned once, in the introduction, as prior first-moment work (citing Vance, Venkatesh, and Gargava's 'Lattice packings through division algebras', Math. Z. 303(1):1-32, 2023). 2510.12893 sharpens the variance estimates for module lattices over arbitrary number fields with emphasis on cyclotomic fields; its text contains no occurrence of 'division', 'quaternion', 'non-commutative' or 'Hurwitz'.

**Source.** https://arxiv.org/abs/2308.15275, https://arxiv.org/pdf/2308.15275, https://arxiv.org/abs/2510.12893, https://arxiv.org/pdf/2510.12893

**Seen or memory.** Seen: abstract pages through the fetch tool, and both PDFs as extracted text (keyword search; the introduction passage of 2308.15275 was read directly).

### (4c) Hughes 2023: exact generality?

**Answer.** Number fields only, S-arithmetic. The paper is Nathan Hughes, 'Mean Values over Lattices in Number Fields and Effective Diophantine Approximation' (arXiv:2306.02499, 4 June 2023). Theorem 4, attributed to Weil, is the first-moment formula on K_S^d for a number field K. Section 2.2, 'Higher Moment Mean Value Formulas on K_S^d', generalises it to sums over linearly independent r-tuples of vectors, with the range 1 < r < d-1 stated in the set-up. No division-algebra or non-commutative setting: a keyword search of the text for 'division', 'quaternion' and 'non-commutative' found nothing.

**Source.** https://arxiv.org/abs/2306.02499 and https://arxiv.org/pdf/2306.02499

**Seen or memory.** Seen: abstract page through the fetch tool, and the PDF as extracted text (keyword search plus the opening of Section 2.2). The full statement of the higher-moment theorem was not read.

### (4d) Is there a Rogers-type higher-moment formula for modules over orders in non-commutative division algebras (Weil 1946, Gargava, Vance)?

**Answer.** None found; only the first moment is available. Gargava, 'Lattice packings through division algebras' (arXiv:2107.04844; Math. Z. 303 (2023)), Theorem 16: for a Q-division algebra D with order O, G = SL_k(D_R), Gamma = SL_k(O), k >= 2, and f in C_c(D_R^k) (extended to Riemann-integrable compactly supported f in Remark 17), the integral over G/Gamma of the sum of f(v) over nonzero v in g*O^k equals the integral of f over D_R^k, with O^k of covolume 1. This is a Siegel mean value theorem only: a keyword search of the paper found no occurrence of 'Rogers', 'higher moment' or 'second moment'. Weil's 1946 result is cited by Hughes and by Abuya-Gargava-Zhao as the general source of first-moment formulas of this kind; I did not open Weil's paper. Vance's Hurwitz-lattice paper (Adv. Math. 2011) is, from memory, also first-moment only. Other higher-moment generalisations that surfaced (Han, S-arithmetic; Alam-Ghosh-Han, affine and congruence; Kelmer-Yu, symplectic; Stroembergsson-Soedergren, lattice and dual) are not for division algebras. Gargava's EPFL thesis (TH10488, 'Mean value theorems for collections of lattices with a prescribed group of symmetries') has a section on Weil's generalisation but was not opened, so a higher-moment statement there cannot be ruled out.

**Source.** https://arxiv.org/pdf/2107.04844; https://arxiv.org/pdf/2306.02499; https://arxiv.org/html/2606.05105; search snippets for https://infosciences.epfl.ch/record/306984/files/EPFL_TH10488.pdf, arXiv:1910.01824, 2111.00848, 1802.09645, 2211.05454

**Seen or memory.** Seen: Gargava 2107.04844 Theorem 16 and Remark 17 in extracted PDF text, with keyword search. Weil 1946 and Vance 2011 are from memory and from citations in the opened papers. The thesis and the other higher-moment papers are search snippets only.

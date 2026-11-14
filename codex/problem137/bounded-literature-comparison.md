# Bounded primary-literature comparison after the exact orbit proof

Status: five promising source/citation-chain checks, deliberately bounded. Four full primary texts were retrieved and their relevant sections inspected; the fifth was verified through its primary university repository but its PDF remained blocked. This is not a systematic novelty review. The exact arithmetic valuation was already known; the possible contribution under assessment is a particular fully combinatorial proof via Sylow actions, micro-depth skeletons, and explicit cycle splicing. Its historical novelty remains unresolved.

The earlier baseline comparison is konvalinka-comparison.md, which actually reads Konvalinka's Proposition 2/Proposition 7 and distinguishes arithmetic exactness from the tree-orbit bound. The official statement is ../sources/problem137-official.md.

## 1. Postnikov–Sagan, weighted binary Catalan numbers (2007)

Alexander Postnikov and Bruce E. Sagan, *What power of two divides a weighted Catalan number?*, J. Combin. Theory Ser. A 114 (2007), 970–977, DOI https://doi.org/10.1016/j.jcta.2006.09.007 ; primary preprint https://arxiv.org/abs/math/0601339 . Full text retrieved to ../sources/postnikov-sagan-2007.pdf and .txt; Sections 1–2 read.

Theorem 2.1 proves equality with the ordinary Catalan 2-adic valuation when b(0) is odd and 2^{r+1} divides Delta^r b(x) for every r>=1 and nonnegative x. Their vertex statistic is the number of *left* edges from the root, equivalent to their Dyck ascent-height model, not ordinary tree depth. Lemmas 2.3–2.4 combine the known minimum-orbit count with an odd, integer orbit-average weight. This is positive historical evidence for the binary minimum-orbit/weighted strategy; it is not merely related background. Our binary mechanism and any claim that invariant/average weights preserve valuations must be placed in this established context.

Specific distinction: their weights vary inside reflection orbits, so they prove integrality and oddness of the average using difference-operator conditions. The elementary ordinary-depth extension proved in depth-weighted-extension.md instead uses weights exactly invariant on the chosen group orbits. It does not subsume their path-weight theorem.

## 2. Gao–Gu, arithmetic of weighted Catalan numbers (2021; preprint 2019)

Yibo Gao and Andrew Gu, *Arithmetic of weighted Catalan numbers*, J. Number Theory 226 (2021), 213–242, DOI https://doi.org/10.1016/j.jnt.2021.03.007 ; primary preprint https://arxiv.org/abs/1908.03914 . Full text retrieved to ../sources/gao-gu-2021.pdf and .txt; Sections 2–3 read, including their minimum-orbit structure theorem and prime-power generalization.

Their Theorem 2.1 weakens Postnikov–Sagan's binary hypotheses. Section 3 uses q=p^k trees, full child-permutation symmetry, and weights b(i) where i counts non-right edges. Theorem 3.1 gives the stronger congruence

C_n^{(q)}(b) congruent C_n^{(q)} modulo p^{xi+k},
xi=(s_p((q-1)n+1)-1)/(p-1),

under b(0)=1 modulo q, q^2 dividing Delta b, and q^r dividing Delta^r b for r>=2. This is meaningful subsequent progress beyond Konvalinka, not a literature gap to ignore. Their technique analyzes orbit-average functions modulo q and contracts complete q-ary subtrees; the argument for the prime-power theorem is stated as a parallel reduction to small sizes.

For the narrow proof-comparison question: xi is identified with the already known unweighted valuation. Taking b identically one in that congruence is an identity, so that specialization alone is not a new independent explanation of unweighted exactness. The inspected Sections 2–3 do not give our Sylow p-subgroup/micro-depth-residue skeleton/cycle-chain enumeration. They do establish substantial overlap in orbit reduction and weighted noncancellation ideas. This bounded inspection does not prove our particular combination absent from all versions or subsequent literature.

## 3. Junkyu An, MIT thesis (2010)

Junkyu An, *Combinatorial enumeration of weighted Catalan numbers*, Ph.D. thesis, MIT, 2010, supervised by Alexander Postnikov. Primary repository record https://hdl.handle.net/1721.1/64609 (current landing page https://dspace.mit.edu/entities/publication/c2830ec4-fd0a-49b1-8d80-d37f383caee2 ).

The primary abstract states that Chapter 3 proves Konvalinka's polynomial-weight conjecture for binary weighted Catalan valuations and extends it to more general weight functions. Gao–Gu's introduction independently describes this as an arithmetic resolution of the polynomial case and discusses a remaining broader combinatorial conjecture. The repository PDF endpoint returned HTTP 429 on the bounded attempts; the thesis full text was **not read**, so no claims about a particular internal theorem, group action, or absence of a typed-skeleton proof are warranted. This is a concrete follow-up source, not negative evidence for novelty.

## 4. Stănică, prime-power Catalan residues (2003; preprint 2000)

Pantelimon Stănică, *p^q-Catalan numbers and squarefree binomial coefficients*, J. Number Theory 100 (2003), 203–216, DOI https://doi.org/10.1016/S0022-314X(02)00125-7 ; primary preprint https://arxiv.org/abs/math/0010148 . Full preprint retrieved to ../sources/stanica-2003.pdf and .txt; preliminary Lucas/Granville/Kummer statements and the initial divisibility/residue theorems read. The arXiv PDF's automatically generated date must not be mistaken for its publication year.

Theorem 4 already gives the prime-arity p-unit repunit classification and the nonzero residue one, using Lucas. Later results treat p^2 arity residues and more general prime-power divisibility with arithmetic tools, including factorial p-unit products and carry counts. Thus our repunit unit patterns and their residue one are not new arithmetic results. This paper is earlier than Konvalinka, but is a useful positive historical anchor from the citation chain. The inspected proof uses arithmetic congruences, not the expanded Sylow action or typed minimum-orbit count.

## 5. Bobrowski–He–Shiue, generalized Catalan/Raney divisibility (2019)

Jacob Bobrowski, Tian-Xiao He and Peter J.-S. Shiue, *Divisibility of Generalized Catalan Numbers and Raney Numbers*, J. Computational Analysis and Applications 26 (2019), no. 6, 1047–1058. Author/institution-hosted full text https://scholars.iwu.edu/ws/portalfiles/portal/39724414/144%20divisibility%20of%20generalized%202019.pdf , retrieved to ../sources/catalan-raney-2019.pdf and .txt; divisibility results and references inspected.

The paper relates generalized Lobb, Catalan and Raney numbers, deriving integer divisors through binomial identities and known arithmetic divisibility results. Its discussion includes the prime-arity exceptional repunit indices via a result of Sury. The inspected arguments and bibliography do not develop the Sylow/micro-depth minimum-orbit proof at issue here. This establishes its different immediate scope; it does not establish that no related proof exists elsewhere.

## Precise present conclusion

There is positive primary evidence that the binary orbit exactness proof, weighted orbit-average strategy, prime-power weighted congruence improvements, and the arithmetic unit-index/residue theorems are established. The current writeup should attribute those facts explicitly. The candidate project result is the complete proof assembly that counts minimum Sylow orbits in the depth-constrained p-ary expansion, with a nonzero minimum-orbit residue and an explicit type-preserving cycle-splicing bijection. This bounded comparison has not settled the historical status of that assembly.

Before a public novelty claim, a targeted expert/literature check of multitype tree enumeration, combinatorial Lagrange inversion and later citations to Gao–Gu/Konvalinka would be appropriate. For this 40-hour research project, a fully checked and understandable reconstruction is already a defensible deliverable even if a precedent is found.

## A concrete extension actually proved here

See depth-weighted-extension.md: for q=p^e and ordinary-depth weights b_h=u modulo p with p not dividing u, the exact valuation remains D and the leading unit becomes u^n(-1)^D/product_j d_j!. In particular any odd depth-weight sequence works for binary trees. The distinction from left-/non-right-edge weighting is essential; it is illustrated by n=2 and weights b_0=1,b_1=3, where ordinary-depth weighting gives 6 but left-edge weighting gives 4. Independent tree-recursion verification passed 254 cases in check_depth_weights.py. This extension too carries no historical novelty claim.

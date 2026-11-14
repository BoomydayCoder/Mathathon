# Official #137 alignment and comparison with Konvalinka's primary paper

Status: official statement and primary paper now retrieved and read. This supersedes the earlier access-blocked status of the initial investigation, while retaining those files as a record of independently developed approaches. No historical novelty claim.

## Exact task alignment

The official source is ../sources/problem137-official.md with raw JSON and retrieval provenance maintained by the parent. It asks for a combinatorial explanation of exact p-adic valuations for q=p^k Fuss–Catalan numbers, beyond the divisibility bound obtained from Konvalinka's tree symmetries; alternative actions, minimum orbits and special cases are expressly within scope. Our independent prime-power investigation therefore matches the exact task, rather than only its earlier informal summary.

The proposed complete result is

L=(q-1)n+1=sum_j d_j p^j, D=(sum_j d_j-1)/(p-1),
v_p(T_q(n))=D,
p^{-D}T_q(n) congruent (-1)^D / product_j d_j! (mod p).

The route is a full Sylow p-group on expanded p-ary trees, exact minimum-orbit skeleton enumeration, and a nonzero residue for their count. See prime-power-sylow-audit.md. The later cycle-chain bijection being developed in the companion review replaces its initial Lagrange–Good counting machinery by an explicit combinatorial encoding.

## Primary source retrieved

Matjaž Konvalinka, *Divisibility of generalized Catalan numbers*, J. Combin. Theory Ser. A 114 (2007), no. 6, 1089–1100, DOI https://doi.org/10.1016/j.jcta.2006.11.003.

Author-hosted primary text: https://users.fmf.uni-lj.si/konvalinka/wcatalan.pdf . Retrieved with ordinary Mozilla/5.0 user agent; local files ../sources/konvalinka-2007.pdf and ../sources/konvalinka-2007.txt. The relevant propositions and proof were actually read, rather than inferred from the abstract or secondary summaries. The author's preprint page numbering below runs 1–12.

## What the paper establishes

* Proposition 2, pp. 2–3, states the exact digit-sum formula and proves it by counting powers of p in factorials.
* Proposition 7, p. 6, gives the corresponding divisibility lower bound for every orbit of its tree action. The paragraph immediately after it explicitly distinguishes this from proving exactness of the total tree count.
* Lemmas 8–10, pp. 7–10, analyze minimum orbits through complete subtrees and equivalence classes of children. Their local orbit factor is q! divided by the factorials of equivalent-child multiplicities. Thus the action uses all permutations of q children, namely the full rooted-tree automorphism group, rather than only cyclic rotations or a p-group.
* At the end of that proof, p. 10, the paper obtains the minimum orbit exponent bound and uses the earlier arithmetic Proposition 2 to identify the exponent exactly. It does not count the minimum contributions modulo p to rule out cancellation in the total count by this group argument alone.
* Its main theorem concerns sufficient conditions on weighted generalized Catalan numbers for divisibility. Our current complete proof is unweighted; weighted extensions remain a further direction and should not be asserted from the unweighted orbit count.

## How the present proof changes the mechanism

1. Expand every q=p^e vertex into a perfect p-ary block of height e. The q-tree set becomes the p-tree set with leaf depths divisible by e. It is stable under the iterated p-wreath group on a finite ambient complete tree.
2. This is a Sylow p-group in the appropriate full symmetry group. Every orbit has p-power cardinality, so a minimum orbit's contribution after dividing by p^D is exactly one, rather than an additional uncontrolled p-unit factor.
3. The standard p-ary induction gives all exponents at least D and identifies equality by carry-free leaf counts. Contract maximal perfect blocks to a labeled skeleton. Label j appears d_j<p times and a leaf at depth i must have i+j=0 modulo e.
4. These digit multiplicities prevent p identical child orbit types anywhere in the skeleton, forcing every equality skeleton orbit to have exactly p^D elements.
5. Count the labeled plane skeletons by depth type. With l_i=sum_{j congruent -i}d_j and k_i internal counts satisfying k_i+l_i=delta_{i0}+p k_{i-1}, their number is N=product_i (p k_i)!/[product_i k_i! product_j d_j!]. The explicit cycle-chain encoding gives a purely combinatorial proof of this formula; the initial independent audit used an explicit Lagrange–Good specialization as a separate check.
6. The minimum-orbit count M=N/p^D is a p-unit: (p k)!/(p^k k!) is exactly the product of integers not divisible by p up to p k, whose residue is (-1)^k by inverse pairing. Since d_j<p, their factorials are units. This proves the noncancellation needed for exactness, without importing factorial valuations.

The remaining factorials count finite permutations and label assignments. The proof does not evaluate v_p((qn)!) or invoke Legendre/Kummer to infer the desired answer; its p-divisibility conclusion comes from explicit orbit sizes and the unit product. This distinction is material to the official request.

## What may and may not be claimed

The above supplies a concrete candidate combinatorial completion of the exactness step for the official project, with independent algebraic and computational audits. The binary base mechanism is already Deutsch–Sagan's established argument and must be credited. The use of p-groups, skeletons, depth-residue counts or cycle-chain bijections may also have precedents; the source comparison does not establish their absence. Before any historical novelty claim, inspect later literature on generalized/weighted Catalan valuations, multitype tree enumeration and Joyal-type functional digraph encodings. A proof meeting the project's mathematical objective can be valuable regardless of historical novelty.

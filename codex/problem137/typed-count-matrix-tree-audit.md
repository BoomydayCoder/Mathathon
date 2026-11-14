# Independent alternative audit of the typed skeleton count

This replaces the Lagrange–Good calculation by a directed Matrix–Tree calculation. The companion cycle-chain bijection is the more elementary main exposition; this provides an independent combinatorial-algebraic check. The Matrix–Tree theorem itself has a standard combinatorial cycle-cancellation proof. No factorial-valuation theorem is used.

Use the notation of prime-power-sylow-audit.md, with D>0, root type 0, internal type counts k_i, leaf type counts l_i and slot constraints k_i+l_i=delta_{i0}+p k_{i-1}. Label every vertex distinctly within its type and internal/leaf status. Fix one type-0 internal label r as the root. Give an internal vertex v variable x_v; leaves have variable zero. Put S_i=sum_{internal v of type i}x_v.

Draw an allowed directed edge from any nonroot vertex of type i to any internal parent of type i-1, with weight x_parent; exclude self-edges when e=1. The directed Matrix–Tree theorem says the root-deleted outgoing Laplacian determinant sums these edge-weight products over all trees pointing toward r. An exponent of x_v records the number of children of v.

**Claim.** That determinant is

x_r S_0^{p k_0-1} product_{i>0} S_i^{p k_i}.

If all k_i>0, write the Laplacian as a diagonal matrix with diagonal S_{i-1} on type-i rows, minus a block-cyclic rank-one matrix: each row subtracts the same vector of allowed parent weights. This expression handles self-edge exclusion too, since in e=1 the diagonal becomes S_0-x_v. The diagonal determinant is product_i S_i^{p k_i}, by the slot constraints after deleting the root. The finite-rank determinant correction has only the identity and full cyclic term; its value is

1-[(S_0-x_r) product_{i>0}S_i]/[product_i S_i]=x_r/S_0.

This proves the claim. If some k_i=0, the corresponding next-type row block is empty and the cyclic determinant term is absent. The nonnegative slot constraints force k_0=1 in this case, so S_0=x_r and the same formula holds. D=0 is separately the single labeled leaf and needs no determinant.

Now extract coefficient product_internal x_v^p. The result is

[(p k_0-1)!/((p-1)!(p!)^{k_0-1})]
  product_{i>0}[(p k_i)!/(p!)^{k_i}]
= product_i(p k_i)!/[k_0 (p!)^D].

These are rooted labeled trees with unordered children. Since labels are distinct, ordering the p children independently at each of D internal vertices multiplies by (p!)^D. The fixed-root labeled plane count is product_i(p k_i)!/k_0. Summing over the k_0 possible internal root labels cancels that denominator. Finally forget the internal labels and distinguish only the desired leaf labels j: divide by product_i k_i! product_j d_j!. This gives exactly

N=product_i(p k_i)!/[product_i k_i! product_j d_j!].

The division is by the number of labelings of a plane tree, whose ordered positions make this action free; no unproved quotient by tree automorphisms is hidden. Thus the typed count is independently consistent with both the explicit Lagrange–Good specialization and the proposed fully bijective functional-digraph encoding.

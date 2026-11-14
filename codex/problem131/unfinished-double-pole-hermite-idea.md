# Unfinished double-pole Hermite construction

Status: preliminary algebraic observation, checkpointed when this agent was redirected to zeta-three proof simplification. No approximation improvement, arithmetic bound, or irrationality claim. No independent review or priority search.

Let D(t) be a squarefree rational polynomial whose roots are negative half-squares x_j=-(j+1/2)^2. Consider rational polynomials P satisfying

    2 P'(x_j)D'(x_j) - P(x_j)D''(x_j) = 0  (every j).

The constant Laurent term of P/D at x_j is zero. Therefore, for any two such polynomials P,Q, the rational function PQ/D^2 has no simple-pole partial fractions. Its only pole terms are c_j(P)c_j(Q)/(t-x_j)^2, where c_j(P)=P(x_j)/D'(x_j).

The already proved positive measure sigma_2 has rational polynomial moments and integral sigma_2(t)/(t+(j+1/2)^2)^2 dt = 2(zeta(5)-H_j^(5)). Consequently its Gram matrix on any linearly independent such rational functions has entries affine over Q in zeta(5), and is positive definite. This is an exact formal way to isolate the pure target without differentiating a quotient of Apéry values.

For P/D = S(t)+sum_j c_j/(t-x_j), the constraints reduce to

    S(x_j)+sum_(l != j) c_l/(x_j-x_l)=0.

The residue coefficient matrix is skew-symmetric. Invertibility, explicit bases, primitive denominators, and useful determinant decay were not established. This construction is separate from the finite-field Hermite ideal in outer-correction-Hermite-ideal.md.

A warning against assuming approximation quality: if S=1, residue cancellation gives exactly

    (P(t)/D(t))^2 = 1+sum_j c_j^2/(t-x_j)^2 >= 1  (t>=0).

Thus the constant-polynomial member, whenever it exists with real residues, is not a small positive test function. A successful determinant use would require an additional mechanism. The investigation stopped here, without numerical testing.

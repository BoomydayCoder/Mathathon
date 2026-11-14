# Independent audit and consequences of the discrete Green formulation

Status: independently checked derivation and exact executable rerun. No error found in discrete-green-motivation.md. The ground-state transformation and resistance identities are standard constructions; no historical novelty is asserted.

## 1. Audit of indexing and assumptions

The recurrence is (n+1)^3u_(n+1)=b_nu_n-n^3u_(n-1), with b_n=(2n+1)(17n²+17n+5). A_0=1,A_1=5 and B_0=0,B_1=6. The conductance is correctly indexed as c_n=(n+1)^3A_nA_(n+1), and r_n=1/c_n.

Substitution u_n=A_nv_n gives c_n(v_(n+1)-v_n)=c_(n-1)(v_n-v_(n-1)) for n>=1. For B, c_0=5 and the first voltage increment is6/5, so the current is6. Thus B_n/A_n=6sum_(j=0)^(n-1)r_j. Growth A_n>=5^n gives convergence.

Identification of the resistance with zeta(3) is not inferred from electrical language. The separate rational-function proof gives F_n=2(A_n zeta(3)-B_n) with |F_n|<=4+2/n. Since A_n grows, this identifies the limiting ratio and hence sum r_j=zeta(3)/6. The normalization factor2 is harmless for the boundedness claim but must not be confused with the network current6.

In the quadratic-form identity, the right side includes the edge from the last nonzero f_N to the grounded value f_(N+1)=0. Its contribution is necessary for the diagonal coefficient b_N. At index0, c_0/A_0²=5=b_0; no fictitious edge to -1 is needed.

For J_N indexed0,...,N, the final edge conductance c_N appears in the energy even though the matrix only has off-diagonal entries up to c_(N-1)/(A_(N-1)A_N). The recurrence supplies the final diagonal correctly. This establishes positive definiteness, not merely nonnegative definiteness.

The finite minimizer has v_n=(sum_(j=n)^N r_j)/(sum_(j=0)^N r_j) for0<=n<=N and v_n=0 afterward. It has f_0=1, energy1/(sum_(j=0)^N r_j), and identifies (J_N^-1)_(0,0)=B_(N+1)/(6A_(N+1)). The determinant formula det J_N=((N+1)!)³A_(N+1) has the same indexing and is correct.

All continued-fraction tail pivots are positive Schur complements of positive definite trailing blocks. This proves positivity for the actual nested fractions; an additional guessed numerical lower bound on every pivot is unnecessary. The stale gap in recurrence.md has been replaced with this argument.

The exact script check_discrete_green.py was independently rerun and passed for N=0,...,9. Its checks cover conductances, the ground-state diagonal coefficients, finite minimizers, inverse entries, and determinants. The general proof, not the finite range, establishes the theorem.

## 2. Full finite Green kernel

The same argument gives every inverse entry, not just the corner:

    (J_N^-1)_(i,j)=A_i A_j sum_(k=max(i,j))^N r_k,
    0<=i,j<=N.                                         (1)

To prove it, fix j and divide the proposed column by A_i. Its voltage is constant for i<=j and then decreases by A_j r_i across each edge i>=j. The net current is zero except at j, where it equals A_j. Conjugating back by A_i makes the matrix source exactly the coordinate vector e_j. The final voltage at N+1 is zero, so the boundary equation is also satisfied.

Let E_j=A_j zeta(3)-B_j. Taking entrywise limits yields

    G_(i,j)=A_i A_j sum_(k=max(i,j))^infinity r_k
           =A_i E_j/6 for i<=j,                         (2)

with the symmetric formula for i>=j. This identifies the growing solution A and decaying solution E as the two fundamental solutions in the Green kernel. No assertion about a chosen self-adjoint operator extension is needed to state these finite-matrix and entrywise-limit identities.

For i,j<=N the exact truncation error is

    G_(i,j)-(J_N^-1)_(i,j)
       =A_i A_j sum_(k=N+1)^infinity r_k.                (3)

Thus the error of the entire leading block is positive semidefinite and rank one. This is a useful stronger consequence of the scalar resistance-tail identity.

## 3. The decaying Apéry form is the infinite energy minimizer

Define the transformed energy

    E(f)=sum_(n>=0)c_n(f_n/A_n-f_(n+1)/A_(n+1))².

Consider sequences with E(f)<infinity, f_0=1, and f_n/A_n tending to zero. Weighted Cauchy–Schwarz gives E(f)>=6/zeta(3). Equality is attained uniquely by

    f_n^*=(A_n zeta(3)-B_n)/zeta(3).

Indeed its voltage increments are -6r_n/zeta(3), giving precisely the equality case. Conversely, equality in Cauchy–Schwarz forces those increments, and the boundary condition determines the voltage uniquely.

This specifies the variational domain. The boundary condition at infinity cannot be dropped: f=A has constant voltage and zero transformed energy, but is not in the stated affine domain. For finite-support sequences the boundary condition is automatic; the original sharp inequality is correct as stated. Its best constant is not attained by a nonzero finite-support sequence, only approached by the displayed finite minimizers.

For f^* the original untransformed quadratic sum also converges absolutely, since the resistance-tail estimate gives E_n=O((n+1)^(-3)/A_n). The equality interpretation is therefore consistent with the original quadratic form for this particular infinite minimizer.

## 4. A sharp inequality at every vertex

For a fixed j and every finitely supported f,

    Q(f)>=6 f_j²/[A_j(A_j zeta(3)-B_j)].                 (4)

The constant is sharp. In voltage coordinates, the energy-minimizing profile is constant to the left of j and drops to zero at infinity in proportion to the remaining resistances. Equivalently, the best evaluation constant is the diagonal Green entry G_(j,j)=A_j E_j/6 from (2). Truncating that profile proves sharpness within finite-support sequences.

This supplies a family of exact discrete inequalities, with the original 6/zeta(3) inequality as j=0. It does not provide an independent arithmetic denominator estimate or identify a new zeta value.

## 5. What the audit establishes

The network viewpoint rigorously replaces separate positivity arguments for the approximation error and continued-fraction pivots once the recurrence, positive A sequence, and bounded-form identification are supplied. It explains why a large denominator solution produces a small companion solution and shows that the latter is an energy minimizer.

It does not independently discover the cubic recurrence, prove binomial integrality, or prove irrationality without arithmetic denominator control. Those dependencies remain explicit in the main note and discovery guide.

# Second independent audit of the original Schur–Apéry bridge

Status: **accepted**. Reviewed in full:
[rank-one-original-apery-heun-rii-proof.md](rank-one-original-apery-heun-rii-proof.md).
This audit independently checks the all-index argument, the normalization of its interpolation input, its small-index cases, and the algebraic identification with Apéry's recurrence. No historical novelty is asserted.

The theorem concerns the original positive measure mu_3/D_(2h):
\[
 q_{h,3}=\frac{5B_h-B_{h-1}}{5A_h-A_{h-1}}\quad(h>=1).
\]
It does not replace that measure by formal negative-parameter Wilson orthogonality. The reviewed [boundary-defect identity](https://github.com/BoomydayCoder/Mathathon/blob/03cc75d694220d0fba888c637e573a29799e6d5b/codex/problem131/rank-one-boundary-independent-review.md) is a genuine input; the new Heun operator cancels its explicit defect.

## 1. Original normalization and interpolation input

The monic minimizer p_h exists uniquely by positive definite orthogonality. Its roots are positive and simple. Therefore Q_h(z)=(-1)^(h-1)p_h(-z) is monic with negative roots, Q_h(0)>0, and R_h=Q_h/Q_h(0) has positive leading coefficient L_h.

I reconstructed the interpolation input directly, without using the previously conjectured Apéry match. At the original Schur root q_h the original Gram matrix has the coefficient vector of p_h as a nullvector. In the right Cauchy basis its coordinates are
\[
 r_b=Q_h(b^2)/C_R'(b^2),\qquad
 C_R(z)=\prod_{b=h+1}^{2h}(z-b^2).
\]
The signs follow from C'(-b^2)=(-1)^(h-1)C_R'(b^2), and monicity gives sum r_b=1. The cross-Cauchy matrix has entries
(f_i(q_h)-f_b(q_h))/(b^2-i^2), so its nullvector condition gives precisely the left interpolation equations. The barycentric construction gives the right ones automatically.

Consequently N_h has degree at most h-1 and
N_h(j^2)=Q_h(j^2)[j^2 q_h+beta_j] at all 2h nodes. Thus
S_h=(N_h-q_h zQ_h)/Q_h(0) has degree at most h, interpolation values beta_j R_h(j^2), and leading coefficient -q_h L_h. These signs and degree bounds are proved all-index; no finite agreement or irrationality assumption is needed.

## 2. The Heun step and its two cutoffs

For the boundary form B(f,g)=ell(f Dg-g Df), applying the accepted defect identity to (Xf,g) and (f,Xg) gives
B_XD+DX=B_D, because X has endpoint values 0 and 1. Thus XD+DX-D is symmetric; the extra multiplication term is symmetric as well.

With n=h-1 and M=2h, the raising coefficient is exactly
2(k-n)(k-n+1). Its zeros at both n-1 and n are essential: the two polynomial spaces of degrees at most n-1 and n are invariant. Symmetry then makes the original orthogonal Q_h an eigenpolynomial. The residual below degree n is both in the lower space and orthogonal to that space, hence is zero by positivity. At h=1 the constant polynomial gives the claimed eigenvalue directly.

I independently derived the top two coefficients of D_M on monomials by expanding A_M(x) at infinity, rather than using the author's Newton action. This gives the same eigenvalue
eta_h=1-2h-2h^3(h+1).

The endpoint at x=0 yields
Hf(0)=(1-M)f(0)-2M f(1), hence
R_h(1)=h^2(h+1)/2. A separate direct evaluation at x=1 gives
\[
 \eta_hR_h(1)
 =A_M(1)[4R_h(4)-R_h(1)]-\kappa_hR_h(1),
 \qquad A_M(1)=4(1-M)/3.
\]
This independently recovers the stated R_h(4), including h=1, without the coefficient e_1 calculation. Both endpoint values are therefore fully normalized.

## 3. The varying-measure recurrence is valid

The interpolation of R_(h+1)/R_h at the two OLD nodes is legitimate because R_h has no positive roots. Dividing the numerator difference by
F_h(z)=[z-(2h-1)^2][z-(2h)^2] gives W of degree at most h-2.

The precise identity
\[
 F_h(-t)=D_{2h}(t)/D_{2h-2}(t)
\]
proves the first measure conversion, with positive sign. The R_h term is orthogonal because its test degree is at most k+1<=h-2. For R_(h+1), converting ell_(2h) to ell_(2h+2) multiplies the test by the TWO NEW factors; its degree is k+2<=h-1, within the next polynomial's orthogonality range. Hence W lies in the one-dimensional orthogonal complement spanned by R_(h-1).

For h=2 the lower orthogonality conditions are empty and W is a constant, which is exactly the required one-dimensional space. No recurrence is asserted at h=1. A possible zero W would still give a proportionality; the explicitly solved coefficient c_h is nonzero for h>=2, so the final recurrence is nondegenerate.

The endpoint coefficient matrix, in unknown order (a,b,c), has determinant
\[
 \frac38h^6(h+1)^4(12h^2-15h+5)>0\quad(h>=2).
\]
I solved this system independently and obtained the stated coefficients. Its quadratic factor has negative discriminant and positive leading coefficient, so no exceptional positive integer is hidden in this step.

## 4. Numerator recurrence, sources, and rescaling

The interpolation-numerator residual has degree at most h+1. At the first 2h-2 nodes all three interpolation values apply; at the last two OLD nodes the F_h term vanishes. It therefore has 2h distinct zeros. Since h+1<2h for every h>=2, it is zero. The smallest case has four zeros against degree at most three, so there is no h=2 exception.

Taking leading coefficients gives the same scalar recurrence for L_h and Y_h=q_h L_h. The sign -q_h L_h in S_h cancels correctly. There is no missing inhomogeneous source: this recurrence starts at h=2, and two separate initial values are supplied.

The factorial rescaling s_h=h^2/[12(2h)!] gives exactly the displayed alpha_h and gamma_h. I separately eliminated u_(h+1) and u_(h-2) from a free Apéry recurrence solution; the adjacent combination 5u_h-u_(h-1) satisfies precisely that rescaled recurrence. The elimination denominator -12h(12h^2-15h+5) is nonzero throughout h>=2.

Independent original-pole moments give
R_1=1, R_2=1+5z, q_1=5/4, q_2=577/480.
They match both rescaled Apéry initial pairs, including
Y_2=577/96 and 5B_2-B_1=1731/4.
Forward recurrence uniqueness proves the two leading-sequence identities. Since L_h>0, their quotient is valid and proves the theorem.

The final constant-term formula for Q_h and the leading-coefficient formula for N_h follow with the stated factorial factors. Equality with the moving-pole optimized error is then a consequence of the common quotient, not an input to this proof.

## Independent exact checks

[check_rank_one_heun_independent.py](computations/check_rank_one_heun_independent.py) imports no author's moment or Heun/RII routines. It checks:

- The eigenvalue from independent Laurent coefficients and both endpoint ratios, including direct evaluation at x=1.
- The free-h endpoint determinant, independently solved recurrence coefficients, factorial rescaling, and Apéry elimination.
- Original Gram moments from derivatives of the complete pole polynomial for h=1,...,6; the affine parameter coefficient and positivity of every leading Gram block.
- Exact negative root counts, all interpolation nodes and the leading numerator sign.
- The RII relation reconstructed by interpolation at the OLD nodes, its numerator counterpart, and both scalar normalizations.

Every check passed. These finite tests corroborate the all-index arguments above; they are not the basis of the acceptance.

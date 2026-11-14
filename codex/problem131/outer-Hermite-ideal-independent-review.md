# Independent review of the explicit outer Hermite ideal

Reviewer: `odd_zeta`, after independently deriving the earlier correction-matrix
kernel and its seven finite checks. Reviewed `outer-correction-Hermite-ideal.md`,
including its later converse and general jet form. No gap found in these claims.

The polynomial F has exactly the required CRT support: order two at removed
nodes, order six at unremoved nodes, and order two at zero. Its degree is
3p−1−4N, so its multiples below degree h have dimension d=K+3N−3p+1.
The local quotient H_f=Wf/D is polynomial and the stated derivative orders
kill all critical zero-weight and old −1/2-weight columns.

The converse is also valid. A fixed unremoved critical column receives a
correction contribution only from its own node. The paired ±a contribution
has the form [u³]ω(u)H_f(u)u^i with ω(0) a nonzero unit; the coordinate-change
factor for s=5 is −4a⁸. For ℓ=4 the order-three column first forces H_f(0)=0,
then the order-two column forces H_f'(0)=0. For ℓ=5 the order-three column
forces H_f(0)=0. Hence the selected-space kernel is exactly the ideal of
multiples of F, not merely a subspace of that ideal. This proves exact
correction rank within the selected space, not exact determinant valuation.

The general jet formula was checked by substituting z=√(a²−u)−a into the
formal-residue change of variables. It gives
ω(u)=−[b(u)(a+b(u))]^{s−1}η(u)/(s−1), in agreement with the displayed
nonzero constant. The square-unit reduction to an anti-diagonal form is valid
in odd characteristic; its isotropic ideals and stated Witt index follow.

For ζ(7), the proposed sufficient ideal requires orders three and eight and
has degree 4p−2−5N. Its degree deficit is negative in the stated band. The note
correctly limits this obstruction to that construction and does not claim
that all higher-odd kernels fail.

The parity lift and coefficientwise-in-X qualification agree with the earlier
independent correction-matrix argument. This review retains the source local
functional and row-bound dependencies; it does not audit the entire preprint.

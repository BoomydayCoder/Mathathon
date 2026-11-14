# Independent recurrence investigation

Status: elementary proofs below independently derived in this session; exact rational spot checks passed. These are background results, **not claims about the unretrieved official problem statement**. No historical novelty is asserted. This note deliberately does not identify the limiting constant with zeta(3) without a separate argument.

## A general discrete reduction-of-order theorem

Let a_n,c_n be nonzero real numbers, n>=1, and let two real sequences A,B obey

    a_n u_(n+1) = b_n u_n - c_n u_(n-1).

Assume A_n is nonzero for every n>=0. Put W_n=A_n B_(n+1)-A_(n+1) B_n. Direct substitution gives

    W_n=(c_n/a_n) W_(n-1),
    B_(n+1)/A_(n+1)-B_n/A_n = W_n/(A_n A_(n+1)).

Consequently, if the series on the right converges absolutely, the ratios tend to a finite L and

    L-B_n/A_n = sum_(j=n)^infinity W_j/(A_j A_(j+1)).             (1)

Conversely, given one nonvanishing solution A, choose an arbitrary W_0 and B_0, define the successive ratios by this sum, and multiply by A_n. The resulting B obeys the recurrence: its Casoratian has exactly the required multiplicative law, and rearranging that law yields the recurrence for B. Thus a second solution is constructed from one solution, rather than guessed.

This is the discrete counterpart of reduction of order for second-order differential equations. The recurrence encodes both a rapidly growing denominator and a small approximation error; (1) makes the mechanism visible.

## Specialization to the classical Apéry recurrence

Define, initially over the rationals,

    (n+1)^3 u_(n+1)=(2n+1)(17n^2+17n+5)u_n-n^3 u_(n-1), n>=1,
    A_0=1, A_1=5; B_0=0, B_1=6.

The preceding theorem gives exactly

    W_n=6/(n+1)^3,
    L-B_n/A_n=6 sum_(j=n)^infinity 1/((j+1)^3 A_j A_(j+1)).     (2)

Here L exists without appealing to zeta values. Indeed A_(n+1)>=5 A_n>0 follows inductively. If A_n/A_(n-1)>=5, the next ratio is at least (b_n-n^3/5)/(n+1)^3 >=5: after multiplication by 5, the last inequality reduces to

    144n^3+180n^2+60n >=0.

The initial ratio is 5. Thus A_n>=5^n and the series converges. More quantitatively, applying the same growth bound starting at n gives

    0 < L-B_n/A_n <= 5/[4 (n+1)^3 A_n^2],
    0 < A_n L-B_n <= 5/[4 (n+1)^3 A_n].                       (3)

Proof of the constant: A_(n+k) A_(n+k+1)>=5^(2k+1) A_n^2 and sum_(k>=0)5^(-2k-1)=5/24; multiply by 6.

These formulas are exact structural consequences of the recurrence. Identifying L=zeta(3), establishing integer A_n and sufficient denominator bounds for B_n, and making the final irrationality inequality are separate tasks. Recurrence alone is not an irrationality proof.

## Continued fractions: a second interpretation

Set b_n=(2n+1)(17n^2+17n+5). For each N>=1 the finite generalized continued fraction satisfies

    B_N/A_N = 6/(b_0 - 1^6/(b_1 - 2^6/(b_2 - ... -(N-1)^6/b_(N-1)))).

Thus its first terms are 6/5, 6/(5-1/117), and 6/(5-1/(117-64/535)).

Proof: multiply A_n and B_n by (n!)^3. The transformed sequences X,Y obey

    X_(n+1)=b_n X_n-n^6 X_(n-1),
    Y_(n+1)=b_n Y_n-n^6 Y_(n-1),
    (X_0,X_1)=(1,5), (Y_0,Y_1)=(0,6).

The numerator and denominator continuants of the displayed fraction satisfy precisely these recurrences and initial conditions (expansion of a tridiagonal determinant proves the continuant recurrence). Every actual nested denominator is positive. To prove this without an auxiliary numerical bound, let J be the tridiagonal matrix indexed 0,...,N-1 with diagonal b_n and neighboring entries -(n+1)^3. For a vector f extended by zero beyond N-1, the ground-state identity is

    f^T J f = sum_(n>=0) (n+1)^3 A_n A_(n+1)
                  (f_n/A_n-f_(n+1)/A_(n+1))^2.

It is strictly positive for every nonzero f, so J is positive definite. Backward elimination produces positive Schur-complement pivots, which are exactly the denominators of the nested continued fraction. Thus the earlier tail-positivity gap is closed. See [discrete-green-motivation.md](discrete-green-motivation.md), Sections 2 and 4, for the boundary check and full argument, and [discrete-green-independent-review.md](discrete-green-independent-review.md) for independent verification. This uses standard ground-state and matrix ideas; no historical novelty is asserted.

## A precise arithmetic gate

If D_n is a positive integer such that D_n A_n and D_n B_n are integers, and the positive expression D_n(A_n L-B_n) tends to zero, then L is irrational. For if L=p/q, its product with q is a positive integer tending to zero, impossible. Neither positivity nor denominator control can be omitted. Bound (3) gives a sufficient condition D_n/((n+1)^3 A_n)->0.

## Computation and next steps

Using Python fractions.Fraction, checked W_n=6/(n+1)^3 for n=0,...,20 and the finite continued-fraction identity for N=1,...,20. Initial A: 1,5,73,1445,33001,819005. Initial B: 0,6,351/4,62531/36,11424695/288,35441662103/36000.

Priority next step: derive A_n=sum_(k=0)^n binom(n,k)^2 binom(n+k,k)^2 from an integral or auxiliary polynomial, and derive its recurrence by a displayed telescoping certificate. That supplies an understandable origin for the coefficients. Independently identify L via a primary-source integral or partial-fraction construction. Read primary papers before attributing these standard identities historically.

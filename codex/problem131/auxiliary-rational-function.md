# Auxiliary rational function: a route toward explaining the recurrence

Status: independently derived algebraic theorem below; connection to the full irrationality proof remains a project with explicit gaps. This is background aligned with the user's topic, not a transcription or solution of the inaccessible official #131 statement. No novelty claim.

For n >= 0 set
R_n(t) = [product_{j=1}^n(t-j) / product_{j=0}^n(t+j)]^2.
Let H_k^(s)=sum_{j=1}^k j^(-s), H_0^(s)=0, H_k=H_k^(1), and
 a_{n,k}=binom(n,k)^2 binom(n+k,k)^2,
 b_{n,k}=2 a_{n,k}(2H_k-H_{n+k}-H_{n-k}).

**Theorem (proved here).**
R_n(t)=sum_{k=0}^n [a_{n,k}/(t+k)^2+b_{n,k}/(t+k)].
Consequently, with A_n=sum_k a_{n,k} and
B_n=sum_k [a_{n,k}H_k^(3)+(b_{n,k}/2)H_k^(2)],
 -1/2 sum_{t=1}^infinity R_n'(t)=A_n zeta(3)-B_n.
The first pairs (A_n,B_n) are (1,0),(5,6),(73,351/4).

**Proof.** At t=-k the coefficient of the double pole is the square of
product_{j=1}^n(-k-j)/product_{0<=j<=n,j!=k}(j-k), giving a_{n,k}.
The logarithmic derivative of (t+k)^2 R_n(t), evaluated at -k, is
2[-(H_{n+k}-H_k)-(H_{n-k}-H_k)], giving b_{n,k}.
Subtracting all principal parts leaves a polynomial that vanishes at infinity,
so it is zero. Since R_n(t)=O(t^-2), comparison of 1/t coefficients gives
sum_k b_{n,k}=0. Differentiate and sum the finite partial fractions over positive
integers. All resulting series converge absolutely. The zeta(2) coefficient
vanishes by sum b=0, leaving the displayed identity. At positive integer zeros
t<=n the derivative is zero, so these terms cause no difficulty. QED.

This makes a useful design principle visible: double poles give zeta(3) after
differentiation; decay at infinity cancels the unwanted zeta(2); integer roots
remove the first n summation terms; squared residue products generate Apéry's
binomial coefficients. The recurrence can be sought from rational-function
identities rather than guessed from numeric approximations.

**Next construction (not yet proved).** Find an explicit rational telescoping
certificate for the n-shift relation with coefficients
(n+1)^3, -(34n^3+51n^2+27n+5), n^3.
One must track boundary terms after differentiation and summation. Merely
checking many n does not prove such an identity.

**Arithmetic gap.** A_n is visibly integral. The needed denominator bound for
B_n is not immediate from the harmonic formula: H_{n+k} introduces denominators
up to 2n. Prove cancellation to an lcm(1,...,n)^3 bound (or a stated adequate
variant) before claiming an irrationality proof. Also supply a rigorous decay
bound for the linear form and prove its nonvanishing. These are separate tasks.

**Concrete generalisation to investigate.** Replace the integer roots 1,...,n
with a chosen root multiset and prescribe double poles at 0,-1,...,-n. Determine
which root patterns retain integral residues and useful small linear forms.
The partial-fraction proof still works for a numerator of degree <=2n; a
recurrence of order two is not guaranteed. This separates universal pole
cancellation from the exceptional arithmetic and recurrence structure.

Primary reference read: R. Apéry, “Irrationalité de ζ2 et ζ3”, Astérisque 61
(1979), pp. 11–13, https://www.numdam.org/item/AST_1979__61__11_0/ .
Page 13 explicitly states the ζ3 recurrence and initial values. The rational
function derivation above was developed independently in this session; it is
not attributed to that short paper and its historical provenance is unresolved.

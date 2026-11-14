# Wilson weight/norm formula retrieval

Retrieved 2026-11-14 through the web tool's open command:
https://dlmf.nist.gov/18.25

Tool result reference: turn750667view0.
This is a transcription of the displayed mathematical formulas, not an
archived raw HTML download. No hash of unretained HTML is claimed.

Equation 18.25.2:

  integral_0^infinity p_n(x) p_m(x) w(x) dx = h_n delta_(n,m).

Equation 18.25.3:

  p_n(x)=W_n(x;a_1,a_2,a_3,a_4).

Equation 18.25.4:

  w(y^2) = (1/(2y)) | product_j Gamma(a_j+iy) / Gamma(2iy) |^2.

Equation 18.25.5:

  h_n = n! * 2*pi * product_(j<ell) Gamma(n+a_j+a_ell)
        / [(2n-1+sum_j a_j) Gamma(n-1+sum_j a_j)].

Table 18.25.1 states x=y^2, y in (0,infinity), with Re(a_j)>0
and nonreal parameters occurring in conjugate pairs.

For a_j=1/2, the six pair sums all equal 1 and sum_j a_j=2.
Thus h_n=2*pi*(n!)^6/(2n+1). Our measure is w(t)dt/(2*pi),
so its Wilson squared norm is (n!)^6/(2n+1).

The DLMF references these classical formulas to Wilson (1980) and,
for this subsection, Andrews et al. (1999), equation (3.8.3).

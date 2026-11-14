# An exact test of the rational Hardy comparison at Fauzan's first parameter set

Status: bounded exact computation plus rigorous Arb interval enclosures. This tests a specific candidate simplification of the analytic estimate, not the arithmetic normalization theorem or its asymptotic regime.

## Parameters and exact matrix

Use K=40,N=3,h=37,r=6, damping order L=51 and damping pole a=41 in `rational-orthogonality.md`. This is the n=1 case of K=40n,N=3n,h=37n,L=52n-1. The choice a=41 avoids coinciding with an existing pole.

The comparison matrix is a rational Hankel matrix with entries

    H[i,j]=(1/π)integral_0^infinity
      y^(2i+2j) D_3(y²)^6/[D_40(y²)(y²+41²)^51]dy.

After cancellation, the rational function in t=y² has denominator

    product_(j=4)^40(t+j²) * (t+41²)^51,

and numerator D_3(t)^5 t^(i+j). It is proper for every required moment i+j<=72.

All 73 moments were computed exactly. At the 37 simple poles the coefficients are ordinary rational residues. At -41², Taylor-expand D_3(t)^5/product_(j=4)^40(t+j²) through order50. The integral of each repeated pole is the exact rational quantity

    (1/π)integral_0^infinity (y²+a²)^(-p)dy
       =binom(2p-2,p-1)/[2*4^(p-1)*a^(2p-1)].

This gives H exactly, and its 37-by-37 determinant was computed with FLINT rational arithmetic. The entire calculation took about three seconds in this environment; no unbounded symbolic search was used.

## Rigorous scalar supremum bound

The closed weight is

    w_5(y)=(2π)^4 y^5/12 * q(1+11q+11q²+q³)/(1-q)^5,
    q=exp(-2πy).

Arb at192-bit precision certifies

    w_5(y)(1+y²/41²)^51 <7/20, y>=0.

The proof in the script splits the line into:

* 0<=y<=1/100: use the analytic bound w_5(y)<=exp(2πy)/π.
* 1/100<=y<=5: 49,900 intervals of length1/10,000; bound each positive monotone factor at the appropriate endpoint. Every comparison is an Arb interval comparison, not a floating-point sample.
* y>=5: use sum ell^4 exp(-2πell y)<=2exp(-2πy), log(1+u²)<=u, and monotonicity of y^5exp(-(2π-51/41)y) past5.

The resulting allowed constant is C=(7/20)41^102. The largest displayed regional upper bound is below0.338953, leaving room beneath0.35. The certificate could be tightened, but the lower bound below shows this would not repair the comparison.

## Result: this particular comparison fails decisively

The exact determinant and rigorous logarithmic enclosures give

    log det H = -11795.14196443...,
    log[(πC)^37 det H] = 2223.39060248....

For comparison with the paper's normalized determinant, include its explicit factor

    S_K=(K!)^(2h) 4^(h-1)
          /[(N!)^(12h) product_(i=1)^(h-1)((2i)!)²].

Then the proved analytic upper bound is

    log[S_K (πC)^37 det H] = 2019.07160734....

It is enormous rather than small. This is an upper bound on the true normalized determinant, so it does **not** prove that the true determinant is large.

There is a stronger diagnosis of this particular bounding method. Since w_5(0)=1/π, the exact optimal scalar constant satisfies

    C_opt >=41^102/π.

Therefore even an exact maximization of the scalar weight would give a comparison quantity with

    log[S_K (πC_opt)^37 det H]
       >=log[S_K 41^(102*37)det H]
       =2015.56002017... >2015.

Thus scalar optimization alone cannot make this fixed-a, fixed-L comparison useful at n=1. The loss is structural: the supremum replaces a varying weight by a uniform constant and discards its relation to the rational orthogonal directions.

This is not an asymptotic no-go theorem, not a test of every damping pole, and not a disproof of the rational comparison strategy. In particular no large-n statement follows from this first instance. The prime-dependent arithmetic factor m_(K,M) was not included: the source's relevant asymptotic integrality theorem has thresholds far above K=40. The calculation is specifically an analytic comparison of S_K times the determinant.

## Reproduction and saved exact data

Run

    PYTHONPATH=/workspace/scratch/odd-zeta-python python codex/computations/rational_comparison_experiment.py

The isolated scratch installation provides python-flint0.9.0. Elsewhere, install that package normally. Exact determinant, S_K and C are saved in `rational_comparison_experiment.exact.txt`; interval results and parameters are in `rational_comparison_experiment.json` in the same computations directory. The exact-data SHA-256 is

    4341d9a66e50b4a2b07871dba6f599ab2efa5dc404d859c5f13ef2c186c1cec4.

The precise lesson is to retain a weighted norm or equilibrium estimate instead of flattening the weight by a global supremum. The rational orthogonality link remains an explanatory theorem, but this elementary scalar comparison is not yet a simplification of Fauzan's successful analytic estimate.

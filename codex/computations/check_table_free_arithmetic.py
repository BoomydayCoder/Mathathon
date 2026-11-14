"""Exact rational certificate for the table-free Fauzan arithmetic bound.

Standard-library only. This checks the final algebra, not the source's local
p-adic lemmas or the prime-number-theorem passage to Proposition 5.2.
"""
from fractions import Fraction as F

alpha = F(3, 40)
lam = 1 - alpha
H5 = sum((F(1, j) for j in range(1, 6)), F(0))

# Positive-part primitives from the two intervals in source (5.8).
outer_left = F(41, 72) - 4 * alpha / 3 + F(5, 4) * (F(1, 3) - alpha)**2 - 5 * alpha**2 / 6
outer_right = F(1, 8) + 3 * (F(1, 2) - alpha)**2 + F(5, 2) * alpha**2
rank_correction = F(5, 2) * alpha**2
floor_contribution = -5 * lam / 3
hinges = sum(((2 * lam - F(j, 3))**2 / (2 * j) for j in range(1, 6)), F(0))
assert hinges == 2 * lam**2 * H5 - 10 * lam / 3 + F(5, 6)
outer = outer_left + outer_right - rank_correction + floor_contribution + hinges
assert outer == F(29, 12) - F(31, 6) * alpha + F(41, 12) * alpha**2 - 5 * lam + 2 * lam**2 * H5
assert outer == F(127751, 96000)
assert rank_correction == F(9, 640)

# log(alpha) = -2 atanh(1/4) - 6 atanh(1/3).
# Omitted terms are strictly negative, so truncation is an upper bound.
log_upper = -2 * (F(1, 4) + F(1, 4)**3 / 3) - 6 * (F(1, 3) + F(1, 3)**3 / 3 + F(1, 3)**5 / 5)
assert log_upper == F(-33559, 12960)
M = 1000
bound = outer + F(13, 144) + lam * (F(5, 3) + 8 * alpha + 12 * alpha * log_upper) + 14 * lam / M
assert bound == F(439979, 320000)
assert F(11, 8) - bound == F(21, 320000)
assert bound < F(11, 8) < F(69, 50)
print('I_out =', outer)
print('log(alpha) <', log_upper)
print('limsup K^-2 log(m_K,1000) <=', bound, '=', float(bound))
print('gap below 11/8 =', F(11, 8) - bound)
print('with analytic bound -69/50, gap below zero >', F(69, 50) - F(11, 8))

# Optional simplification of the outer local lemma: discard its extra rank
# saving, hence add back the integral of d. Local justification is separate.
unrefined = bound - 14 * lam / 1000 + 14 * lam / 10000 + rank_correction
assert unrefined == F(2203747, 1600000)
assert F(69, 50) - unrefined == F(4253, 1600000) > F(1, 400)
print('unrefined outer rank, M=10000:', unrefined, '=', float(unrefined))

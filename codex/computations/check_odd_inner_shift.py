"""Exact arithmetic for the general odd-index envelope obstruction.

Checks scalar algebra and class-weight tables; does not audit local p-adic
completion or prove sharpness of the sufficient valuation envelopes.
"""
from fractions import Fraction as F

alpha = F(3, 40)
lam = 1 - alpha

def class_cost(s, ell, removed):
    return sum(max(0, ell + s - 1 - 6 * removed - 2 * i)
               for i in range(ell - 2))

assert [class_cost(5, ell, 1) for ell in range(2, 7)] == [0, 1, 2, 4, 6]
assert [class_cost(7, ell, 1) for ell in range(2, 7)] == [0, 3, 6, 9, 12]
for s in (5, 7, 9, 11):
    for ell in range(2, 7):
        assert class_cost(s, ell, 0) == (s + 2) * (ell - 2)

log_lower = -2 * F(1, 4) / (1 - F(1, 4)**2) - 6 * F(1, 3) / (1 - F(1, 3)**2)
assert log_lower == F(-167, 60)
J_lower = lam * (F(4, 3) + 10 * alpha + 12 * alpha * log_lower)
assert J_lower == F(-9361, 24000)
rank_saving = (1 + 4 * alpha - F(2, 3))**2 / 4
assert rank_saving == F(361, 3600)
outer_plain = F(129101, 96000)
lower = outer_plain + J_lower - F(1, 36) - F(13, 16 * 40**2) + 2 * lam / 3 - rank_saving
assert lower == F(1662179, 1152000)
assert lower > F(7, 5) > F(69, 50)
print('s=7 extra small+inner cost:', 2 * lam / 3)
print('maximum saving from deleting old plain rank penalty:', rank_saving)
print('direct arithmetic envelope >', lower, '=', float(lower))

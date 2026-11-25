import exact_enum, sys
# p=17, mu=4: exact PSD test only for |s|>=6 (the decisive tail)
exact_enum.run(17,4,only_abs_s_ge=6)

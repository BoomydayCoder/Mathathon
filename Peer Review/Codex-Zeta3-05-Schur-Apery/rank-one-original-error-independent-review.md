# Independent review: original error tail, arithmetic, and sharp rate

Verdict: accepted. Root audited rank-one-original-error-tail-and-rate.md.
The original quotient identification and the separately reviewed Apéry
asymptotics are the explicit inputs; no finite extrapolation is used.

For V_(h+1)=a_h V_h-c_h V_(h-1), the oriented Casoratian
Y_h L_(h-1)-Y_(h-1)L_h multiplies by c_h. The initial value -23/96 and
the factorial ratio give exactly -d_h/[2(h-1)(2h)!(2h-2)!]. Dividing by
the two L factors yields the increment with coefficient72 and denominator
h²(h-1)³C_hC_(h-1). Positivity follows since C_h,L_h,d_h are positive.
The established limit q_h->zeta(3) justifies the positive infinite tail.

The classical Casoratian A_(h-1)B_h-A_hB_(h-1)=6/h³ gives the exact bracket
width6/[h³A_hC_h]; h=1 supplies width1/20 and fixes normalization. The ratio
of the reviewed asymptotic constants is C_f/C_A=4pi³/lambda. Since
f_(h-1)/f_h->lambda and A_(h-1)/A_h->lambda^(-1), the original upper error
is asymptotic to the lower one times (lambda-5)/(5-lambda^(-1))=sqrt(lambda).
This coefficient is nonzero and positive. Hence the full leading constant
is4pi³/sqrt(lambda)=4pi³(3-2sqrt2), with rate lambda^(-2h) and no power ofh.

The existing measure-comparison theorem sandwiches every fixed higher-odd
rank-one optimized error between this one times polynomial factors. Taking
h-th roots therefore proves the same actual rate limit for each fixed odd
index. The result is not uniform in an index growing withh and provides no
new denominator estimate. This closes the previously unresolved analytic
rate, not the higher-odd irrationality problem.

Finally the lcm-cube theorem for B_h also clears5B_h-B_(h-1), while
5A_h-A_(h-1) is integral. Its positive numerator error is less than f_(h-1).
The established elementary lcm rate and Apéry decay make the cleared error
tend to zero. Thus the arithmetic transfer to the original Schur family is
correct and is explicitly dependent on the classical denominator theorem.

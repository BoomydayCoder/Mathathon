# Independent limited audit of the Apéry CM-character bridge

Status: checked the elementary arithmetic/Hasse bridge in Sections 2–3 and the concluding implication in Section 6 of apery-cm-character.md. No gap found in those parts. This review does not independently certify every modular-function pole bound or eta-product identity in Section 4; those require the separate modular audit. Standard elliptic specialization/Cartier facts remain explicitly identified machinery.

1. For b_n=sum_k binom(n,k)^2 binom(2k,k), write n=pN+r,k=pK+s. Lucas splits the first binomial. If 2s>=p, the low digit of 2k is 2s-p<s, so the central binomial vanishes modulo p. Otherwise it splits into binom(2K,K)binom(2s,s). Zero extension handles s>r and summing gives b_{pN+r}=b_N b_r modulo p.

2. Assuming the previously proved characteristic-zero pullback A(T(x))=(1-x)B(x)^2, division by its x^p version gives P_p(T(x))=(1-x)/(1-x^p) B_p(x)^2=(1-x)^{1-p}B_p(x)^2. This verifies the power and direction of the factor in equation (3), with no numerical extrapolation.

3. For Lambda=(1+X+Y)(1+X^{-1}+Y^{-1}), the constant term of Lambda^n is the sum of squared multinomial coefficients. Grouping by k=b+c and using sum_b binom(k,b)^2=binom(2k,k) yields exactly b_n. Expanding (1-x Lambda)^{p-1}, the coefficient (-1)^n binom(p-1,n) is one modulo p, proving the finite constant-term formula.

4. With 1-x Lambda=a-bY-c/Y, its Y-constant term has summands binom(p-1,2j)binom(2j,j)a^{p-1-2j}(bc)^j. Since binom(p-1,2j)=1 modulo p and binom((p-1)/2,j)(-4)^j=binom(2j,j) modulo p, this equals (a^2-4bc)^{(p-1)/2}. Multiplication by X^{p-1} converts the X-constant term to coefficient X^{p-1} of the stated quartic power. All signs and degrees agree.

5. For V^2=Phi(X), omega=dX/V=Phi^{(p-1)/2}dX/V^p. The numerator degree is 2p-2, so Cartier retains only its X^{p-1} coefficient. Consequently Cartier(omega)=B_p(x)^{1/p}omega over an algebraic closure. Given the standard genus-one Cartier/supersingularity criterion, this proves the exact vanishing equivalence in Section 3.

6. The critical polynomial D_x=1-18x+9x^2 satisfies Delta(T(x))=D_x^2/(1-x)^2. For p>=5 its roots are distinct and avoid 0,1,1/9. In particular the pullback factor (1-x)^{p-1} is nonzero at either critical point. Therefore P_p(T(x_0)) vanishes exactly when B_p(x_0) does. Given the CM identification and reduction criterion in Sections 4–5, the final equivalence to (-6/p)=-1 follows, and the earlier square-factor theorem identifies epsilon correctly.

The direct CM argument in Section 5 was also sanity-checked: for inert p the CM generator acting on ordinary geometric p-torsion would solve a^2=-6 in F_p, impossible. For split p a principal power beta of one prime ideal has p-power norm and trace prime to p. In a supersingular reduction its endomorphism and its conjugate would both be purely inseparable, hence have zero tangent map, contradicting their sum being multiplication by that trace. This is sound assuming the explicitly stated good-reduction specialization preserves endomorphism relations and degree, and the standard description of geometric p-torsion.

No historical novelty assertion. This review should be read as an audit of specified steps, not a replacement for verification of all modular/CM inputs.

# Riemann-Explicit-Formula-for-Primes

Author: Daniel Hutama, McGill University Desautels Faculty of Management

The form of Riemann's R(x) formula for primes using the nontrivial zeros of the zeta function assumes the Riemann Hypothesis and uses the approximation li(x)~(x/logx) * (asymptotic term),
where the asymptotic term is a divergent infinite sum. 
The code allows the user to determine the cutoff term of the infinite sum, depending on desired precision.


There are 3 infinite sums involved in this form of the explicit formula.

The first sum over n affects the "main term" and the correction term involving the zeta zeros.

The second sum over the (positive) imaginary parts of the nontrivial zeros of the zeta function (denoted gamma) affets the "primary correction term."

The third sum is over k in the asymptotic approximation li(x)~x/log(x)*sum_k(k!/(log(x))^k). For very large x, the cutoff term is significant as the error is on the order of sqrt(x)/log(x)^k.

Note: n starts at 1, gamma starts at gamma_1 = 14.1347..., and k starts at 0.

The user should call Riemann_explicit(x, numbersums, numberzeros, numberasymterms) to obtain an approxiation for the number of primes up to x, 
where "x" is a real number > 2. "numbersums", "numberzeros", and "numberasymterms" refer to the cutoff terms for n, gamma, and k respectively.

Special note about k (numberasymterms): If you enter 4 for numberasymterms, the term for k=1 (i.e. 1) is already included in the code, so entering 4 will cause the code to sum over k in the range [0,4] inclusive.
This is not the case for n (numbersums) or gamma (numberzeros). If you enter 4 for numbersums, the code will sum over n in the range [1,4] inclusive.
If you enter 4 for numberzeros, the code will sum over the first 4 zeta zeros. I.e. gamma_1 up to gamma_4, inclusive.

# Riemann Explicit Formula for Primes up to x

Author: Daniel Hutama, McGill University Desautels Faculty of Management

This form of Riemann's R(x) formula for primes using the nontrivial zeros of the zeta function assumes the Riemann Hypothesis.
The code allows the user to determine the cutoff terms of each infinite sum, depending on desired precision.


There are 2 infinite sums involved in this form of the explicit formula.

The first sum over n affects the "main term" and the correction term involving the zeta zeros.

The second sum over the (positive) imaginary parts of the nontrivial zeros of the zeta function (denoted gamma) affects the "primary correction term."


Note: n starts at 1, and gamma starts at gamma_1 = 14.1347... .

The user should call Riemann_explicit(x, numbersums, numberzeros) to obtain an approxiation for the number of primes up to x, 
where "x" is a real number > 2. The inputs "numbersums" and "numberzeros" refer to the cutoff terms for n, and gamma respectively.

Note: If you enter 4 for numbersums, the code will sum over n in the range [1,4] inclusive.
If you enter 4 for numberzeros, the code will sum over the first 4 zeta zeros. I.e. gamma_1 up to gamma_4, inclusive.

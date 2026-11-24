# Signed autocorrelation constant C_1d: upper bound 1.4575790574398573

Certificate for an upper bound on the signed constant C_1d of Tao's optimization-constants table (https://teorth.github.io/optimizationproblems/constants/1d.html).

- `c4800.txt`: heights v_0..v_4799 of a step function f on [-1/4, 1/4] with 4800 equal cells (h = 1/9600), one decimal per line, read exactly.
- Ratio max_m |(f*f)(m h)| / (int f)^2 = 1749094868927836737655798170753876661351000000 / 1200000000000005456633880000006203094437575803 = 1.4575790574398573... (exact rational arithmetic; f*f is piecewise linear, so its max is at a grid point).
- `checker.py`: exact Fraction check (`python3 checker.py c4800.txt`, add `--binary64` to read heights as doubles; same verdict). `indep.py`: independent big-integer convolution written from the problem statement alone; agrees.

Found by Alexander LaVerdi (Independent), built with an AI research agent I developed (LSE-annealed FFT L-BFGS, refined 600 -> 1200 -> 2400 -> 4800 cells). MIT license.

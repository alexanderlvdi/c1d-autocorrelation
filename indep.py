#!/usr/bin/env python3
"""Independent verifier-5 checker (written from the statement, not from checker.py).

Statement: for step function f = sum_i v_i * 1_[i h, (i+1) h), h = 1/(2N), on [-1/4,1/4],
  ratio = max_t |(f*f)(t)| / (int f)^2.
(f*f) is piecewise linear with knots at multiples of h, so the max is at knots:
  (f*f)(m h) = h * sum_{i+j=m-1} v_i v_j   (m = 1 .. 2N-1).
  int f = h * S, S = sum v_i.
  ratio = h*M / (h*S)^2 with M = max_m |sum_{i+j=m-1} v_i v_j|  =  2N * M / S^2.
All arithmetic is exact INTEGER arithmetic: every heights is scaled to an integer by a common
power-of-ten (decimal) or power-of-two (binary64) denominator D; the scale cancels in the ratio:
  ratio = 2N * max|conv(w,w)| / (sum w)^2,  w_i = v_i * D integers.
The convolution is done with Python big ints via Kronecker substitution (independent of the
Fraction loop in checker.py).

Usage: indep.py <file> [--binary64]     (file: whitespace/comma separated numbers or JSON {"values":[...]})
Prints N, ratio (as float and exact integer fraction), gate, and VERDICT.
"""
import sys, re, json
from fractions import Fraction
from math import lcm

GATE = Fraction("1.457854884921989")  # decimal gate from the notebook, fixed

def read_values(path):
    txt = open(path).read()
    if txt.lstrip().startswith("{"):
        m = re.search(r'"values"\s*:\s*\[(.*?)\]', txt, re.S)
        if m is None:
            raise SystemExit("JSON input without a \"values\" array")
        body = m.group(1)
    else:
        body = txt
    return re.findall(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?', body)

def to_fraction(tok, binary64):
    if binary64:
        return Fraction(float(tok))          # exact value of the binary64 nearest to tok
    return Fraction(tok)                     # exact decimal

def integerize(vals):
    """Scale rationals to integers with a common denominator D (exact)."""
    D = 1
    for x in vals:
        D = lcm(D, x.denominator)
    return [int(x * D) for x in vals], D

def conv_max_abs(w):
    """max_k |sum_{i+j=k} w_i w_j| with exact big-int Kronecker substitution.
    Splits w into positive and negative parts so every packed coefficient is >= 0."""
    n = len(w)
    pos = [x if x > 0 else 0 for x in w]
    neg = [-x if x < 0 else 0 for x in w]
    B = max(max(abs(x) for x in w), 1)
    bits = 2 * B.bit_length() + n.bit_length() + 2          # enough for any coefficient
    def pack(a):
        acc = 0
        for x in reversed(a):
            acc = (acc << bits) | x
        return acc
    P, Q = pack(pos), pack(neg)
    PP, QQ, PQ = P * P, Q * Q, P * Q
    # conv = PP + QQ - 2*PQ  (since conv(w,w) = conv(p,p)+conv(q,q)-2 conv(p,q)), coefficientwise.
    # Each product is packed, so unpack coefficient by coefficient with an offset-free method:
    # compute the signed polynomial directly as integers and unpack with a signed digit walk.
    total = PP + QQ - 2 * PQ
    mask = (1 << bits) - 1
    coeffs = []
    carry = 0
    x = total
    # Signed base-2^bits digit extraction (digits in (-2^(bits-1), 2^(bits-1)]).
    half = 1 << (bits - 1)
    for _ in range(2 * n - 1):
        d = x & mask
        x >>= bits
        if d >= half:
            d -= (1 << bits)
            x += 1
        coeffs.append(d)
    assert x == 0 or x == -1, "unpack did not consume the whole polynomial"
    return max(abs(c) for c in coeffs), coeffs

def ratio_exact(vals):
    w, _ = integerize(vals)
    n = len(w)
    S = sum(w)
    if S == 0:
        return None
    M, _ = conv_max_abs(w)
    return Fraction(2 * n * M, S * S)

if __name__ == "__main__":
    path = sys.argv[1]
    binary64 = "--binary64" in sys.argv
    toks = read_values(path)
    vals = [to_fraction(t, binary64) for t in toks]
    N = len(vals)
    r = ratio_exact(vals)
    if r is None:
        print(f"N={N} zero integral"); print("VERDICT: FAIL"); sys.exit(1)
    print(f"N={N} mode={'binary64' if binary64 else 'decimal'} ratio={float(r):.16f} exact={r.numerator}/{r.denominator}")
    print(f"gate={float(GATE):.16f}")
    print("VERDICT: " + ("PASS" if r < GATE else "FAIL"))

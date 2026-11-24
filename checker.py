#!/usr/bin/env python3
"""Frozen checker: C1d signed autoconvolution step-function upper bound.
Usage: checker.py <file> [--binary64]   file = JSON ({"values":[...]}) or whitespace/comma separated numbers.
Heights read EXACTLY as decimals (default) or as exact binary64 values (--binary64).
ratio = max_m |(f*f)(m h)| / (h*sum v)^2, h = 1/(2N), (f*f)(m h) = h*sum_{i+j=m-1} v_i v_j.
PASS iff ratio < GATE (exact Fraction of 1.457854884921989) and sum != 0."""
import sys, json, re
from fractions import Fraction
GATE = Fraction("1.457854884921989")  # previous best listed bound (SCiEL2026)

def load(path, b64):
    txt = open(path).read()
    if txt.lstrip().startswith("{"):
        toks = re.findall(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?', re.search(r'"values"\s*:\s*\[(.*?)\]', txt, re.S).group(1))
    else:
        toks = re.findall(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?', txt)
    return [Fraction(float(t)) if b64 else Fraction(t) for t in toks]

def ratio(v):
    n = len(v); h = Fraction(1, 2 * n)
    s = sum(v)
    if s == 0: return None
    best = Fraction(0)
    for m in range(1, 2 * n):                       # knots; endpoints m=0,2n give 0
        c = sum(v[i] * v[m - 1 - i] for i in range(max(0, m - n), min(n, m)))
        c = abs(c * h)
        if c > best: best = c
    return best / (h * s) ** 2

if __name__ == "__main__":
    v = load(sys.argv[1], "--binary64" in sys.argv)
    r = ratio(v)
    if r is None: print("N=%d zero integral" % len(v)); print("VERDICT: FAIL"); sys.exit(1)
    print("N=%d ratio=%.16f exact=%d/%d" % (len(v), r, r.numerator, r.denominator))
    print("gate=%.16f" % float(GATE)); print("VERDICT: " + ("PASS" if r < GATE else "FAIL"))

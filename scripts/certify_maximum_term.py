#!/usr/bin/env python3
"""Global, interval-certified bound in the He--Tang entire-function family.

The squared trigonometric modulus is a polynomial in x=cos(theta)^2.
Bernstein subdivision covers the whole interval [0,1]; no sampled maximum
or unproved list of stationary points is accepted as evidence.
"""
import argparse
import json
from fractions import Fraction as Q
from math import comb
from pathlib import Path

from mpmath import iv
from certify_turan import rational, interval, endpoints, show

if not __debug__:
    raise RuntimeError("run without -O: assertions are part of certificate verification")


def odd_chebyshev(nmax):
    """Exact coefficients of cos((2n+1)t)/cos(t) in x=cos(t)^2."""
    polys = [[1]]
    if nmax:
        polys.append([-3, 4])
    for n in range(1, nmax):
        a, b = polys[-1], polys[-2]
        new = [0]*(len(a)+1)
        for k, v in enumerate(a):
            new[k] -= 2*v
            new[k+1] += 4*v
        for k, v in enumerate(b):
            new[k] -= v
        polys.append(new)
    return polys


def coefficients(kq, aq, nmax):
    k, a = interval(kq), interval(aq)
    return [iv.mpc(iv.cos(a*t), iv.sin(a*t))/k**t
            for t in (n*(n+1)//2 for n in range(nmax+1))]


def tail(k, nmax):
    # Beyond N, successive triangular exponents differ by at least N+2.
    return 2*k**(-((nmax+1)*(nmax+2)//2))/(1-k**(-(nmax+2)))


def squared_modulus_polynomial(cs):
    r = [iv.mpc(0) for _ in cs]
    for c, poly in zip(cs, odd_chebyshev(len(cs)-1)):
        for k, v in enumerate(poly):
            r[k] += c*v
    q = [iv.mpf(0) for _ in range(2*len(cs))]
    for i, a in enumerate(r):
        for j, b in enumerate(r):
            q[i+j+1] += 4*(a.real*b.real+a.imag*b.imag)
    return q


def bernstein_upper_coefficients(power_coefficients):
    d = len(power_coefficients)-1
    return [endpoints(sum((power_coefficients[j]*interval(Q(comb(k, j), comb(d, j)))
                          for j in range(k+1)), iv.mpf(0)))[1]
            for k in range(d+1)]


def split_bernstein(bs):
    # Exact de Casteljau subdivision at the midpoint.
    left, right = [bs[0]], [bs[-1]]
    while len(bs) > 1:
        bs = [(a+b)/2 for a, b in zip(bs, bs[1:])]
        left.append(bs[0])
        right.append(bs[-1])
    return left, list(reversed(right))


def global_polynomial_upper(bs, threshold, max_depth=70):
    stack = [(bs, 0)]
    nodes, leaves, deepest = 0, 0, 0
    upper = Q(0)
    while stack:
        values, depth = stack.pop()
        nodes += 1
        bound = max(values)
        if bound <= threshold:
            upper = max(upper, bound)
            leaves += 1
            deepest = max(deepest, depth)
        else:
            assert depth < max_depth and nodes < 200_000, "bound not proved by complete subdivision"
            a, b = split_bernstein(values)
            stack.extend([(a, depth+1), (b, depth+1)])
    assert leaves > 0
    return upper, {"nodes": nodes, "leaves": leaves, "maximum_depth": deepest}


def certify(data):
    k, a = rational(data["K"]), rational(data["alpha"])
    claim = rational(data["lower_bound"])
    nmax = data["N"]
    assert k > 1 and isinstance(nmax, int) and 1 <= nmax <= 30
    assert Q("0.585078819653") < claim < Q("0.6366")
    cs = coefficients(k, a, nmax)
    t = tail(k, nmax)
    available = 1/claim-t
    assert available > 0
    q = squared_modulus_polynomial(cs)
    # Early rejection if the endpoint already contradicts the target.
    assert endpoints(sum(q, iv.mpf(0)))[0] < available**2
    bs = bernstein_upper_coefficients(q)
    upper, stats = global_polynomial_upper(bs, available**2)
    maxmod = iv.sqrt(interval(upper))+interval(t)
    beta = 1/maxmod
    assert endpoints(beta)[0] >= claim

    # An upper bound on the ACTUAL beta of the old function, obtained from
    # a lower bound for its modulus at theta=0. This proves a change of
    # construction, rather than extra displayed digits of the old witness.
    old_k, old_a = Q("3.568182317714"), Q("3.961543335688")
    old_p = 2*sum(coefficients(old_k, old_a, nmax), iv.mpc(0))
    old_modulus = iv.sqrt(old_p.real**2+old_p.imag**2)-interval(tail(old_k, nmax))
    assert endpoints(old_modulus)[0] > 0
    old_beta_upper = endpoints(1/old_modulus)[1]
    gain = endpoints(beta)[0]-old_beta_upper
    assert gain > 0
    headline_gain = claim-old_beta_upper
    assert headline_gain > 0, "headline must beat the actual previous witness"
    return {
        "certified_lower_bound": str(claim),
        "decimal_lower_bound": data["lower_bound"],
        "K": data["K"], "alpha": data["alpha"], "N": nmax,
        "interval_dps": iv.dps,
        "polynomial_degree": len(q)-1,
        "subdivision": stats,
        "tail_bound": show(interval(t)),
        "certified_maximum_upper_enclosure": show(maxmod),
        "reciprocal_of_maximum_upper_enclosure": show(beta),
        "old_function_beta_upper_enclosure": show(interval(old_beta_upper)),
        "certified_gain_over_old_function_at_least": show(interval(gain)),
        "headline_gain_over_old_function_at_least": show(interval(headline_gain)),
        "status": "PASS: global maximum bounded on the entire circle, with infinite tail included",
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificate", nargs="?", default="certificates/maximum_term.json")
    ap.add_argument("--dps", type=int, default=100)
    ap.add_argument("--output")
    args = ap.parse_args()
    assert args.dps >= 40
    iv.dps = args.dps
    result = certify(json.loads(Path(args.certificate).read_text()))
    rendered = json.dumps(result, indent=2)+"\n"
    print(rendered, end="")
    if args.output:
        Path(args.output).write_text(rendered)

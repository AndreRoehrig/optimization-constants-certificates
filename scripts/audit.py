#!/usr/bin/env python3
"""Additional checks of formulas, finite constructions, and rejection paths.

Direct quadrature here is an independent numerical check, not the rigorous
certificate. The two certify_* scripts supply the interval certificates.
"""
import copy
import json
from fractions import Fraction as Q
from math import comb
from pathlib import Path

import mpmath as mp
from mpmath import iv
import certify_turan as turan
import certify_maximum_term as maximum


def rejects(label, f):
    try:
        f()
    except (AssertionError, ValueError):
        print("PASS rejected:", label)
        return
    raise AssertionError("accepted invalid case: "+label)


def real(q):
    q = Q(q)
    return mp.mpf(q.numerator)/q.denominator


def midpoint(x):
    a, b = turan.endpoints(x)
    return real((a+b)/2)


def quad_integral(alpha, x):
    # u=x exp(-v) removes the integrable singularity at u=0.
    return x**alpha*mp.quad(lambda v: mp.exp(-alpha*v)/(1-x*mp.exp(-v)),
                            [0, 1, 4, 16, 64, 256, mp.inf])


def quad_pair(alpha, tq, rq):
    if tq+rq >= 1:
        return mp.mpc(0)
    t, r, l = real(tq), real(rq), real(1-tq-rq)
    def f(v):
        z = l*mp.exp(-v)
        return mp.exp(-alpha*v)/(1-z)*mp.log((1-t-z)*(1-r-z)/(t*r))
    return l**alpha*mp.quad(f, [0, 1, 4, 16, 64, 256, mp.inf])


def audit_turan(data):
    alpha = mp.mpc(*map(real, data["alpha"]))
    ai = iv.mpc(*map(turan.interval, map(Q, data["alpha"])))
    cuts = list(map(Q, data["breaks"]))
    errors = []
    for x in sorted(set(cuts+[1-t for t in cuts])):
        q = quad_integral(alpha, real(x))
        enclosed = turan.integral(ai, x, 480)
        errors.append(abs(q-mp.mpc(midpoint(enclosed.real), midpoint(enclosed.imag))))
    for i, t in enumerate(cuts[:-1]):
        for r in cuts[:i+1]:
            q = quad_pair(alpha, t, r)
            enclosed = turan.pair_integral(ai, t, r, 480)
            errors.append(abs(q-mp.mpc(midpoint(enclosed.real), midpoint(enclosed.imag))))
    assert max(errors) < mp.mpf("1e-45")
    print("PASS Turan series vs independent integral quadrature; max error", mp.nstr(max(errors), 5))

    # Compute the finite coefficient cancellation using the exact formula
    # for Y_n, then use a DIFFERENT triangular recurrence for b_n.
    s = 1-alpha
    eta = [mp.mpc(*map(real, z)) for z in data["eta"]]
    jumps = [eta[0]-s]+[b-a for a, b in zip(eta, eta[1:])]
    for n in [30, 31, 32, 72]:
        a = n//3
        middle = list(range(a+1, n-a))
        final = list(range(n-a, n+1))
        beta = [mp.mpc(1)]
        for l in range(1, n+1):
            beta.append(beta[-1]*(alpha+l-1)/l)
        w = {m: sum((d for d, t in zip(jumps, cuts) if m > n*t), mp.mpc(0))
             for m in middle}
        p = {m: real(Q(n, m))*beta[n-m] for m in final}
        y = n*beta[n]-n*sum((w[m]*beta[n-m]/m for m in middle), mp.mpc(0))
        y += real(Q(n, 2))*sum((w[u]*w[v]*beta[n-u-v]/(u*v)
                               for u in middle for v in middle if u+v <= n), mp.mpc(0))
        y += s*sum(p.values(), mp.mpc(0))
        denom = sum(map(abs, p.values()))
        powers = [mp.mpc(0)]+[s]*a+[s+w[m] for m in middle]
        powers += [y/denom*mp.conj(p[m])/abs(p[m]) for m in final]
        assert len(powers) == n+1
        b = [mp.mpc(1)]
        for l in range(1, n+1):
            b.append(sum((b[j]*(1-powers[l-j]) for j in range(l)), mp.mpc(0))/l)
        assert abs(b[n]) < mp.mpf("1e-45")
    print("PASS finite coefficient cancellation for n=30,31,32,72 (all floor residue classes)")


def bernstein(bs, x):
    d = len(bs)-1
    return sum(b*comb(d, k)*x**k*(1-x)**(d-k) for k, b in enumerate(bs))


def audit_maximum(data):
    bs = [Q(1, 7), Q(-4, 5), Q(13, 3), Q(2), Q(-1, 2)]
    a, b = maximum.split_bernstein(bs)
    for x in [Q(0), Q(1, 7), Q(3, 5), Q(1)]:
        assert bernstein(a, x) == bernstein(bs, x/2)
        assert bernstein(b, x) == bernstein(bs, (1+x)/2)
    _, stats = maximum.global_polynomial_upper([Q(0), Q(2), Q(0)], Q(1))
    assert stats["leaves"] == 2
    rejects("interior maximum above target", lambda: maximum.global_polynomial_upper(
        [Q(0), Q(2), Q(0)], Q(1, 2), max_depth=8))
    print("PASS exact Bernstein subdivision identities and known interior maximum")

    k, alpha = real(data["K"]), real(data["alpha"])
    nmax = data["N"]
    polynomials = maximum.odd_chebyshev(nmax)
    cs = [k**(-t)*mp.exp(1j*alpha*t) for t in (n*(n+1)//2 for n in range(nmax+1))]
    for i in range(33):
        theta = mp.pi*i/32
        x = mp.cos(theta)**2
        trig = 2*sum(c*mp.cos((2*n+1)*theta) for n, c in enumerate(cs))
        poly = 2*mp.cos(theta)*sum(c*sum(v*x**j for j, v in enumerate(r))
                                 for c, r in zip(cs, polynomials))
        assert abs(trig-poly) < mp.mpf("1e-45")
    print("PASS polynomial reduction vs direct trigonometric evaluation")


def main():
    iv.dps = 100
    mp.mp.dps = 60
    t = json.loads(Path("certificates/turan.json").read_text())
    m = json.loads(Path("certificates/maximum_term.json").read_text())
    turan.certify(t)
    maximum.certify(m)
    print("PASS both interval certificates")
    audit_turan(t)
    audit_maximum(m)
    bad = copy.deepcopy(t); bad["eta"] = [["0", "0"] for _ in bad["eta"]]
    rejects("Turan profile with admissible radii but failed limiting inequality", lambda: turan.certify(bad))
    bad = copy.deepcopy(t); bad["C"] = "0.6889"
    rejects("Turan radius below profile values", lambda: turan.certify(bad))
    bad = copy.deepcopy(t); bad["breaks"][0] = "1/4"
    rejects("Turan correction threshold requiring cubic terms", lambda: turan.certify(bad))
    bad = copy.deepcopy(m); bad["K"] = "3.568182317714"; bad["alpha"] = "3.961543335688"
    rejects("previous maximum-term parameters cannot support the new bound", lambda: maximum.certify(bad))
    print("ALL AUDITS PASSED")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("run without -O: assertions are part of this audit")
    main()

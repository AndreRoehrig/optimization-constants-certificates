# A small parameter improvement for the Erdős maximum-term constant

Author: Andre Roehrig. September 9, 2026.

We propose the computer-assisted bound

$$
 C_{51}\ge 0.5850788196744.
$$

This is a **small** improvement. Nat Sothanaphan's
[February 27, 2026 note](https://drive.google.com/file/d/1wZnzui_eeBE32HnkrnSB7YhfcTOiYolp/view)
certifies $0.585078819653$, although the
[constant page](https://teorth.github.io/optimizationproblems/constants/51a.html)
only displays $0.5850788$. The new certificate also beats an upper bound
on the actual value of Sothanaphan's specified entire function, so the gain
is a change of parameters, not simply additional displayed digits.

## Construction and mathematical dependency

We use the family and Theorem 2.8 of Yixin He and Quanyu Tang,
[*Generalizing the Clunie–Hayman construction in an Erdős maximum-term
problem*](https://arxiv.org/abs/2602.12217):

$$
 f_{K,\epsilon}(z)=\sum_{n=0}^\infty
 \epsilon^{n(n-1)/2}K^{-n(n+1)/2}z^n,
 \qquad \epsilon=e^{i\alpha},\quad K>1.
$$

This is transcendental entire: all coefficients are nonzero and their
$n$-th roots in modulus tend to zero. The cited theorem gives

$$
 \beta(f_{K,\epsilon})=\frac1{\max_\theta|P(\theta)|},\qquad
 P(\theta)=2\sum_{n=0}^\infty K^{-T_n}e^{i\alpha T_n}
                         \cos((2n+1)\theta),\qquad T_n=n(n+1)/2.
$$

For context, the theorem follows from the scaling identity of the bilateral
Laurent series, its negligible negative-index tail as $|z|\to\infty$,
and concavity of $\log(\mu/M)$ between successive radii where adjacent
terms tie. We rely on that stated theorem; this computation does not
establish a new general reduction theorem.

Use the following exact rational parameters (terminating decimals):

$$
 K=3.5681854310698435,\qquad
 \alpha=3.9615431126097391.
$$

These exact parameters specify the entire function to be verified. The
checker proves a bound on the whole circle directly from them.

## A certificate for the maximum on the whole circle

Truncate at $N=9$. The infinite tail is bounded uniformly by the exact
rational number

$$
 T=\frac{2K^{-(N+1)(N+2)/2}}{1-K^{-(N+2)}}.
$$

Set $x=\cos^2\theta$. Define integer polynomials

$$
 R_0(x)=1,\quad R_1(x)=4x-3,\quad
 R_{n+1}(x)=(4x-2)R_n(x)-R_{n-1}(x).
$$

The cosine addition identity gives
$\cos((2n+1)\theta)=\cos\theta\,R_n(\cos^2\theta)$; the formula also
holds at $\cos\theta=0$ by continuity. With
$c_n=K^{-T_n}e^{i\alpha T_n}$, form the real polynomial of degree 19

$$
 Q(x)=4x\left|\sum_{n=0}^9c_nR_n(x)\right|^2.
$$

For every real $\theta$, $|P_9(\theta)|^2=Q(\cos^2\theta)$.
Thus an upper bound for $Q$ on the *entire* interval $[0,1]$ suffices.

If $Q(x)=\sum_{j=0}^d q_jx^j$, its Bernstein coefficients on $[0,1]$ are

$$
 b_k=\sum_{j=0}^k q_j\frac{\binom{k}{j}}{\binom{d}{j}},\qquad
 Q(x)=\sum_{k=0}^d b_k\binom{d}{k}x^k(1-x)^{d-k}.
$$

The basis functions are nonnegative and sum to one. Hence replacing every
$b_k$ by a rigorous rational upper endpoint gives a polynomial which
bounds $Q$ from above, and the largest control coefficient bounds that
polynomial. De Casteljau subdivision at the midpoint is implemented by
exact rational averages, preserving these inequalities.

The checker starts with $[0,1]$, accepts a cell only when all its upper
control coefficients meet the target bound, and otherwise splits it into
its two halves. Every cell must be accepted; exhaustion of a depth or node
budget fails the verification. The successful run covers $[0,1]$ by 19
leaves, with 37 total nodes and maximum depth 18.

After adding the infinite tail, it proves

$$
 \max_\theta |P(\theta)|
 \le 1.709171425067816128872617735945.
$$

The reciprocal of the computed upper bound is enclosed by

$$
 [0.585078819674464320841496314106,
  0.585078819674464320841496314107].
$$

The rational headline $0.5850788196744$ is below the lower endpoint.
Both 100- and 200-decimal-place interval calculations pass.

## Comparison with the actual previous witness

The previous exact parameters are

$$
 K_{\rm old}=3.568182317714,\qquad
 \alpha_{\rm old}=3.961543335688.
$$

Evaluate $P_{9,\rm old}(0)$ with intervals, subtract the old uniform tail
bound, and take a reciprocal. Since a global maximum is at least its value
at any specified point, this proves

$$
 \beta(f_{\rm old})
 <0.585078819653352950854632621542.
$$

The new certified lower endpoint exceeds this upper bound by more than

$$
 0.000000000021111369986863692564.
$$

Even the shorter claimed bound $0.5850788196744$ exceeds that upper bound
by more than $2.1047\times10^{-11}$. No assertion that either parameter
pair is optimal is needed or made.

## Replay and status

```sh
python3 scripts/certify_maximum_term.py
python3 scripts/certify_maximum_term.py --dps 200
```

Inputs: `certificates/maximum_term.json`. Transcript:
`certificates/maximum_term_verification.json`. All certificate inputs are
exact rational strings. Only transcendental evaluations use intervals;
control-coefficient subdivision and final comparisons use exact rationals.

This result, its verification code, and the exposition were developed with
AI assistance using Codex. Andre Roehrig reviewed the mathematical arguments,
references, and supporting verification material before publication.
The verification commands were executed and checked during preparation.
The result has not been peer reviewed.

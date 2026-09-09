# Verification certificates for C42 and C51

Author: [André Röhrig](https://github.com/AndreRoehrig). September 9, 2026.

This repository contains the proofs, exact inputs, and verification programs
for two proposed improvements to the
[Optimization Constants in Mathematics](https://teorth.github.io/optimizationproblems/)
collection.

| Constant | Previous bound | New bound |
| --- | --- | --- |
| [C42: Turán's pure power-sum constant](https://teorth.github.io/optimizationproblems/constants/42a.html) | Upper: 0.6906538 | **Upper: 0.688983** |
| [C51: Erdős maximum-term constant](https://teorth.github.io/optimizationproblems/constants/51a.html) | Lower: 0.585078819653 in the cited paper | **Lower: 0.5850788196744** |

C42 uses an eight-step generalization of Griego's coefficient-cancellation
construction. Exact rational inputs and interval arithmetic verify the
profile radii and a strict limiting inequality. The proof handles integer
rounding and the singular endpoint in the limiting argument. No explicit
finite threshold N is supplied.

C51 uses a new parameter pair in the He–Tang family. Bernstein subdivision
bounds a degree-19 polynomial over the entire interval [0,1], and an explicit
tail bound covers the infinite series. This is a small improvement: the
claimed value exceeds a verified upper bound on the actual previous
specified function by more than 2.1047e-11. The comparison uses the full
prior result, which is more precise than the site's displayed 0.5850788.

## Run the verification

From the repository root, using Python 3:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/certify_turan.py
.venv/bin/python scripts/certify_maximum_term.py
.venv/bin/python scripts/audit.py
```

The only third-party dependency is the pinned mpmath package. Run without
Python's `-O` flag: assertions are part of the verification, and the scripts
explicitly reject that mode.

Repeat the interval computations at higher precision:

```sh
.venv/bin/python scripts/certify_turan.py --dps 200 --terms 480
.venv/bin/python scripts/certify_maximum_term.py --dps 200
```

Both interval certificates pass at 100 and 200 decimal places. Each script
exits with status zero only if its checks pass. The audit additionally tests
the mathematical formulas, finite coefficient cancellation, exact Bernstein
subdivision, and rejection of invalid claims. Its direct quadrature is a
separate numerical check; the interval scripts provide the rigorous bounds.

## Proofs and data

- [C42 proof](notes/turan.md), [exact parameters](certificates/turan.json),
  [expected verification result](certificates/turan_verification.json).
- [C51 proof](notes/maximum_term.md), [exact parameters](certificates/maximum_term.json),
  [expected verification result](certificates/maximum_term_verification.json).
- [Additional mathematical checks](scripts/audit.py).
- [SHA-256 manifest](MANIFEST.sha256), checked with `sha256sum -c MANIFEST.sha256`.

The proofs credit the preceding constructions of Griego, Biró, He–Tang,
and Sothanaphan, and link to their sources. The certificates use exact
rational data; floating-point search results are not used as proof.

## Review and AI assistance

The mathematical exploration, verification programs, and exposition were
developed with AI assistance using Codex. André Röhrig reviewed the proofs,
references, and supporting verification material before publication. The
verification commands were executed and checked during preparation.
These are proposed computer-assisted results; they have not been peer reviewed.

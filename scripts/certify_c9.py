"""Replay the C7 code construction with integer counts from the 367 base words.

Run python3 scripts/certify_c9.py. The optional argument is an offline copy of Tang's
pinned source tarball. The SHA-256 is checked before loading any of its code.
"""

import argparse
import hashlib
import io
import json
import sys
import tarfile
import tempfile
from fractions import Fraction
from pathlib import Path
from urllib.request import urlopen

from c9_patterns import extend_patterns


COMMIT = "b5d6c2c2c981f89f5c62f3db9d26dab045464f81"
ARCHIVE = f"https://codeload.github.com/gomyway/c9-shannon-c7-lower-bound-certificate/tar.gz/{COMMIT}"
SHA256 = "ed617f4062c1c2b5fe39b306875b1655f848fc0cb75452f2d73550cadb3e90cf"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def extract_source(raw, target):
    require(hashlib.sha256(raw).hexdigest() == SHA256, "Unexpected source archive")
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as archive:
        for entry in archive:
            if entry.isfile():
                parts = Path(entry.name).parts
                require(parts[0] == f"c9-shannon-c7-lower-bound-certificate-{COMMIT}"
                        and ".." not in parts, "Unexpected archive path")
                path = target.joinpath(*parts[1:])
                path.parent.mkdir(exist_ok=True, parents=True)
                path.write_bytes(archive.extractfile(entry).read())


def evaluate(t, data, patterns, base):
    leaves = {"G5": (base, 5), "sG5": (t.sigma(base), 5)}
    for st in data["gadgets"]:
        require(st["name"] not in leaves, "Repeated gadget name")
        if st["op"] == "sigma":
            r, dim = leaves[st["of"]]
            value = t.sigma(r), dim
        elif st["op"] == "gao":
            left, dl = leaves[st["L"]]
            right, dr = leaves[st["R"]]
            value = t.gao(left, right), dl+dr
        else:
            require(st["op"] == "hetgao", "Unknown construction")
            left, kw, dl = patterns[st["pattern"]]
            right, dr = leaves[st["R"]]
            value = t.hetgao(left, right, **kw), dl+dr
        leaves[st["name"]] = value
    vectors = {"Phi(" + k + ")": (t.Phi(p), d) for k, (p, d) in leaves.items()}
    top = data["top"]
    rule = t.bpz.RULES[top["rule"]]
    code = t.bpz.CODES[top["terminal_code"]]
    require(not t.bpz.rule_problems(rule) and not t.bpz.code_problems(code),
            "Invalid BPZ combining rule or terminal code")
    arity = len(code[0])
    branches = top["subtrees"] if "subtrees" in top else [top["rule_inputs"]] * arity
    require(len(branches) == arity, "Terminal code arity disagrees")
    subtrees = []
    dim = 0
    for branch in branches:
        values = [vectors[n] for n in branch]
        dim += sum(d for _, d in values)
        subtrees.append(t.bpz.apply_rule(rule, [p for p, _ in values]))
    size = t.bpz.apply_code(code, subtrees)
    require(size == int(data["M"]) and dim == data["dimension"], "Size or dimension disagrees")
    return size, dim


def check(t, data, upstream):
    ledger = t.C.load_certificates(recomputed=True)
    theorem = t.C.theorem2(ledger)
    patterns = t.tan_patterns(ledger, theorem)
    base, _ = t.C.base_profiles(ledger)
    for step in upstream["gadgets"]:
        name = step.get("pattern")
        if name is not None and name not in patterns:
            t.chain_pattern(ledger, theorem, patterns, name)
    old, old_d = evaluate(t, upstream, patterns, base)
    requested = set(patterns) | {s["pattern"] for s in data["gadgets"] if "pattern" in s}
    extend_patterns(t, theorem, base, patterns, requested)
    size, dim = evaluate(t, data, patterns, base)
    require(size**old_d > old**dim, "No strict improvement over Tang")
    q = Fraction(data["lower_bound"])
    require(size * q.denominator**dim > q.numerator**dim and
            q.numerator**old_d > old*q.denominator**old_d, "Decimal bound fails")
    print(f"Verified dimension {dim}, cardinality {size}")
    print(f"Candidate root = {t.C.nth_root_digits(size, dim, 24)}...")
    print(f"Tang's root = {t.C.nth_root_digits(old, old_d, 24)}...")
    print(f"Theta(C7) > {data['lower_bound']}")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Use Python without -O so upstream assertions are checked.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", nargs="?", type=Path)
    args = parser.parse_args()
    if args.archive:
        raw = args.archive.read_bytes()
    else:
        with urlopen(ARCHIVE, timeout=60) as response:
            raw = response.read()
    data = json.loads((Path(__file__).resolve().parents[1] / "certificates" / "c9.json").read_bytes())
    with tempfile.TemporaryDirectory() as directory:
        target = Path(directory)
        extract_source(raw, target)
        sys.path.insert(0, str(target))
        import check_c9 as t
        check(t, data, json.loads((target/"certificate.json").read_bytes()))

"""Count side codebooks for new left-associated products of Gao gadgets.

The two bits of a flag tell whether a word is in the closed neighborhood of
the target gadget's neutral set and its union of one-sided sets. For a Gao
product, two equal input classes give a neutral output and two different
classes give a one-sided output. Thus both neighborhood flags propagate by
the same parity rule, including words belonging to both neighborhoods.
"""

from functools import lru_cache
from itertools import product


# For each class of the source Gao product, list its disjoint product blocks.
# Rows and indices are B, R, P^H, P^V, X^0, X^H, X^V.
BLOCKS = (
    ((0, 0), (2, 5), (3, 6), (5, 3), (6, 2)),
    ((1, 4), (4, 1)),
    ((2, 4), (4, 3)),
    ((3, 4), (4, 2)),
    ((4, 4), (5, 5), (5, 6), (6, 5), (6, 6)),
    ((5, 4), (4, 6)),
    ((6, 4), (4, 5)),
)

BASE_LEFT = {"P10": "G10", "P10A": "G10A", "P10D": "G10D", "P15X": "G15X",
             "P15A8": "G15AX", "P15Ah": "G15AX", "P15A5": "G15AX",
             "P15D": "G15DX", "P15w": "G15h"}
TYPES = {"G5": "5", "G10": "10", "G10A": "10", "G10D": "10", "G15X": "15X",
         "G15AX": "15X", "G15DX": "15X", "G15h": "15h", "G15Ah": "15A", "G15Dh": "15A"}


def factors(name):
    if name.startswith("S<"):
        _, right, pattern = name[2:-1].split(",")
        return pattern, right
    if name.startswith("N<"):
        inner, right = name[2:-1].rsplit(">,", 1)
        return inner + ">", right
    raise ValueError(f"Unknown product {name}")


def named_profile(t, theorem, base, name):
    original = name.removeprefix("s")
    kind = TYPES[original]
    dim = 5 if kind == "5" else 10 if kind == "10" else 15
    profile = base if dim == 5 else theorem[original]
    return (t.sigma(profile) if original != name else profile), dim


def combine(first, second):
    rows = []
    for blocks in BLOCKS:
        counts = [0, 0, 0, 0]
        for i, j in blocks:
            for a, b in product(range(4), repeat=2):
                neutral = bool(a & b)
                one_sided = bool((a & 1 and b & 2) or (a & 2 and b & 1))
                counts[int(neutral) + 2*int(one_sided)] += first[i][a] * second[j][b]
        rows.append(counts)
    return rows


def base_tables():
    from shannon import certificates as C

    g, transform, exchanges, _ = C.load_base()
    fifteen, _, _ = C.build(g, transform, exchanges)
    atoms = {a.name: a for blocks in fifteen.values() for row in blocks for a in row}
    def one(*atoms):
        return [(a,) for a in atoms]
    five = [one(getattr(g, k)) for k in ("B", "R", "PH", "PV", "X0", "XH", "XV")]
    ten = [one(atoms[k]) for k in ("B2", "R2", "PH2", "PV2", "X2_0", "X2_H", "X2_V")]
    main15 = [fifteen[k] for k in ("B15", "R15", "PH15", "PV15", "X15_0", "X15_H", "X15_V")]
    data = {"5": (five, one(g.X0), one(g.XH, g.XV)),
            "10": (ten, one(atoms["X2_0"]), one(atoms["X2_r"])),
            "15X": (main15, fifteen["X15X_0"], fifteen["X15X_r"]),
            "15h": (main15, fifteen["X15_0"], fifteen["X15_H"] + fifteen["X15_V"]),
            "15A": (main15, fifteen["X15A_0"], fifteen["X15A_r"])}
    result = {}
    for name, (classes, neutral, side) in data.items():
        rows = []
        for source in classes:
            words = C.apply_T_blocks(source, transform)
            neither = C.count_outside(words, neutral + side)
            outside_neutral = C.count_outside(words, neutral)
            outside_side = C.count_outside(words, side)
            n = C.size(words)
            row = [neither, outside_side-neither, outside_neutral-neither,
                   n-outside_neutral-outside_side+neither]
            if min(row) < 0 or sum(row) != n:
                raise ValueError("Invalid neighborhood partition")
            rows.append(row)
        result[name] = rows
    return result


def extend_patterns(t, theorem, base, patterns, requested):
    tables = base_tables()
    @lru_cache(maxsize=None)
    def walk(name):
        inner, right_name = factors(name)
        if name.startswith("S<"):
            initial = name[2:-1].split(",")[0]
            if BASE_LEFT.get(inner) != initial:
                raise ValueError("Mismatched first factor")
            table = tables[TYPES[initial]]
            source_name = {"5": "G5", "10": "G10", "15": "G15h"}[str(patterns[inner][2])]
            reference, _ = named_profile(t, theorem, base, source_name)
        else:
            table, reference = walk(inner)
        target, kw, dl = patterns[inner]
        right, dr = named_profile(t, theorem, base, right_name)
        kind = TYPES[right_name.lstrip("s")]
        table = combine(table, tables[kind])
        source_name = {5: "G5", 10: "G10", 15: "G15h"}[dr]
        source_right, _ = named_profile(t, theorem, base, source_name)
        reference = t.gao(reference, source_right)
        # The source main code is B union R. Flags 0 and 2 lie outside N(X^0).
        j = sum(sum(row) for row in table[:2])
        q = sum(row[0] + row[2] for row in table[:2])
        if j != reference.a:
            raise ValueError("Count does not match the source main code")
        left = t.gao(target, right)
        sibling = t.hetgao(target, right, **kw)
        if left.as_tuple()[:2] != sibling.as_tuple()[:2]:
            raise ValueError("The two ordered products do not share code and pairs")
        side = t.SideCodebook(j, q)
        opts = dict(J0=t.neutral_from(sibling), JH=side, JV=side)
        if name in patterns:
            prev, params, dim = patterns[name]
            if (left.as_tuple() != prev.as_tuple() or dim != dl+dr or
                    (j, q) != (params["JH"].j, params["JH"].q)):
                raise ValueError("Disagreement with Tang's direct block counter")
        patterns[name] = left, opts, dl+dr
        print(f"Verified {name}, dimension {dl+dr}, side size {j}, outside count {q}")
        return table, reference

    for name in sorted(requested):
        if name.startswith(("S<", "N<")):
            walk(name)

"""Verify a compressed 54-point family and its finite-chain capacity bound.

Run python3 scripts/certify_c49.py. All masks, sets, and counts are explicit in
certificates/c49.json. Each half uses seven core bits and two disjoint ten-point
pendants. A box is [core mask, first pendant family, second pendant family].
One top-level product pairs a union of left boxes with a union of right boxes.

To glue two blocks, keep pairs with opposite bits at the joined ports and
erase those bits. Suppose two projected members are comparable. If their
deleted bits differ, restoring them makes the containment strict in one of
the two original blocks, contradicting its antichain property. If both
deleted bits agree, containment forces equality in both blocks. This also
gives uniqueness of the representation used to count the chain.

For three intersection-equal projected members, the restored port on at
least one side has zero, one or three ones, so the triple on that side is
still intersection-equal. In a sunflower-free antichain every such triple
is constant. Its three port bits therefore agree, and so do their opposites
on the other side. That other triple is constant as well. Thus gluing
preserves both properties. The same coordinate argument applies to tensor
powers of the finite chain and proves the claimed capacity lower bound.
"""

import json
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import cache
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def subset(x, y):
    return x & y == x


def flower(a, b, c):
    return a & b == a & c == b & c


def verify(data):
    pendants = data["pendant_families"]
    for family in pendants:
        require(family and all(0 <= x < 1024 for x in family), "Invalid pendant mask")
        require(all(not (subset(a, b) or subset(b, a)) for a, b in combinations(family, 2)),
                "A pendant family contains comparable members")
        require(all(not flower(*triple) for triple in combinations(family, 3)),
                "A pendant family contains a sunflower")

    ids = range(len(pendants))
    le = {(i, j): any(subset(a, b) for a, b in product(pendants[i], pendants[j]))
          for i, j in product(ids, repeat=2)}
    @cache
    def possible_thirds(a, b):
        # A third mask must contain a & b, omit a ^ b, and may use any
        # submask of the coordinates outside their union.
        free = 1023 ^ (a | b)
        sub = free
        bits = 0
        while True:
            bits |= 1 << ((a & b) | sub)
            if sub == 0:
                return bits
            sub = (sub-1) & free

    memberships = [sum(1 << a for a in fam) for fam in pendants]
    can_flower = {}
    for i, j in combinations_with_replacement(ids, 2):
        thirds = 0
        for a, b in product(pendants[i], pendants[j]):
            thirds |= possible_thirds(min(a, b), max(a, b))
        for k, members in enumerate(memberships):
            can_flower[i, j, k] = can_flower[j, i, k] = bool(thirds & members)

    def boxes_le(u, v):
        return subset(u[0], v[0]) and le[u[1], v[1]] and le[u[2], v[2]]

    def boxes_flower(u, v, w):
        return (flower(u[0], v[0], w[0]) and can_flower[u[1], v[1], w[1]]
                and can_flower[u[2], v[2], w[2]])

    @cache
    def check_half(boxes, forbidden):
        require(boxes, "Empty half")
        for core, i, j in boxes:
            require(0 <= core < 1024 and core & forbidden == 0 and i in ids and j in ids,
                    "Invalid core or pendant index")
        for u, v in combinations(boxes, 2):
            require(not (boxes_le(u, v) or boxes_le(v, u)), "Overlapping or comparable boxes")
        for i, j, k in combinations_with_replacement(range(len(boxes)), 3):
            if i != k:
                require(not boxes_flower(boxes[i], boxes[j], boxes[k]), "Sunflower inside half")

    # Inside a single box the product of Sperner sunflower-free factors is
    # again Sperner and sunflower-free. An intersection-equal triple in any
    # Sperner family is either constant or consists of three distinct sets.
    products = [(tuple(map(tuple, x)), tuple(map(tuple, y))) for x, y in data["products"]]
    left, right = zip(*products)
    for halves, forbidden in ((left, 2 | 4 | 8), (right, 1 | 4 | 8)):
        for half in halves:
            check_half(half, forbidden)

    def comparable(first, second):
        return any(boxes_le(u, v) for u in first for v in second)

    for i, j in combinations(range(len(products)), 2):
        for a, b in ((i, j), (j, i)):
            require(not (comparable(left[a], left[b]) and comparable(right[a], right[b])),
                    f"Equal or comparable members in products {a}, {b}")

    def completion_table(halves):
        # A bit represents a literal box occurrence. All three coordinate
        # conditions must hold for that same box. Collapse to product indices
        # only after taking their intersection.
        n = len(halves)
        slots = []
        ranges = {}
        next_slot = n
        for i, half in enumerate(halves):
            if len(half) == 1:
                slots.append((i, half[0]))
            else:
                ranges[i] = ((1 << len(half)) - 1) << next_slot
                slots.extend((next_slot + j, box) for j, box in enumerate(half))
                next_slot += len(half)

        cores = {box[0] for _, box in slots}
        def select(predicate):
            return sum(1 << i for i, box in slots if predicate(box))
        core_bits = {(a, b): select(lambda box: flower(a, b, box[0]))
                     for a, b in product(cores, repeat=2)}
        pendant_bits = {(axis, i, j): select(lambda box: can_flower[i, j, box[axis]])
                        for axis in (1, 2) for i, j in product(ids, repeat=2)}
        def completions(i, j):
            bits = 0
            for a, b, c in halves[i]:
                for d, e, f in halves[j]:
                    bits |= core_bits[a, d] & pendant_bits[1, b, e] & pendant_bits[2, c, f]
            labels = bits & ((1 << n) - 1)
            for label, mask in ranges.items():
                if bits & mask:
                    labels |= 1 << label
            return labels
        return completions

    left_completions = completion_table(left)
    right_completions = completion_table(right)
    for i, j in combinations_with_replacement(range(len(products)), 2):
        hits = left_completions(i, j) & right_completions(i, j)
        # Triples with all three labels equal have already been checked.
        require(hits >> (j + int(i == j)) == 0, f"Sunflower involving products {i} and {j}")

    def port_counts(halves, mask):
        return [sum(len(pendants[a]) * len(pendants[b]) for c, a, b in halves
                    if bool(c & mask) == bit) for bit in (0, 1)]
    def counts(family):
        N = [[0, 0], [0, 0]]
        for lhalf, rhalf in family:
            l, r = port_counts(lhalf, 1), port_counts(rhalf, 2)
            for p, q in product(range(2), repeat=2):
                N[p][q] += l[p] * r[q]
        return N
    N = counts(products)
    old = [[4397570400, 7211840834], [7211840834, 11504366352]]
    require(counts(products[:4]) == old and N == data["counts"], "Port counts disagree")
    require(all(N[i][j] >= old[i][j] for i, j in product(range(2), repeat=2)) and N != old,
            "No strict Perron improvement")
    q = Fraction(data["finite_lower_bound"])
    b, a = N[0]
    d, c = N[1]
    lhs = 2 * q**52 - a - d
    require(lhs < 0 or lhs**2 < (a-d)**2 + 4*b*c, "Claim exceeds the new rate")
    old_a, old_b, old_c = old[0][1], old[0][0], old[1][1]
    require(q**52 > old_a and (q**52 - old_a)**2 > old_b*old_c,
            "Claim does not exceed Canli's rate")

    # DEGKM's opposite-port glue gives v(k+1) = M v(k) and n(k) = 52k+2.
    # Directed integer rounding certifies one finite chain at the stated q.
    M = [[a, b], [c, d]]
    scale = 10**90
    H = [[int(scale * x / q**52) for x in row] for row in M]
    for _ in range(data["squarings"]):
        H = [[sum(H[i][k] * H[k][j] for k in range(2)) // scale
              for j in range(2)] for i in range(2)]
    initial = [sum(row) for row in N]
    size = sum(H[i][j] * initial[j] for i, j in product(range(2), repeat=2))
    require(size > scale * q**54, "Finite chain too short to prove the claim")
    with localcontext() as ctx:
        ctx.prec = 40
        rho = (a + d + Decimal((a-d)**2 + 4*b*c).sqrt())/2
        print(f"Rate = {rho**(Decimal(1)/52)}")
    print(f"Verified {len(products)} disjoint products on 54 points, {sum(map(sum, N))} sets")
    print(f"Port counts = {N}. Added sets = {sum(map(sum,N)) - sum(map(sum,old))}")
    print(f"Finite chain of {2**data['squarings']+1} blocks proves C49 > {data['finite_lower_bound']}")


if __name__ == "__main__":
    verify(json.loads((Path(__file__).resolve().parents[1] / "certificates" / "c49.json").read_bytes()))

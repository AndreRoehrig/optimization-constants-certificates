"""Check the optimized C47 construction exactly with Python's standard library.

Run python3 scripts/certify_c47.py. The gzipped certificate is self-contained.
Decompression also produces Pinho's c47-witness-rectangles-v1 input format.
"""

import gzip
import json
import math
from bisect import bisect_left, bisect_right
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def largest_rectangle(box, radius, px, py):
    left, right, bottom, top = box
    return (max(0, right - radius), min(px, left + radius),
            max(0, top - radius), min(py, bottom + radius))


def scaled_integer(value, unit):
    return value.numerator * (unit // value.denominator)


def check_masses(points, px, py, boxes):
    """Count every periodic copy with coordinates inside a closed box."""
    coord_unit = math.lcm(px.denominator, py.denominator,
                          *(p[i].denominator for p in points for i in (0, 1)),
                          *(v.denominator for box in boxes for v in box))
    mass_unit = math.lcm(*(m.denominator for _, _, m in points))
    period_x, period_y = (scaled_integer(v, coord_unit) for v in (px, py))
    int_boxes = [tuple(scaled_integer(v, coord_unit) for v in b) for b in boxes]
    x_min = min(b[0] for b in int_boxes)
    x_max = max(b[1] for b in int_boxes)
    y_min = min(b[2] for b in int_boxes)
    y_max = max(b[3] for b in int_boxes)
    copies = []
    for x, y, mass in points:
        x, y = scaled_integer(x, coord_unit), scaled_integer(y, coord_unit)
        weight = scaled_integer(mass, mass_unit)
        for k in range(-((x - x_min) // period_x), (x_max - x) // period_x + 1):
            for j in range(-((y - y_min) // period_y), (y_max - y) // period_y + 1):
                copies.append((x + k * period_x, y + j * period_y, weight))
    copies.sort()
    xs = [x for x, _, _ in copies]
    for left, right, bottom, top in int_boxes:
        lo, hi = bisect_left(xs, left), bisect_right(xs, right)
        total = sum(m for _, y, m in copies[lo:hi] if bottom <= y <= top)
        yield F(total, mass_unit)


def union_area(rectangles):
    """Sweep vertical edges while tracking covered length on the y-axis."""
    unit = math.lcm(*(v.denominator for rect in rectangles for v in rect))
    rects = [tuple(scaled_integer(v, unit) for v in rect) for rect in rectangles]
    ys = sorted({y for _, _, bottom, top in rects for y in (bottom, top)})
    index = {y: i for i, y in enumerate(ys)}
    intervals = len(ys) - 1
    cover = [0] * (4 * intervals)
    length = [0] * (4 * intervals)

    def update(node, lo, hi, start, stop, delta):
        if start <= lo and hi <= stop:
            cover[node] += delta
        else:
            mid = (lo + hi) // 2
            if start < mid:
                update(2 * node, lo, mid, start, stop, delta)
            if stop > mid:
                update(2 * node + 1, mid, hi, start, stop, delta)
        if cover[node] > 0:
            length[node] = ys[hi] - ys[lo]
        elif hi == lo + 1:
            length[node] = 0
        else:
            length[node] = length[2 * node] + length[2 * node + 1]

    events = []
    for left, right, bottom, top in rects:
        events.append((left, index[bottom], index[top], 1))
        events.append((right, index[bottom], index[top], -1))
    events.sort()
    area = 0
    previous = events[0][0]
    for x, start, stop, delta in events:
        area += (x - previous) * length[1]
        update(1, 0, intervals, start, stop, delta)
        previous = x
    require(length[1] == 0, "Unbalanced rectangle edges")
    return F(area, unit**2)


def check(data):
    m = data["measure"]
    px, py, level = (F(m[k]) for k in ("Px", "Py", "lam"))
    points = [tuple(map(F, p)) for p in m["pts"]]
    require(points and min(px, py, level) > 0 and all(0 <= x < px and 0 <= y < py and v > 0 for x, y, v in points),
            "Invalid periodic measure")
    require(len({(x, y) for x, y, _ in points}) == len(points), "Repeated atom")
    witnesses = data["witnesses"]
    require(witnesses, "Missing witness rectangles")
    boxes = [tuple(map(F, w["box"])) for w in witnesses]
    rects, radii, slacks = [], [], []
    for w, box, mass in zip(witnesses, boxes, check_masses(points, px, py, boxes)):
        r = F(w["radius"])
        rect = tuple(map(F, w["rectangle"]))
        slack = mass - 4 * level * r**2
        require(box[0] <= box[1] and box[2] <= box[3] and r > 0 and slack > 0,
                "Witness square fails the strict mass test")
        require(largest_rectangle(box, r, px, py) == rect and rect[0] < rect[1]
                and rect[2] < rect[3], "Rectangle of centers disagrees")
        rects.append(rect)
        radii.append(r)
        slacks.append(slack)
    ratio = level * union_area(rects) / sum(v for _, _, v in points)
    target = F(data["finite_lower_bound"])
    # Periodic ratio of Pinho’s round-nine construction in PR #350.
    prior = F("22878843460258009274861600625288681562329176839147945011303/"
              "13414024367090776613897918866165655414442557440000000000000")
    require(ratio == F(data["periodic_ratio"]) and ratio > target > prior,
            "Exact ratio, claim, or comparison with Pinho failed")

    # An h-cell margin keeps every witness square inside the support block.
    # Smearing each atom on a square of radius eps enlarges every required
    # witness square by eps, preserving the strict mass comparison.
    r_max = max(radii)
    h = 1 + r_max // min(px, py)
    n = max(2*h+1, 1 + math.ceil(4*h*ratio/(ratio-target)))
    eps = min(F(1, 2), min(slacks)/(16*level*(2*r_max+1)))
    finite = ratio * F((n-2*h)**2, n**2)
    require(h*min(px, py) > r_max and eps > 0 and
            4*level*(2*r_max+1)*eps < min(slacks) and finite > target,
            "Finite truncation or atom smearing failed")
    with localcontext() as ctx:
        ctx.prec = 35
        print(f"Periodic ratio = {Decimal(ratio.numerator)/Decimal(ratio.denominator)}")
    print(f"Verified {len(points)} atoms, {len(witnesses)} rectangles")
    print(f"Exact periodic ratio = {ratio}")
    print(f"N = {n}, h = {h}, epsilon = {eps}")
    print(f"Finite construction proves C47 > {data['finite_lower_bound']}")


if __name__ == "__main__":
    with gzip.open(Path(__file__).resolve().parents[1] / "certificates" / "c47.json.gz", "rt") as f:
        check(json.load(f))

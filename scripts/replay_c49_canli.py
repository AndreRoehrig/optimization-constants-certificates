"""Check the 54-point block and its four original products with Canli's referee."""

import argparse
import hashlib
import json
import sys
from pathlib import Path


SOURCE_HASHES = {
    "c49_base.py": "2ffc555f7126e6cb618cf3f0b139886b0818f1cf06c771fb46fbae84b0b6969a",
    "c49_referee.py": "44ff70f566efcb7f35f7b09c2cd0aef31871af747dec33db71f06a830136c894",
}


def literal_products(products, pendants):
    """Compare the two encodings after expanding each pendant label to its masks."""
    return [
        [sorted((core, tuple(sorted(pendants[i])), tuple(sorted(pendants[j])))
                for core, i, j in half) for half in product]
        for product in products
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("referee_dir", type=Path, help="Canli's improvements/c49/referee directory")
    directory = parser.parse_args().referee_dir.resolve()
    if not __debug__:
        raise SystemExit("Use Python without -O so upstream assertions are checked.")
    for name, expected in SOURCE_HASHES.items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Unexpected source file {name}")
    sys.path.insert(0, str(directory))
    import c49_referee as referee

    path = Path(__file__).resolve().parents[1] / "certificates" / "c49.json"
    data = json.loads(path.read_bytes())
    products = data["products"]
    pendants = data["pendant_families"]
    registry = referee.Fams()
    original, valid, _, _ = referee.build54(registry, xstar_left=12, xstar_right=12)
    if not valid or literal_products(products[:4], pendants) != literal_products(original, registry.fams):
        raise ValueError("The first four products do not reproduce Canli's construction")
    tables = referee.PendTables([tuple(masks) for masks in pendants])
    verdict = referee.decide(products, tables, want_witness=False)
    if any(verdict.values()):
        raise ValueError(f"The referee found a counterexample: {verdict}")
    counts = referee.counts(products, tables, referee.PORT_P, referee.PORT_Q)
    if counts != data["counts"]:
        raise ValueError("Port counts disagree")
    print(f"Canli's original block and all {len(products)} products verified.")
    print(f"Port counts = {counts}")


if __name__ == "__main__":
    main()

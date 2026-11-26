"""Generate feasible subtour points and bound their worst metric ratios."""

import argparse
import ctypes
import heapq
import json
import math
import random
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parent
SCALE = 10**8
INT64_ARRAY = np.ctypeslib.ndpointer(dtype=np.int64, flags="C_CONTIGUOUS")


class Exact:
    def __init__(self):
        filename = "exact.dll" if hasattr(ctypes, "windll") else "exact.so"
        self.library = ctypes.CDLL(str(ROOT / filename))
        for function in ("minimum_cut", "shortest_tour"):
            routine = getattr(self.library, function)
            routine.argtypes = (ctypes.c_int, INT64_ARRAY, ctypes.POINTER(ctypes.c_int))
            routine.restype = ctypes.c_int64

    def minimum_cut(self, n, edges, weights):
        matrix = np.zeros((n, n), dtype=np.int64)
        for (u, v), weight in zip(edges, weights, strict=True):
            matrix[u, v] = matrix[v, u] = weight
        mask = ctypes.c_int()
        value = self.library.minimum_cut(n, matrix, ctypes.byref(mask))
        return value, mask.value

    def shortest_tour(self, distance):
        n = len(distance)
        path = (ctypes.c_int * n)()
        length = self.library.shortest_tour(n, distance, path)
        return length, list(path)


def triangle_paths(lengths):
    edges = {(u, v): Fraction(1, 2) for offset in (0, 3)
             for u in range(offset, offset + 3) for v in range(u + 1, offset + 3)}
    n = 6
    for endpoint, length in enumerate(lengths):
        path = [endpoint, *range(n, n + length - 1), endpoint + 3]
        n += length - 1
        for u, v in zip(path, path[1:]):
            edges[tuple(sorted((u, v)))] = Fraction(1)
    return n, dict(sorted(edges.items()))


def subtour_point(exact, n, edges, objective):
    degrees = np.zeros((n, len(edges)))
    for index, (u, v) in enumerate(edges):
        degrees[u, index] = degrees[v, index] = 1
    cuts = {}
    while True:
        result = linprog(
            objective, A_ub=list(cuts.values()) or None, b_ub=[-2]*len(cuts) or None,
            A_eq=degrees, b_eq=np.full(n, 2), bounds=(0, 1), method="highs-ds",
        )
        if not result.success:
            raise RuntimeError(result.message)
        capacities = np.rint(result.x * SCALE).astype(np.int64)
        minimum, mask = exact.minimum_cut(n, edges, capacities)
        if minimum >= 2 * SCALE - len(edges):
            point = {e: Fraction(float(x)).limit_denominator(10_000)
                     for e, x in zip(edges, result.x, strict=True) if x > 1e-8}
            denominator = math.lcm(*(x.denominator for x in point.values()))
            if denominator > SCALE:
                raise RuntimeError("Unexpected denominator in subtour point")
            numerators = [int(x * denominator) for x in point.values()]
            if exact.minimum_cut(n, point, numerators)[0] != 2 * denominator:
                raise RuntimeError("Rational reconstruction violated a subtour cut")
            for v in range(n):
                if sum(x for e, x in point.items() if v in e) != 2:
                    raise RuntimeError("Rational reconstruction violated a degree equation")
            tight = [mask for mask in cuts
                     if sum(x for (u, v), x in point.items()
                            if ((mask >> u) ^ (mask >> v)) & 1) == 2]
            return point, tight
        if mask in cuts:
            raise RuntimeError("Previously separated cut still violated")
        cuts[mask] = [-(((mask >> u) ^ (mask >> v)) & 1) for u, v in edges]


def shortest_paths(n, edges, weights):
    graph = [[] for _ in range(n)]
    for index, ((u, v), weight) in enumerate(zip(edges, weights, strict=True)):
        graph[u].append((v, int(weight), index))
        graph[v].append((u, int(weight), index))
    infinity = 1 << 60
    distance = np.full((n, n), infinity, dtype=np.int64)
    paths = [[None] * n for _ in range(n)]
    for s in range(n):
        distance[s, s] = 0
        queue = [(0, s)]
        while queue:
            cost, u = heapq.heappop(queue)
            if cost != distance[s, u]:
                continue
            for v, weight, index in graph[u]:
                candidate = cost + weight
                if candidate < distance[s, v]:
                    distance[s, v] = candidate
                    paths[s][v] = (u, index)
                    heapq.heappush(queue, (candidate, v))
    return distance, paths


def graphic_tour(exact, n, edges, weights):
    distance, paths = shortest_paths(n, edges, weights)
    length, tour = exact.shortest_tour(distance)
    usage = [0] * len(edges)
    for s, v in zip(tour, tour[1:] + tour[:1]):
        while v != s:
            v, index = paths[s][v]
            usage[index] += 1
    if sum(int(w) * a for w, a in zip(weights, usage, strict=True)) != length:
        raise RuntimeError("Expanded tour length did not match the solver")
    return length, tuple(usage)


def exact_lower_bound(exact, n, edges, point, weights):
    divisor = math.gcd(*weights)
    weights = [w // divisor for w in weights]
    length, _ = graphic_tour(exact, n, edges, weights)
    cost = sum(x * w for x, w in zip(point, weights, strict=True))
    return Fraction(length) / cost, weights


def bound_gap(exact, n, point, max_rounds):
    edges = list(point)
    fractions = list(point.values())
    objective = list(map(float, fractions))
    _, initial = graphic_tour(exact, n, edges, [1] * len(edges))
    rows = [initial]
    known = {initial}
    best_lower = Fraction(0)
    best_upper = Fraction(2)
    best_weights = None
    best_mixture = None
    for rounds in range(1, max_rounds + 1):
        result = linprog(objective, A_ub=-np.asarray(rows), b_ub=-np.ones(len(rows)),
                         bounds=(0, None), method="highs-ds")
        if not result.success:
            raise RuntimeError(result.message)
        raw_mass = [max(0, -float(y)) for y in result.ineqlin.marginals]
        exact_mass = [Fraction(y).limit_denominator(100_000) for y in raw_mass]
        total_mass = sum(exact_mass)
        normalized = [mass / total_mass for mass in exact_mass]
        upper = max(sum(mass * row[i] for mass, row in zip(normalized, rows, strict=True))
                    / x for i, x in enumerate(fractions))
        if upper <= best_upper or best_mixture is None:
            best_upper = upper
            best_mixture = [[str(p), list(row)] for p, row in zip(normalized, rows, strict=True)
                            if p > 0]

        weights = [max(0, int(round(w * SCALE))) for w in result.x]
        length, usage = graphic_tour(exact, n, edges, weights)
        approximate_cost = sum(x * w for x, w in zip(fractions, weights, strict=True))
        lower = Fraction(length) / approximate_cost
        if lower > best_lower:
            best_lower = lower
            best_weights = weights
        if length >= SCALE - 5 * n or usage in known:
            # Recover simple integral optimal weights where possible.
            w = [Fraction(float(c)).limit_denominator(100_000) for c in result.x]
            denominator = math.lcm(*(c.denominator for c in w))
            if denominator < SCALE:
                lower, candidate = exact_lower_bound(exact, n, edges, fractions,
                                                     [int(c * denominator) for c in w])
                if lower >= best_lower:
                    best_lower, best_weights = lower, candidate
            break
        known.add(usage)
        rows.append(usage)
    if best_weights is None or best_lower > best_upper:
        raise RuntimeError("Inconsistent cost and mixture bounds")
    return {
        "edges": [[*e, str(x)] for e, x in point.items()],
        "lower": str(best_lower), "upper": str(best_upper),
        "lengths": best_weights, "mixture": best_mixture,
        "rounds": rounds,
    }


def run(args):
    exact = Exact()
    n, control = triangle_paths([4, 5, 5])
    rng = random.Random(args.seed)
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n) if (u, v) not in control]
    with args.output.open("x", encoding="utf-8") as output:
        result = bound_gap(exact, n, control, args.rounds)
        output.write(json.dumps({"case": "control", "cities": n, **result}) + "\n")
        output.flush()
        print(f"Control: {result['lower']} to {result['upper']}", flush=True)
        if args.guided:
            reference, _ = shortest_paths(n, control, result["lengths"])
        seen = set()
        for attempt in range(args.attempts):
            allowed = pairs
            if args.local:
                # Leave the unselected interior vertices at degree two.
                # The original path segments incident to them retain weight one.
                active = set(range(6)) | set(rng.sample(range(6, n), rng.randint(5, 8)))
                allowed = [(u, v) for u, v in pairs if u in active and v in active]
            chords = rng.sample(allowed, rng.randint(10, min(19, len(allowed))))
            edges = sorted(set(control) | set(chords))
            objective = [rng.randrange(-1000, 1001) - (attempt % 4) * 350 * (e in control)
                         for e in edges]
            if args.guided:
                noise = 3 * (1 + attempt % 8)
                objective = [1000 * int(reference[e]) + noise * value
                             for e, value in zip(edges, objective, strict=True)]
            point, tight = subtour_point(exact, n, edges, objective)
            if len(point) < n + 7 or all(x.denominator <= 2 for x in point.values()):
                continue
            key = tuple(point.items())
            if key in seen:
                continue
            seen.add(key)
            result = bound_gap(exact, n, point, args.rounds)
            record = {"case": attempt, "cities": n, "chords": chords,
                      "tight_cuts": tight, **result}
            output.write(json.dumps(record) + "\n")
            output.flush()
            print(f"{attempt}: support {len(point)}, "
                  f"{result['lower']} to {result['upper']} ({result['rounds']} rounds)", flush=True)
            if Fraction(result["lower"]) > Fraction(4, 3) or len(seen) == args.cases:
                break
        print(f"Checked {len(seen)} distinct candidate vectors.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--seed", type=int, default=61125)
    parser.add_argument("--attempts", type=int, default=600)
    parser.add_argument("--cases", type=int, default=30)
    parser.add_argument("--rounds", type=int, default=200)
    parser.add_argument("--local", action="store_true", help="restrict chord endpoints")
    parser.add_argument("--guided", action="store_true", help="bias objectives toward the control metric")
    arguments = parser.parse_args()
    if arguments.attempts < 0 or arguments.cases < 1 or arguments.rounds < 1:
        parser.error("expected nonnegative attempts and positive cases and rounds")
    run(arguments)

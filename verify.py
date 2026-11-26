"""Check a result ledger without NumPy, SciPy or the search library."""

import argparse
import json
from array import array
from collections import deque
from fractions import Fraction
from math import lcm
from pathlib import Path


def check(condition, failure):
    if not condition:
        raise ValueError(failure)


def check_subtour(n, edges, point):
    divisor = lcm(*(x.denominator for x in point))
    capacity = [[0] * n for _ in range(n)]
    for (u, v), x in zip(edges, point, strict=True):
        capacity[u][v] = capacity[v][u] = int(divisor * x)
    required = 2 * divisor
    check(all(sum(row) == required for row in capacity), "Degree equation failed")
    # Every cut separates vertex 0 from some terminal. An integral flow
    # of the required value to each terminal proves all cut inequalities.
    for target in range(1, n):
        residual = [row[:] for row in capacity]
        value = 0
        while value < required:
            prev = [-1] * n
            prev[0] = 0
            queue = deque([0])
            while queue and prev[target] < 0:
                u = queue.popleft()
                for v, cap in enumerate(residual[u]):
                    if cap > 0 and prev[v] < 0:
                        prev[v] = u
                        queue.append(v)
            check(prev[target] >= 0, "Cut inequality failed")
            path = []
            v = target
            delta = required - value
            while v > 0:
                u = prev[v]
                path.append((u, v))
                delta = min(delta, residual[u][v])
                v = u
            for u, v in path:
                residual[u][v] -= delta
                residual[v][u] += delta
            value += delta


def check_extreme(n, edges, point, masks):
    # At a point with this support, coordinates outside it are already
    # fixed at zero. Find independent equalities on the remaining edges.
    rows = [[int(v in edge) for edge in edges] for v in range(n)]
    for mask in masks:
        check(type(mask) is int and 0 < mask < (1 << n) and mask % 2 == 0,
              "Invalid cut mask")
        row = [((mask >> u) ^ (mask >> v)) & 1 for u, v in edges]
        check(sum(a * x for a, x in zip(row, point, strict=True)) == 2,
              "Claimed tight cut is slack")
        rows.append(row)
    for i, x in enumerate(point):
        if x == 1:
            rows.append([int(i == j) for j in range(len(edges))])
    pivots = {}
    for row in rows:
        reduced = list(map(Fraction, row))
        for i, coefficient in enumerate(reduced):
            if not coefficient:
                continue
            if i not in pivots:
                pivots[i] = [value / coefficient for value in reduced]
                break
            reduced[i:] = [a - coefficient * b
                           for a, b in zip(reduced[i:], pivots[i][i:], strict=True)]
    check(len(pivots) == len(edges), "Active constraints do not certify a vertex")


def check_mixture(n, edges, point, mixture, upper):
    degrees = [0] * n
    average = [Fraction(0) for _ in edges]
    total = Fraction(0)
    for mass_string, row in mixture:
        mass = Fraction(mass_string)
        check(mass > 0 and len(row) == len(edges), "Invalid mixture term")
        total += mass
        degrees[:] = [0] * n
        graph = [[] for _ in range(n)]
        for i, ((u, v), mult) in enumerate(zip(edges, row, strict=True)):
            check(type(mult) is int and mult >= 0, "Invalid edge multiplicity")
            if mult:
                degrees[u] += mult
                degrees[v] += mult
                graph[u].append(v)
                graph[v].append(u)
                average[i] += mass * mult
        check(all(d > 0 and d % 2 == 0 for d in degrees), "Tour is not Eulerian")
        reached = {0}
        queue = deque([0])
        while queue:
            for v in graph[queue.popleft()]:
                if v not in reached:
                    reached.add(v)
                    queue.append(v)
        check(len(reached) == n, "Tour is disconnected")
    check(total == 1, "Mixture masses do not sum to 1")
    check(max(a / x for a, x in zip(average, point, strict=True)) == upper,
          "Incorrect mixture upper bound")


def minimum_tour(n, edges, lengths):
    infinity = (n * sum(lengths)) + 1
    check(infinity < 2**62, "Tour lengths exceed verifier capacity")
    distances = [[0 if u == v else infinity for v in range(n)] for u in range(n)]
    for (u, v), weight in zip(edges, lengths, strict=True):
        distances[u][v] = distances[v][u] = weight
    for via in range(n):
        for u in range(n):
            for v in range(n):
                distances[u][v] = min(distances[u][v], distances[u][via] + distances[via][v])
    # Forward relaxation on subsets. Every cycle is a path from 0
    # visiting each other vertex once, followed by its return to 0.
    width = n - 1
    end = 1 << width
    all_vertices = end - 1
    dp = array("q", [infinity]) * (end * width)
    for vertex in range(width):
        dp[(1 << vertex) * width + vertex] = distances[0][vertex + 1]
    for mask in range(1, all_vertices):
        endpoints = mask
        while endpoints:
            endpoint = endpoints & -endpoints
            endpoints ^= endpoint
            vertex = endpoint.bit_length() - 1
            cost = dp[mask * width + vertex]
            row = distances[vertex + 1]
            missing = all_vertices ^ mask
            while missing:
                bit = missing & -missing
                missing ^= bit
                next_vertex = bit.bit_length() - 1
                index = (mask | bit) * width + next_vertex
                candidate = cost + row[next_vertex + 1]
                if candidate < dp[index]:
                    dp[index] = candidate
    return min(dp[all_vertices * width + v] + distances[v + 1][0]
               for v in range(width))


def verify(record):
    n = record["cities"]
    check(type(n) is int and 3 <= n <= 18, "Unsupported city count")
    edges = [(u, v) for u, v, _ in record["edges"]]
    point = [Fraction(x) for _, _, x in record["edges"]]
    check(all(type(u) is int and type(v) is int and 0 <= u < v < n for u, v in edges),
          "Invalid endpoint")
    check(len(set(edges)) == len(edges) and all(0 < x <= 1 for x in point), "Invalid support")
    check_subtour(n, edges, point)
    check_extreme(n, edges, point, record.get("tight_cuts", []))
    lower, upper = (Fraction(record[k]) for k in ("lower", "upper"))
    check(0 < lower <= upper, "Inconsistent interval")
    check_mixture(n, edges, point, record["mixture"], upper)
    lengths = record["lengths"]
    check(len(lengths) == len(edges) and all(type(w) is int and w >= 0 for w in lengths),
          "Invalid edge lengths")
    cost = sum(x * w for x, w in zip(point, lengths, strict=True))
    check(cost > 0 and minimum_tour(n, edges, lengths) == cost * lower,
          "Incorrect tour lower bound")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger", type=Path)
    args = parser.parse_args()
    count = 0
    seen = set()
    with args.ledger.open(encoding="utf-8") as ledger:
        for line in ledger:
            record = json.loads(line)
            verify(record)
            key = tuple(sorted((u, v, Fraction(x)) for u, v, x in record["edges"]))
            check(key not in seen, "Duplicate vector")
            seen.add(key)
            count += 1
            print(f"{record['case']}: certified {record['lower']} to {record['upper']}", flush=True)
    check(count > 0, "Empty ledger")
    print(f"Verified {count} distinct subtour vectors.")

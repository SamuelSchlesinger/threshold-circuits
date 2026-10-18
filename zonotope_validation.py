#!/usr/bin/env python3
"""Finite checks for the finite-state and counting sections (Sections 6, 9, 11).

Checked statements:

* Lemma 9.2 (integer additive images): |S(R)| is at most the sum of lattice
  indices over independent column sets (computed as gcds of maximal minors),
  at most the Hadamard-weighted binomial sum, and at most (e k W / sqrt r)^r.
* Proposition 9.5 (totally unimodular summaries): |S(R)| equals the number of
  independent column subsets, which is at most sum_{j <= r} C(k, j).
* The forest count for reduced graph incidence matrices.
* The Fibonacci identity (P1) for s = 1..100.
* Theorem 9.3: the multiplicity-preserving dictionary algorithm returns the
  exact satisfying-input count of rational threshold circuits whose input
  rows lie in a supplied integer row span, compared with exhaustive
  enumeration, and every summary class has one transcript.
* Theorem 6.1: the pending-bias update reproduces circuit evaluation, with
  state counts within the bounds (9) and (11).
* Proposition 6.3: the width-two path encoding of 3-CNF preserves the
  satisfying-assignment count.
* Theorem 11.1: the quadratic potential has a unique sink on every face of
  the gate cube, equal to the circuit transcript on the full cube.

These are semantic checks of small instances, not proofs for all sizes.
"""

import math
import random
from fractions import Fraction
from itertools import combinations, product

SEED = 20261018


# --------------------------------------------------------------------------
# Exact linear algebra over the rationals
# --------------------------------------------------------------------------

def rank(rows):
    M = [[Fraction(v) for v in row] for row in rows]
    r = 0
    cols = len(M[0]) if M else 0
    for c in range(cols):
        pivot = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if pivot is None:
            continue
        M[r], M[pivot] = M[pivot], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                factor = M[i][c] / M[r][c]
                M[i] = [a - factor * b for a, b in zip(M[i], M[r])]
        r += 1
        if r == len(M):
            break
    return r


def det(M):
    n = len(M)
    A = [[Fraction(v) for v in row] for row in M]
    d = Fraction(1)
    for c in range(n):
        pivot = next((i for i in range(c, n) if A[i][c] != 0), None)
        if pivot is None:
            return Fraction(0)
        if pivot != c:
            A[c], A[pivot] = A[pivot], A[c]
            d = -d
        d *= A[c][c]
        for i in range(c + 1, n):
            if A[i][c] != 0:
                factor = A[i][c] / A[c][c]
                A[i] = [a - factor * b for a, b in zip(A[i], A[c])]
    return d


def columns(R, J):
    return [[row[j] for j in J] for row in R]


def boolean_image(R):
    k = len(R[0])
    return {tuple(sum(row[i] * u[i] for i in range(k)) for row in R)
            for u in product((0, 1), repeat=k)}


def independent_sets(R):
    r, k = len(R), len(R[0])
    sets = []
    for j in range(0, min(r, k) + 1):
        for J in combinations(range(k), j):
            if j == 0 or rank(columns(R, J)) == j:
                sets.append(J)
    return sets


def lattice_index(R, J):
    """gcd of the maximal minors of R_J: the index [span(R_J) cap Z^r : R_J Z^j]."""
    j = len(J)
    if j == 0:
        return 1
    g = 0
    for I in combinations(range(len(R)), j):
        minor = det([[R[i][c] for c in J] for i in I])
        assert minor.denominator == 1
        g = math.gcd(g, abs(int(minor)))
    assert g > 0
    return g


def is_totally_unimodular(R):
    r, k = len(R), len(R[0])
    for j in range(1, min(r, k) + 1):
        for I in combinations(range(r), j):
            for J in combinations(range(k), j):
                if abs(det([[R[i][c] for c in J] for i in I])) > 1:
                    return False
    return True


# --------------------------------------------------------------------------
# Lemma 9.2 and Proposition 9.5
# --------------------------------------------------------------------------

def random_full_rank(rng, r, k, W):
    while True:
        R = [[rng.randint(-W, W) for _ in range(k)] for _ in range(r)]
        if rank(R) == r:
            return R


def interval_matrix(rng, r, k):
    """A consecutive-ones matrix (totally unimodular)."""
    R = []
    for _ in range(r):
        a = rng.randrange(k)
        b = rng.randrange(a, k)
        R.append([1 if a <= j <= b else 0 for j in range(k)])
    return R


def incidence_matrix(rng, v, k):
    """A vertex-edge incidence matrix of a random oriented multigraph."""
    edges = [tuple(rng.sample(range(v), 2)) for _ in range(k)]
    R = [[0] * k for _ in range(v)]
    for j, (a, b) in enumerate(edges):
        R[a][j] = 1
        R[b][j] = -1
    return R, edges


def row_basis(R):
    basis = []
    for row in R:
        if rank(basis + [row]) > len(basis):
            basis.append(row)
    return basis


def check_image_bounds(rng):
    matrices = 0
    tu_cases = 0
    for _ in range(400):
        r = rng.randint(1, 4)
        k = rng.randint(r, 9)
        W = rng.randint(1, 3)
        kind = rng.random()
        if kind < 0.6:
            R = random_full_rank(rng, r, k, W)
        elif kind < 0.8:
            R = row_basis(interval_matrix(rng, r, k))
            if not R:
                continue
            r, W = len(R), 1
        else:
            R, _ = incidence_matrix(rng, r + 1, k)
            R = row_basis(R)
            r, W = len(R), 1
        image = boolean_image(R)
        indep = independent_sets(R)
        index_sum = sum(lattice_index(R, J) for J in indep)
        hadamard = sum(math.comb(k, j) * (math.sqrt(j) * W) ** j for j in range(r + 1))
        simple = (math.e * k * W / math.sqrt(r)) ** r
        assert len(image) <= 2 ** k
        assert len(image) <= index_sum, (R, len(image), index_sum)
        assert index_sum <= hadamard * (1 + 1e-9), (R, index_sum, hadamard)
        assert hadamard <= simple * (1 + 1e-9), (R, hadamard, simple)
        matrices += 1
        if k <= 8 and is_totally_unimodular(R):
            tu_cases += 1
            assert len(image) == len(indep), (R, len(image), len(indep))
            assert len(indep) <= sum(math.comb(k, j) for j in range(r + 1))
    return matrices, tu_cases


def count_forests(v, edges):
    forests = 0
    for mask in range(1 << len(edges)):
        parent = list(range(v))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        acyclic = True
        for j, (a, b) in enumerate(edges):
            if mask >> j & 1:
                ra, rb = find(a), find(b)
                if ra == rb:
                    acyclic = False
                    break
                parent[ra] = rb
        forests += acyclic
    return forests


def check_forests(rng):
    graphs = 0
    for _ in range(150):
        v = rng.randint(2, 6)
        k = rng.randint(1, 9)
        R, edges = incidence_matrix(rng, v, k)
        basis = row_basis(R)
        if not basis:
            continue
        assert len(boolean_image(basis)) == count_forests(v, edges)
        graphs += 1
    return graphs


def fibonacci(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def check_fibonacci():
    for s in range(1, 101):
        total = 2 + sum(math.comb(k, j) for k in range(s) for j in range(s - k + 1))
        assert total == fibonacci(s + 3), s
    return 100


# --------------------------------------------------------------------------
# Threshold circuits with rational weights
# --------------------------------------------------------------------------

def random_fraction(rng, lo=-3, hi=3, den=3):
    return Fraction(rng.randint(lo * den, hi * den), rng.randint(1, den))


def evaluate(gates, x):
    z = []
    for a, b, theta in gates:
        score = sum(ai * xi for ai, xi in zip(a, x)) + sum(bj * zj for bj, zj in zip(b, z))
        z.append(1 if score >= theta else 0)
    return z


def random_summary_circuit(rng, N, R, s):
    gates = []
    for i in range(s):
        coeffs = [random_fraction(rng) for _ in R]
        if rng.random() < 0.25:
            coeffs = [Fraction(0)] * len(R)
        a = [sum(c * row[j] for c, row in zip(coeffs, R)) for j in range(N)]
        b = [random_fraction(rng) if rng.random() < 0.7 else Fraction(0) for _ in range(i)]
        gates.append((a, b, random_fraction(rng)))
    return gates


def dictionary_count(R, gates, N):
    r = len(R)
    table = {tuple([0] * r): (1, ())}
    stages = 0
    for j in range(N):
        new = {}
        col = tuple(row[j] for row in R)
        for key, (count, rep) in table.items():
            for bit in (0, 1):
                nk = tuple(k + bit * c for k, c in zip(key, col))
                if nk in new:
                    new[nk] = (new[nk][0] + count, new[nk][1])
                else:
                    new[nk] = (count, rep + (bit,))
        table = new
        stages = max(stages, len(table))
    accepted = sum(count for key, (count, rep) in table.items()
                   if evaluate(gates, rep)[-1] == 1)
    return accepted, table, stages


def check_counting(rng):
    circuits = 0
    nonconstant = 0
    for _ in range(100):
        N = rng.randint(2, 8)
        r = rng.randint(1, min(3, N))
        W = rng.randint(1, 2)
        R = random_full_rank(rng, r, N, W)
        s = rng.randint(1, 6)
        gates = random_summary_circuit(rng, N, R, s)
        accepted, table, stages = dictionary_count(R, gates, N)
        by_key = {}
        truth = 0
        for x in product((0, 1), repeat=N):
            key = tuple(sum(row[j] * x[j] for j in range(N)) for row in R)
            z = tuple(evaluate(gates, x))
            by_key.setdefault(key, []).append(z)
            truth += z[-1]
        assert set(by_key) == set(table)
        for key, transcripts in by_key.items():
            assert len(set(transcripts)) == 1, "equal summaries must share a transcript"
            assert table[key][0] == len(transcripts)
        assert accepted == truth, (R, gates, accepted, truth)
        assert stages <= len(boolean_image(R))
        circuits += 1
        nonconstant += 0 < truth < 2 ** N
    return circuits, nonconstant


# --------------------------------------------------------------------------
# Pending biases, the 3-CNF path encoding, and the quadratic potential
# --------------------------------------------------------------------------

def check_pending_bias(rng):
    circuits = 0
    for _ in range(200):
        N = rng.randint(1, 5)
        s = rng.randint(1, 6)
        W = rng.randint(1, 3)
        gates = []
        for i in range(s):
            a = [rng.randint(-2, 2) for _ in range(N)]
            b = [rng.randint(-W, W) for _ in range(i)]
            gates.append((a, b, rng.randint(-3, 3)))
        states = [set() for _ in range(s + 1)]
        for x in product((0, 1), repeat=N):
            z = evaluate(gates, x)
            # Simulate (8): the pending vector p_k over gates k+1..s.
            p = [0] * s
            for k in range(s):
                a, b, theta = gates[k]
                zk = 1 if sum(ai * xi for ai, xi in zip(a, x)) >= theta - p[k] else 0
                assert zk == z[k]
                for later in range(k + 1, s):
                    p[later] += gates[later][1][k] * zk
                states[k + 1].add(tuple(p[k + 1:]))
        for k in range(1, s + 1):
            rows = [gates[later][1][:k] for later in range(k, s)]
            r_k = rank(rows) if rows else 0
            assert len(states[k]) <= min(2 ** k, (k * W + 1) ** r_k)
            if rows:
                cols = {tuple(row[j] for row in rows) for j in range(k)}
                mult = {c: sum(1 for j in range(k) if tuple(row[j] for row in rows) == c)
                        for c in cols}
                assert len(states[k]) <= math.prod(m + 1 for m in mult.values())
        circuits += 1
    return circuits


def check_cnf_paths(rng):
    formulas = 0
    for _ in range(100):
        n = rng.randint(1, 7)
        c = rng.randint(1, 6)
        clauses = [[(rng.randrange(n), rng.randrange(2)) for _ in range(3)] for _ in range(c)]
        gates = []
        for i, clause in enumerate(clauses):
            a = [Fraction(0)] * n
            const = Fraction(0)
            for var, neg in clause:
                if neg:
                    a[var] -= Fraction(1, 3)
                    const += Fraction(1, 3)
                else:
                    a[var] += Fraction(1, 3)
            b = [Fraction(0)] * i
            if i > 0:
                b[i - 1] = Fraction(1)
            theta = (Fraction(1, 3) if i == 0 else Fraction(4, 3)) - const
            gates.append((a, b, theta))
        states = [set() for _ in range(c + 1)]
        count = 0
        for x in product((0, 1), repeat=n):
            z = evaluate(gates, x)
            sat = all(any(x[v] != neg for v, neg in clause) for clause in clauses)
            assert z[-1] == sat
            count += sat
            for k in range(1, c + 1):
                states[k].add(tuple(sum(gates[later][1][j] * z[j] for j in range(k))
                                    for later in range(k, c)))
        assert all(len(S) <= 2 for S in states)
        formulas += 1
    return formulas


def check_potential(rng):
    checked = 0
    for _ in range(60):
        N = rng.randint(1, 3)
        s = rng.randint(1, 4)
        gates = []
        for i in range(s):
            a = [rng.randint(-2, 2) for _ in range(N)]
            b = [rng.randint(-2, 2) for _ in range(i)]
            gates.append((a, b, rng.randint(-2, 2)))
        # Strict integer margins: Q = 2(score - theta) + 1 is odd, |Q| >= 1.
        lam = [0] * s
        for i in reversed(range(s)):
            lam[i] = 1 + sum(lam[j] * abs(2 * gates[j][1][i]) for j in range(i + 1, s))

        def Q(i, x, z):
            a, b, theta = gates[i]
            score = sum(ai * xi for ai, xi in zip(a, x)) + sum(b[j] * z[j] for j in range(i))
            return 2 * (score - theta) + 1

        def energy(x, z):
            return -sum(lam[i] * z[i] * Q(i, x, z) for i in range(s))

        for x in product((0, 1), repeat=N):
            transcript = tuple(evaluate(gates, x))
            for face in product((0, 1, None), repeat=s):
                free = [i for i in range(s) if face[i] is None]
                points = []
                for bits in product((0, 1), repeat=len(free)):
                    z = list(face)
                    for i, bit in zip(free, bits):
                        z[i] = bit
                    points.append(tuple(z))
                sinks = []
                for z in points:
                    E = energy(x, z)
                    is_sink = True
                    for i in free:
                        w = list(z)
                        w[i] ^= 1
                        E2 = energy(x, tuple(w))
                        assert E2 != E, "tied edge"
                        if E2 < E:
                            is_sink = False
                            break
                    if is_sink:
                        sinks.append(z)
                assert len(sinks) == 1, (gates, x, face, sinks)
                if not any(v is not None for v in face):
                    assert sinks[0] == transcript
            checked += 1
    return checked


def main():
    rng = random.Random(SEED)
    matrices, tu_cases = check_image_bounds(rng)
    print(f"integer additive images: {matrices} full-rank matrices, "
          f"including {tu_cases} totally unimodular cases")
    graphs = check_forests(rng)
    print(f"graph incidence summaries: {graphs} oriented graphs, image size = forest count")
    print(f"Fibonacci identity (P1): s = 1..{check_fibonacci()}")
    circuits, nonconstant = check_counting(rng)
    print(f"exact counting: {circuits} rational circuits with supplied summaries, "
          f"{nonconstant} with nonconstant output, all counts match enumeration")
    print(f"pending-bias program: {check_pending_bias(rng)} circuits simulated")
    print(f"width-two 3-CNF encodings: {check_cnf_paths(rng)} formulas")
    print(f"quadratic unique-sink potential: {check_potential(rng)} (circuit, input) pairs, "
          f"every face checked")
    print("zonotope_validation: all checks passed")


if __name__ == "__main__":
    main()

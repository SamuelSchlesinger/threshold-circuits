#!/usr/bin/env python3
"""Finite checks for the leaf-count rectangle theorem and the energy section.

The checks exercise the certificate extracted by the proof of the leaf-count
rectangle theorem (Theorem 2.2) and its bi-monotone extension (Theorem 2.6):
at every node the quantile pivot keeps at least the chosen child's leaf
fraction on both sides, and the fractions telescope to 1/L along the selected
path.  Queries are separable comparisons alpha(x) + beta(y) >= theta or the
bi-monotone variants max(alpha, beta) >= theta and min(alpha, beta) >= theta.  They also check the
Hadamard rectangle bound for inner product on small sizes, the depth-three
block construction of Section 3.3, and the activation-pattern bound
P(C) <= B(q, e) on random small threshold circuits.

It also checks the affine leaf scores of Theorem 2.9 (the reached leaf beats
every other leaf by at least one, in exact integer arithmetic) and the
Hamming-ball lists of Proposition 3.7 (exact mass, bias, facet description,
and list agreement against the stated bounds).

These are semantic checks of finite instances.  They do not prove the
theorems for unbounded parameters.
"""

import math
import random
from fractions import Fraction
from itertools import product

SEED = 42
TREES = 10_000
MAX_LEAVES = 39
MAX_DOMAIN = 49
SCORE_RANGE = 6  # small integer scores produce many ties


# --------------------------------------------------------------------------
# Separable-comparison trees
# --------------------------------------------------------------------------

class Leaf:
    __slots__ = ("label", "leaves", "ident")

    def __init__(self, label, ident):
        self.label = label
        self.leaves = 1
        self.ident = ident


QUERY_KINDS = ("sum", "sum", "max", "min")


class Node:
    __slots__ = ("alpha", "beta", "theta", "kind", "zero", "one", "leaves")

    def __init__(self, alpha, beta, theta, kind, zero, one):
        self.alpha = alpha
        self.beta = beta
        self.theta = theta
        self.kind = kind
        self.zero = zero
        self.one = one
        self.leaves = zero.leaves + one.leaves

    def query(self, x, y):
        a, b = self.alpha[x], self.beta[y]
        score = a + b if self.kind == "sum" else (max(a, b) if self.kind == "max" else min(a, b))
        return 1 if score >= self.theta else 0


def random_tree(rng, leaves, nx, ny, counter):
    """A random binary separable-comparison tree with exactly `leaves` leaves."""
    if leaves == 1:
        counter[0] += 1
        return Leaf(rng.randrange(2), counter[0])
    left = rng.randrange(1, leaves)
    alpha = [rng.randrange(-SCORE_RANGE, SCORE_RANGE + 1) for _ in range(nx)]
    beta = [rng.randrange(-SCORE_RANGE, SCORE_RANGE + 1) for _ in range(ny)]
    theta = rng.randrange(-2 * SCORE_RANGE, 2 * SCORE_RANGE + 1)
    kind = rng.choice(QUERY_KINDS)
    return Node(alpha, beta, theta, kind,
                random_tree(rng, left, nx, ny, counter),
                random_tree(rng, leaves - left, nx, ny, counter))


def route(tree, x, y):
    while isinstance(tree, Node):
        tree = tree.one if tree.query(x, y) else tree.zero
    return tree.ident


def extract_rectangle(tree, A, B):
    """The rectangle supplied by the proof of the leaf-count theorem.

    Returns (A', B', leaf identifier).  Both lists are returned sorted by
    the current node's scores, which the recursion re-sorts as needed.
    """
    if isinstance(tree, Leaf):
        return A, B, tree.ident
    L = tree.leaves
    L0 = tree.zero.leaves
    u, v = len(A), len(B)
    # i = ceil(a u) with a = L0 / L, computed exactly in integers.
    i = -((-L0 * u) // L)
    j = -((-L0 * v) // L)
    A_sorted = sorted(A, key=lambda x: tree.alpha[x])
    B_sorted = sorted(B, key=lambda y: tree.beta[y])
    xi, yj = A_sorted[i - 1], B_sorted[j - 1]
    if tree.query(xi, yj) == 0:
        A0, B0 = A_sorted[:i], B_sorted[:j]
        # Quantile lemma, 0-branch: sizes at least a|A| and a|B|.
        assert A0 and B0
        assert len(A0) * L >= L0 * u and len(B0) * L >= L0 * v
        for x in A0:
            for y in B0:
                assert tree.query(x, y) == 0
        return extract_rectangle(tree.zero, A0, B0)
    A1, B1 = A_sorted[i - 1:], B_sorted[j - 1:]
    L1 = L - L0
    assert len(A1) * L >= L1 * u and len(B1) * L >= L1 * v
    for x in A1:
        for y in B1:
            assert tree.query(x, y) == 1
    return extract_rectangle(tree.one, A1, B1)


def check_leaf_count_theorem():
    rng = random.Random(SEED)
    entries = 0
    for _ in range(TREES):
        leaves = rng.randrange(1, MAX_LEAVES + 1)
        nx = rng.randrange(1, MAX_DOMAIN + 1)
        ny = rng.randrange(1, MAX_DOMAIN + 1)
        tree = random_tree(rng, leaves, nx, ny, [0])
        A, B = list(range(nx)), list(range(ny))
        A2, B2, ident = extract_rectangle(tree, A, B)
        L = tree.leaves
        assert len(A2) * L >= nx, (len(A2), L, nx)
        assert len(B2) * L >= ny, (len(B2), L, ny)
        for x in A2:
            for y in B2:
                assert route(tree, x, y) == ident
        entries += len(A2) * len(B2)
    return entries


# --------------------------------------------------------------------------
# Inner product: Hadamard rectangle bound and the linear decision tree bound
# --------------------------------------------------------------------------

def inner_product(x, y):
    return bin(x & y).count("1") & 1


def check_hadamard_rectangles(max_m=3):
    """Every monochromatic rectangle for IP_m has at most 2^m entries (m <= 3)."""
    rectangles = 0
    for m in range(1, max_m + 1):
        M = 1 << m
        rows = list(range(M))
        for amask in range(1, 1 << M):
            A = [r for r in rows if amask >> r & 1]
            # The largest B on which A x B is monochromatic, for each color.
            for color in (0, 1):
                B = [c for c in rows if all(inner_product(a, c) == color for a in A)]
                if B:
                    rectangles += 1
                    assert len(A) * len(B) <= M, (m, A, B)
    return rectangles


def all_ldt_leaf_bound(m=2):
    """Exhaustively check, for IP_2, that the fraction (2) forces L >= 2^{m/2}.

    This is not an exhaustive search over trees (infinitely many weights);
    it checks the corollary's arithmetic: with |X| = |Y| = 2^m and rectangle
    area at most 2^m, the rectangle theorem forces L^2 >= 2^m.
    """
    M = 1 << m
    L = 1
    while (M // L if M % L == 0 else M / L) ** 2 > M:
        L += 1
    return L


# --------------------------------------------------------------------------
# Threshold circuits: activation patterns and the block construction
# --------------------------------------------------------------------------

def evaluate_circuit(gates, x):
    """gates: list of (input_weights, gate_weights, theta).  Returns transcript."""
    z = []
    for a, b, theta in gates:
        score = sum(ai * xi for ai, xi in zip(a, x)) + sum(bj * zj for bj, zj in zip(b, z))
        z.append(1 if score >= theta else 0)
    return z


def random_circuit(rng, n, s):
    gates = []
    for i in range(s):
        a = [rng.choice([0, 0, 1, -1, 2, -2]) for _ in range(n)]
        if rng.random() < 0.3:
            a = [0] * n
        b = [rng.choice([0, 0, 1, -1, 2]) for _ in range(i)]
        theta = rng.randrange(-3, 4)
        gates.append((a, b, theta))
    return gates


def check_pattern_bound(rng, trials=300):
    """P(C) <= B(q, e(C)) for the input-reading transcript of random circuits."""
    checked = 0
    for _ in range(trials):
        n = rng.randrange(1, 7)
        s = rng.randrange(1, 8)
        gates = random_circuit(rng, n, s)
        reading = [i for i, (a, _, _) in enumerate(gates) if any(a)]
        q = len(reading)
        patterns = set()
        e = 0
        for x in product((0, 1), repeat=n):
            z = evaluate_circuit(gates, x)
            restricted = tuple(z[i] for i in reading)
            patterns.add(restricted)
            e = max(e, sum(restricted))
        ball = sum(math.comb(q, j) for j in range(e + 1))
        assert len(patterns) <= ball, (gates, patterns, e)
        assert len(patterns) <= 2 ** n
        checked += 1
    return checked


def block_circuit(m, k):
    """The depth-three construction of Section 3.3 for IP_m with block size k.

    Returns gates over 2m inputs ordered x_1..x_m, y_1..y_m, with the output
    gate last.  Weights are in {-1, 0, 1} and thresholds are integers.
    """
    blocks = [list(range(i, min(i + k, m))) for i in range(0, m, k)]
    gates = []
    minterm_indices = []
    for block in blocks:
        kk = len(block)
        for local in product((0, 1), repeat=2 * kk):
            xs, ys = local[:kk], local[kk:]
            if sum(a & b for a, b in zip(xs, ys)) & 1 == 0:
                continue
            a = [0] * (2 * m)
            ones = 0
            for idx, coord in enumerate(block):
                for side, bit in ((0, xs[idx]), (m, ys[idx])):
                    a[side + coord] = 1 if bit else -1
                    ones += bit
            # sum of literals >= 2kk  <=>  sum_{bit=1} x - sum_{bit=0} x >= ones
            gates.append((a, [0] * len(gates), ones))
            minterm_indices.append(len(gates) - 1)
    b = len(blocks)
    base = len(gates)
    for j in range(1, b + 1):
        w = [0] * base
        for idx in minterm_indices:
            w[idx] = 1
        gates.append(([0] * (2 * m), w, j))
    w = [0] * len(gates)
    for j in range(1, b + 1):
        w[base + j - 1] = 1 if j % 2 == 1 else -1
    gates.append(([0] * (2 * m), w, 1))
    return gates, len(minterm_indices), b


def check_block_construction(max_m=5):
    circuits = 0
    for m in range(1, max_m + 1):
        for k in range(1, m + 1):
            gates, q, b = block_circuit(m, k)
            assert len(gates) == q + b + 1
            sizes = [min(k, m - i) for i in range(0, m, k)]
            assert q == sum((4 ** kk - 2 ** kk) // 2 for kk in sizes)
            max_energy = 0
            max_reading = 0
            for x in product((0, 1), repeat=2 * m):
                z = evaluate_circuit(gates, x)
                ip = sum(x[i] & x[m + i] for i in range(m)) & 1
                assert z[-1] == ip, (m, k, x, z)
                max_energy = max(max_energy, sum(z))
                max_reading = max(max_reading, sum(z[:q]))
            assert max_reading <= b
            assert max_energy <= 2 * b + 1
            circuits += 1
    return circuits


# --------------------------------------------------------------------------
# Theorem 2.9: affine leaf scores
# --------------------------------------------------------------------------

def leaves_with_paths(tree, depth=0, path=()):
    if isinstance(tree, Leaf):
        yield tree.ident, path
        return
    yield from leaves_with_paths(tree.zero, depth + 1, path + ((tree, depth, 0),))
    yield from leaves_with_paths(tree.one, depth + 1, path + ((tree, depth, 1),))


def check_affine_scores(rng, trials=400):
    """Scores S_p = sum over one-edges of B^(d-1-depth) q_v, with q_v = 2(alpha+beta-theta)+1."""
    checked = 0
    for _ in range(trials):
        leaves = rng.randrange(1, 13)
        nx, ny = rng.randrange(1, 9), rng.randrange(1, 9)
        tree = random_tree(rng, leaves, nx, ny, [0])
        # Restrict to separable (sum) queries, which are affine.
        def make_sum(t):
            if isinstance(t, Node):
                t.kind = "sum"
                make_sum(t.zero); make_sum(t.one)
        make_sum(tree)
        paths = list(leaves_with_paths(tree))
        d = max((len(path) for _, path in paths), default=0)
        if d == 0:
            continue
        def q(node, x, y):
            return 2 * (node.alpha[x] + node.beta[y] - node.theta) + 1
        M = 1
        for x in range(nx):
            for y in range(ny):
                for _, path in paths:
                    for node, _, _ in path:
                        M = max(M, abs(q(node, x, y)))
        B = 2 * M + 1
        for x in range(nx):
            for y in range(ny):
                scores = {}
                for ident, path in paths:
                    scores[ident] = sum(B ** (d - 1 - dep) * q(node, x, y)
                                        for node, dep, edge in path if edge == 1)
                    for node, _, _ in path:
                        assert abs(q(node, x, y)) >= 1 and (q(node, x, y) > 0) == bool(node.query(x, y))
                reached = route(tree, x, y)
                for ident, sc in scores.items():
                    if ident != reached:
                        assert scores[reached] >= sc + 1, (x, y, reached, ident)
                checked += 1
    return checked


# --------------------------------------------------------------------------
# Proposition 3.7: Hamming-ball cells and lists
# --------------------------------------------------------------------------

def check_ball_lists(max_m=8):
    from itertools import combinations
    cases = 0
    for m in range(2, max_m + 1):
        pts = list(product((0, 1), repeat=2 * m))
        total = len(pts)
        for k in range(0, (m - 1) // 4 + 1):
            if not 4 * k < m:
                continue
            facets = list(combinations(range(m), k + 1))
            r = Fraction(3 * k, m - k + 1)
            P_k = Fraction(math.comb(m, k) * 3 ** (m - k), 4 ** m)
            mass = 0
            signed = 0
            agreement = 0
            cell_counts = {}
            for pt in pts:
                x, y = pt[:m], pt[m:]
                inter = sum(a & b for a, b in zip(x, y))
                ip = inter & 1
                inside = all(sum(x[i] + y[i] for i in U) <= 2 * k + 1 for U in facets)
                assert inside == (inter <= k)
                if inside:
                    mass += 1
                    signed += 1 - 2 * ip
                    agreement += ip == (k % 2)
                else:
                    first = next(j for j, U in enumerate(facets)
                                 if sum(x[i] + y[i] for i in U) >= 2 * k + 2)
                    cell_counts.setdefault(first, [0, 0])[ip] += 1
            agreement += sum(max(c) for c in cell_counts.values())
            assert Fraction(mass, total) >= P_k
            assert Fraction(abs(signed), mass) >= (1 - r) ** 2
            assert Fraction(agreement, total) >= Fraction(1, 2) + (1 - r) ** 2 * P_k / 2
            cases += 1
    return cases


# --------------------------------------------------------------------------
# Analytic constants used in the energy corollaries
# --------------------------------------------------------------------------

def h2(p):
    if p in (0, 1):
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def h2_inverse(y, lo=0.0, hi=0.5):
    for _ in range(200):
        mid = (lo + hi) / 2
        if h2(mid) < y:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def check_constants():
    quarter = h2_inverse(0.25)
    half = h2_inverse(0.5)
    assert abs(quarter - 0.04169269027) < 1e-10, quarter
    assert abs(half - 0.1100278644) < 1e-9, half
    # Hamming-ball estimate (7) on a grid, and the (5) form.
    checks = 0
    for v in range(1, 61):
        for u in range(0, v // 2 + 1):
            ball = sum(math.comb(v, j) for j in range(u + 1))
            assert math.log2(ball) <= v * h2(Fraction(u, v)) + 1e-9
            checks += 1
        for e in range(1, v + 1):
            ball = sum(math.comb(v, j) for j in range(e + 1))
            assert ball <= (math.e * v / e) ** e * (1 + 1e-12)
            checks += 1
    return quarter, half, checks


def main():
    entries = check_leaf_count_theorem()
    print(f"leaf-count rectangle theorem: {TREES} random separable and bi-monotone "
          f"trees, seed {SEED}, {entries} rectangle entries routed to a single leaf")
    rect = check_hadamard_rectangles()
    print(f"Hadamard rectangle bound: {rect} maximal monochromatic rectangles "
          f"checked for m <= 3")
    rng = random.Random(SEED)
    patterns = check_pattern_bound(rng)
    print(f"activation-pattern Hamming-ball bound: {patterns} random circuits")
    circuits = check_block_construction()
    print(f"depth-three block construction for IP_m: {circuits} (m, k) pairs, "
          f"m <= 5, exhaustive inputs")
    scores = check_affine_scores(rng)
    print(f"affine leaf scores: {scores} (tree, input) pairs, reached leaf wins by at least one")
    balls = check_ball_lists()
    print(f"Hamming-ball cells and lists: {balls} (m, k) pairs with m <= 8, exact masses, "
          f"biases, and list agreements")
    quarter, half, checks = check_constants()
    print(f"constants: h2^-1(1/4) = {quarter:.11f}, h2^-1(1/2) = {half:.10f}, "
          f"{checks} Hamming-ball estimates")
    print("weighted_pivot_validation: all checks passed")


if __name__ == "__main__":
    main()

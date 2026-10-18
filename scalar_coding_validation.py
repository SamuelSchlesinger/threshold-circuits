#!/usr/bin/env python3
"""Finite checks for the exact scalar model (Section 4).

Checked statements:

* Theorem 4.1, upper direction: for q <= 6 and every L between 1 and B(q, e),
  the recursive construction yields an ordered comparison tree with exactly
  L leaves, depth at most q, and at most e right branches on every path; the
  one-hot prefix indicator gate p_h recognizes exactly its prefix.
* Theorem 4.3 (enumerative coding): for every (q, e) with q <= 10 and every
  rank t < B(q, e), the comparison procedure recovers the word of
  lexicographic rank t among words of weight at most e, with the stated
  invariant at every stage.
* Theorem 4.3 as a circuit: for N <= 9 and every feasible (q, e) with
  q <= 12, the gate-level realization (flags u, d, b; the q scalar queries;
  the N decoder gates) recovers every input x with at most e active queries.
* Corollary 4.4: the sharp N-versus-N+1 capacity thresholds for N <= 9.

These are finite semantic checks, not proofs for unbounded parameters.
"""

import math
from itertools import product


def B(q, e):
    if e < 0:
        return 0
    return sum(math.comb(q, j) for j in range(min(q, e) + 1))


# --------------------------------------------------------------------------
# Theorem 4.1: ordered comparison trees with a prescribed leaf count
# --------------------------------------------------------------------------

def build_tree(q, e, L):
    """A tree over L consecutive runs 0..L-1 with depth <= q and <= e right branches."""
    assert 1 <= L <= B(q, e)
    if L == 1:
        return ("leaf", 0)
    assert e >= 1 and q >= 1
    cap0, cap1 = B(q - 1, e), B(q - 1, e - 1)
    L0 = min(cap0, L - 1)
    L1 = L - L0
    assert 1 <= L1 <= cap1
    return ("node", L0, build_tree(q - 1, e, L0), build_tree(q - 1, e - 1, L1))


def tree_paths(tree, depth=0, rights=0, offset=0):
    """Yields (leaf run index, depth, right branches)."""
    if tree[0] == "leaf":
        yield offset, depth, rights
        return
    _, L0, left, right = tree
    yield from tree_paths(left, depth + 1, rights, offset)
    yield from tree_paths(right, depth + 1, rights + 1, offset + L0)


def check_trees(max_q=6):
    trees = 0
    for q in range(0, max_q + 1):
        for e in range(0, q + 1):
            for L in range(1, B(q, e) + 1):
                tree = build_tree(q, e, L)
                leaves = sorted(tree_paths(tree))
                assert [leaf for leaf, _, _ in leaves] == list(range(L))
                assert all(d <= q and r <= e for _, d, r in leaves)
                trees += 1
    return trees


def check_prefix_indicators(max_len=6):
    checks = 0
    for n in range(0, max_len + 1):
        for h in product((0, 1), repeat=n):
            for z in product((0, 1), repeat=n):
                fired = sum((2 * hj - 1) * zj for hj, zj in zip(h, z)) >= sum(h)
                assert fired == (h == z)
                checks += 1
    return checks


# --------------------------------------------------------------------------
# Theorem 4.3: enumerative coding by scalar comparisons
# --------------------------------------------------------------------------

def words_of_weight_at_most(q, e):
    return [w for w in product((0, 1), repeat=q) if sum(w) <= e]


def decode_rank(q, e, t):
    """Runs the comparison procedure; returns (word, comparisons, invariants)."""
    z = []
    r, a = 0, 0
    comparisons = invariants = 0
    for i in range(1, q + 1):
        assert 0 <= r <= e and a <= t < a + B(q - i + 1, e - r)
        invariants += 1
        zi = 1 if t >= a + B(q - i, e - r) else 0
        comparisons += 1
        if r == e:
            assert zi == 0
        z.append(zi)
        a += B(q - i, e - r) * zi
        r += zi
    assert a == t and sum(z) <= e
    invariants += 1
    return tuple(z), comparisons, invariants


def check_ranks(max_q=10):
    pairs = ranks = comparisons = invariants = 0
    for q in range(0, max_q + 1):
        for e in range(0, q + 1):
            words = words_of_weight_at_most(q, e)  # lexicographic order
            assert len(words) == B(q, e)
            for t, expected in enumerate(words):
                z, c, inv = decode_rank(q, e, t)
                assert z == expected, (q, e, t, z, expected)
                ranks += 1
                comparisons += c
                invariants += inv
            pairs += 1
    return pairs, ranks, comparisons, invariants


# --------------------------------------------------------------------------
# Theorem 4.3 at gate level
# --------------------------------------------------------------------------

def gate_level_circuit(N, q, e, x):
    """Evaluates the auxiliary flags, the q scalar queries and the decoders.

    Returns (z, decoded x, number of decoder gates evaluated).
    Every gate below is a threshold gate whose only raw-input term is t.
    """
    t = sum(x[i] << i for i in range(N))
    z = []
    b = {}  # (j, r) -> b_{j,r}
    for i in range(1, q + 1):
        prefix_sum = sum(z)
        u = {0: 1}
        for r in range(1, e + 2):
            u[r] = 1 if prefix_sum >= r else 0
        assert u[e + 1] == 0
        d = {r: (1 if u[r] - u[r + 1] >= 1 else 0) for r in range(e + 1)}
        assert sum(d.values()) == 1 and d[prefix_sum] == 1
        offset = sum(B(q - j, e - r) * b[(j, r)] for j in range(1, i) for r in range(e + 1))
        offset += sum(B(q - i, e - r) * d[r] for r in range(e + 1))
        zi = 1 if t - offset >= 0 else 0
        z.append(zi)
        for r in range(e + 1):
            b[(i, r)] = 1 if zi + d[r] >= 2 else 0
    T = sum(B(q - j, e - r) * b[(j, r)] for j in range(1, q + 1) for r in range(e + 1))
    assert T == t
    decoded = [None] * N
    decoders = 0
    for ell in range(N, 0, -1):
        partial = sum((1 << (k - 1)) * decoded[k - 1] for k in range(ell + 1, N + 1))
        decoded[ell - 1] = 1 if T - partial >= (1 << (ell - 1)) else 0
        decoders += 1
    return z, decoded, decoders


def check_gate_level(max_N=9, max_q=12):
    triples = encodings = decoders = 0
    for N in range(1, max_N + 1):
        for q in range(N, max_q + 1):
            for e in range(1, q + 1):
                if B(q, e) < 2 ** N:
                    continue
                for x in product((0, 1), repeat=N):
                    z, decoded, dec = gate_level_circuit(N, q, e, x)
                    assert decoded == list(x), (N, q, e, x, decoded)
                    assert sum(z) <= e
                    encodings += 1
                    decoders += dec
                triples += 1
    return triples, encodings, decoders


def check_sharp_thresholds(max_N=9):
    for N in range(1, max_N + 1):
        half = -(-N // 2)
        assert B(N + 1, half) >= 2 ** N > B(N + 1, half - 1), N
        for e in range(N):
            assert B(N, e) < 2 ** N
    return max_N


def main():
    trees = check_trees()
    print(f"ordered comparison trees: {trees} (q, e, L) triples with q <= 6")
    print(f"one-hot prefix indicators: {check_prefix_indicators()} (h, z) pairs")
    pairs, ranks, comparisons, invariants = check_ranks()
    print(f"enumerative coding: {pairs} parameter pairs with q <= 10, {ranks} ranks, "
          f"{comparisons} comparisons, {invariants} prefix/final invariants")
    triples, encodings, decoders = check_gate_level()
    print(f"gate-level normalization: {triples} feasible (N, q, e) triples with N <= 9 "
          f"and q <= 12, {encodings} encodings and decodings, {decoders} decoder gates")
    print(f"sharp capacity thresholds: N = 1..{check_sharp_thresholds()}")
    print("scalar_coding_validation: all checks passed")


if __name__ == "__main__":
    main()

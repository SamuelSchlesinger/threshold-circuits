#!/usr/bin/env python3
"""Finite checks for the one-sided entropy arguments (Sections 5 and 5.1).

Checked certificates:

* Ordered inner-product runs (Lemma 5.1): for every column order of {0,1}^m,
  the expected logarithm of the maximal monochromatic run through a fixed
  position is below two, and the anchor estimate Pr[A > 2^k] <= 2^{-k-1}
  holds exactly for every position and k.
* The dyadic certificate of Lemma 5.6: for every sign sequence of length
  1, 2, 4, 8, every interval of absolute average at least eta contains a
  maximal dyadic block of size at least eta*len/8 and average bias at least
  eta/2, and the number of dyadic blocks of that size near any anchor
  position is at most 16/eta + 2.
* Exact tail probabilities of L_eta(u, j) for m <= 5 against the bound
  min(1, 1152 / (eta^4 t)) and the expectation bound 12 + 4 log2(1/eta).
* The elementary inequalities used in the corollaries: the quadratic
  entropy bound, the Hamming-ball entropy estimate, and the arithmetic
  that turns (E8) into (E1).

These are finite semantic checks; they do not prove the lemmas for all m.
"""

import math
import random
from fractions import Fraction
from itertools import permutations

SEED = 20261018
ETAS = [Fraction(k, 8) for k in range(1, 9)]


def popcount_parity(x):
    return bin(x).count("1") & 1


# --------------------------------------------------------------------------
# Lemma 5.1: ordered inner-product runs
# --------------------------------------------------------------------------

def run_lengths(seq):
    """For each position, the length of the maximal constant run containing it,
    plus the left extension A (including the position)."""
    M = len(seq)
    left = [0] * M
    for j in range(M):
        left[j] = left[j - 1] + 1 if j > 0 and seq[j] == seq[j - 1] else 1
    right = [0] * M
    for j in reversed(range(M)):
        right[j] = right[j + 1] + 1 if j + 1 < M and seq[j] == seq[j + 1] else 1
    return [left[j] + right[j] - 1 for j in range(M)], left


def check_order(m, order):
    """Returns the conditional run entropy E_{u,j} log2 R(u,j) for one order."""
    M = 1 << m
    total = 0.0
    per_position = [0.0] * M
    anchor_counts = [[0] * (m + 1) for _ in range(M)]
    for u in range(M):
        seq = [popcount_parity(u & v) for v in order]
        R, A = run_lengths(seq)
        for j in range(M):
            lg = math.log2(R[j])
            total += lg
            per_position[j] += lg
            for k in range(m):
                if A[j] > (1 << k):
                    anchor_counts[j][k] += 1
    anchors = 0
    for j in range(M):
        assert per_position[j] / M < 2, (m, order, j, per_position[j] / M)
        for k in range(m):
            # Pr[A > 2^k] <= 2^{-k-1}, exactly in integers.
            assert anchor_counts[j][k] * (1 << (k + 1)) <= M
            anchors += 1
    return total / (M * M), anchors


def check_ordered_runs(rng):
    orders_checked = 0
    anchors = 0
    worst = 0.0
    for m in (1, 2):
        for order in permutations(range(1 << m)):
            value, a = check_order(m, list(order))
            worst = max(worst, value)
            anchors += a
            orders_checked += 1
    for m in range(3, 8):
        for _ in range(160):
            order = list(range(1 << m))
            rng.shuffle(order)
            value, a = check_order(m, order)
            worst = max(worst, value)
            anchors += a
            orders_checked += 1
    return orders_checked, anchors, worst


# --------------------------------------------------------------------------
# Lemma 5.6: dyadic certificate
# --------------------------------------------------------------------------

def maximal_dyadic_blocks(lo, hi):
    """Maximal aligned dyadic blocks partitioning the positions lo..hi-1 (0-indexed)."""
    blocks = []
    a = lo
    while a < hi:
        size = 1
        while a % (2 * size) == 0 and a + 2 * size <= hi:
            size *= 2
        blocks.append((a, size))
        a += size
    return blocks


def aligned_blocks_meeting(lo, hi, h):
    """Number of aligned blocks of size h (on the integer line) meeting [lo, hi]."""
    first = math.floor(lo / h)
    last = math.floor(hi / h)
    return last - first + 1


def check_dyadic_certificates():
    certificates = 0
    locations = 0
    for M in (1, 2, 4, 8):
        for bits in range(1 << M):
            X = [1 if bits >> i & 1 else -1 for i in range(M)]
            prefix = [0]
            for x in X:
                prefix.append(prefix[-1] + x)
            for lo in range(M):
                for hi in range(lo + 1, M + 1):
                    ell = hi - lo
                    total = abs(prefix[hi] - prefix[lo])
                    blocks = maximal_dyadic_blocks(lo, hi)
                    sizes = [s for _, s in blocks]
                    assert all(sizes.count(s) <= 2 for s in set(sizes))
                    for eta in ETAS:
                        if total < eta * ell:
                            continue
                        retained = [(a, s) for a, s in blocks if s >= eta * ell / 8]
                        discarded = sum(s for a, s in blocks if s < eta * ell / 8)
                        assert discarded < eta * ell / 2
                        witnesses = [(a, s) for a, s in retained
                                     if abs(prefix[a + s] - prefix[a]) >= (eta / 2) * s]
                        assert witnesses, (X, lo, hi, eta)
                        certificates += 1
                        a, h = witnesses[0]
                        for j in range(lo, hi):
                            assert all(abs(p - j) <= ell for p in range(a, a + h))
                            count = aligned_blocks_meeting(j - ell, j + ell, h)
                            assert count <= 16 / eta + 2, (count, eta, h, ell)
                            locations += 1
    return certificates, locations


# --------------------------------------------------------------------------
# Exact tails of L_eta(u, j)
# --------------------------------------------------------------------------

def biased_interval_lengths(X, eta):
    """L_eta(j) for every position j, for one sign vector X."""
    M = len(X)
    prefix = [0]
    for x in X:
        prefix.append(prefix[-1] + x)
    L = [1] * M
    for lo in range(M):
        for hi in range(lo + 1, M + 1):
            if abs(prefix[hi] - prefix[lo]) >= eta * (hi - lo):
                ell = hi - lo
                for j in range(lo, hi):
                    if ell > L[j]:
                        L[j] = ell
    return L


def check_tails(rng):
    tails = 0
    biases = 0
    for m in range(1, 6):
        M = 1 << m
        orders = [list(range(M))]
        for _ in range(3):
            order = list(range(M))
            rng.shuffle(order)
            orders.append(order)
        for order in orders:
            for eta in ETAS:
                lengths = []
                for u in range(M):
                    X = [1 - 2 * popcount_parity(u & v) for v in order]
                    L = biased_interval_lengths(X, eta)
                    lengths.append(L)
                    biases += M
                for j in range(M):
                    column = [lengths[u][j] for u in range(M)]
                    for t in range(1, M + 1):
                        prob = Fraction(sum(1 for L in column if L >= t), M)
                        assert prob <= min(1, Fraction(1152) / (eta ** 4 * t)), (m, eta, t, prob)
                        tails += 1
                    expectation = sum(math.log2(L) for L in column) / M
                    assert expectation < 12 + 4 * math.log2(1 / eta)
    return tails, biases


# --------------------------------------------------------------------------
# Elementary inequalities
# --------------------------------------------------------------------------

def h2(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def check_inequalities():
    checks = 0
    for i in range(1, 1000):
        p = i / 1000
        assert 1 - h2(p) >= 2 / math.log(2) * (p - 0.5) ** 2 - 1e-12
        checks += 1
    for q in range(1, 41):
        for mu10 in range(1, 10 * q + 1):
            mu = mu10 / 10
            assert q * h2(mu / q) <= mu * math.log2(math.e * q / mu) + 1e-12
            checks += 1
    for m in range(2, 2001):
        for g in range(1, 100):
            gamma = g / 100
            if gamma <= 1 / m:
                continue
            lhs = (gamma - 1 / m) / (1 - 1 / m) * m
            assert lhs >= gamma * m - 2 - 1e-9
            checks += 1
    assert math.log2(1152) + math.log2(math.e) < 12
    # Tail integration constant: E log2 L <= log2 K + log2 e when
    # Pr[L >= t] <= min(1, K/t).
    for K in (1.0, 2.5, 1152.0, 1152.0 * 8 ** 4):
        integral = 0.0
        y, dy = 0.0, 1e-3
        while y < math.log2(K) + 60:
            integral += min(1.0, K / 2 ** y) * dy
            y += dy
        assert integral <= math.log2(K) + math.log2(math.e) + 1e-2
        checks += 1
    return checks


def main():
    rng = random.Random(SEED)
    orders, anchors, worst = check_ordered_runs(rng)
    print(f"ordered inner-product runs: {orders} column orders "
          f"(all 26 for m = 1, 2; 160 sampled for m = 3..7), "
          f"largest conditional run entropy {worst:.9f} bits, "
          f"{anchors} exact anchor probabilities")
    certificates, locations = check_dyadic_certificates()
    print(f"dyadic certificates: {certificates} interval/eta pairs over every sign "
          f"sequence of length 1, 2, 4, 8; {locations} anchor/location conditions")
    tails, biases = check_tails(rng)
    print(f"biased-interval tails: {tails} exact tail probabilities, "
          f"{biases} interval/bias maximizations, m <= 5")
    checks = check_inequalities()
    print(f"elementary inequalities: {checks} grid checks")
    print("interval_entropy_validation: all checks passed")


if __name__ == "__main__":
    main()

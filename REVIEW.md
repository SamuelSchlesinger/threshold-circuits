# Review and validation record

## Second-pass review, October 18, 2026

Reviewed the draft `main.tex` as received (1,778 lines, 107 KB) and prepared
the repository around it. The source now lives in
`threshold_circuit_representations.tex`; its mathematical content is
unchanged except for the four edits listed under "Changes made".

### Build

The draft compiled on the first attempt. After the edits below and the
priority-check additions, three `pdflatex` passes produce a 39-page PDF with
no warnings, no overfull or underfull boxes, all 54 labels resolving, and
every one of the 23 bibliography keys cited. `validation.txt` records the
final build and check output.

### Bibliography

Three recent references were checked against their sources:

- Håstad, ECCC TR25-209 (2025), "Efficiently finding small representations
  for LTFs": title, author, and year match.
- Chen, Tal, and Wang, ECCC TR26-039 (2026), "Super-quadratic lower bounds
  for depth-2 linear threshold circuits": title, authors, and the
  `n^{2.5 - eps}` bound for a function in E^NP match the draft's description.
- Dahiya, Vignesh K, Mahajan, and Sreenivasaiah, IPL 183 (2024): the
  manuscript indeed records an `exp(Omega(sqrt n))` size lower bound for
  linear decision trees computing inner product, credited there to the
  rank-depth tradeoff of Turán and Vatan. The draft's statement that its
  `2^{m/2}` bound is stronger than that recorded statement is accurate.

One suggestion, not applied: the sentence in Section 6.2 citing TR25-209 for
the existence of `O(v log v)`-bit integer weights states a classical fact
(Muroga; tight by Håstad 1994), whereas TR25-209 concerns finding such
representations efficiently. Adding the classical reference next to it would
make the attribution precise.

### Proof reading

Every proof was read line by line. No mathematical error was found. Points
checked in detail:

- Section 2. The quantile lemma with pivots `i = ceil(a u)`, `j = ceil(a v)`:
  the inclusive suffixes have sizes `u - i + 1 >= (1 - a) u` because
  `i <= a u + 1`. The leaf-count induction telescopes correctly, and the
  Hadamard rectangle bound follows from `|1_A^T H 1_B| <= 2^{m/2} sqrt(|A||B|)`.
- Section 3. The generating-function estimate `t^e B(q, e) <= (1 + t)^q <= e^e`
  at `t = e / q` gives (5). The constants `h2^{-1}(1/4) = 0.04169269027...`
  and `h2^{-1}(1/2) = 0.1100278644...` were recomputed to the printed digits.
  The depth-three block construction was evaluated exhaustively for `m <= 5`.
- Section 4. The exact capacity "iff": transcript fibers are intervals because
  fixing the transcript makes every remaining condition a half-line in the
  score; the converse builds an ordered tree with any leaf count in
  `[1, B(q, e)]` by Pascal's identity. The enumerative-coding invariant
  `a_i <= t < a_i + B(q - i + 1, e - r_i)` is preserved because the zero and
  one subblocks have sizes `B(q - i, e - r_i)` and `B(q - i, e - r_i - 1)`.
  The coefficient bound uses `q >= N`, which follows from `B(q, e) >= 2^N`.
  Corollary 4.4's two binomial inequalities hold for both parities of N.
- Section 5. In Lemma 5.1, `2^k + 1` distinct columns span an F2-space of
  dimension at least `k + 1`, giving `Pr[A > 2^k] <= 2^{-k-1}`; summing over
  `k < m` gives `E log2 A < 1`. In Lemma 5.6 the constants compose as
  `18 * 4 * 16 = 1152`, and `log2 1152 + log2 e = 11.61 < 12`. The count of
  dyadic blocks of size h meeting an interval of `2 l + 1` positions is at
  most `2 l / h + 2`, as used. The passage from (E8) to (E1) with
  `eta = 1 / m` is `gamma m - (1 - gamma) m / (m - 1) >= gamma m - 2`. The
  optimality example decodes Bob's bits with m collinear queries.
- Section 6. The Fibonacci identity (P1) was verified symbolically and for
  `s <= 100`. The width-two 3-CNF path encoding and its totally unimodular
  internal matrix are correct.
- Section 7. The capacity proof reproduces the statement formalized as
  `TwoSidedRectangleFree.two_pow_lt_mul_cost` in `complexitylib`. The menu
  theorem is a valid block decomposition: each precomputed table is one
  scalar block, after which gates of that direction are input-free updates.
- Section 8. The cutting lemma relies on Chazelle's theorem with an
  unspecified constant `C_r`; the perturbation argument for degeneracies is
  sound and the draft claims nothing about the constant.
- Section 9. The lattice-point formula (Z2) is Stanley's zonotope formula in
  the arithmetic-Tutte form of D'Adderio and Moci; the index bound
  `m(J) <= |det M|` and the Hadamard step are correct, and the final
  `(e k W / sqrt r)^r` estimate uses `t = r / (k a) <= 1`, which holds since
  `sqrt r <= k W`. The total-unimodularity integrality argument via vertices
  of `{R x = v, 0 <= x <= 1}` is standard.
- Sections 10 and 11. The big-M inequalities (25) force `Q_i` into the
  correct sign range for each `z_i`, and the priorities `lambda_i` make each
  gate's own margin dominate all later terms, giving a unique sink on every
  face.
- Section 12. The `P/poly = L/poly` barrier and the width-two linear-query
  program for inner product are correct as stated.

### Changes made

1. The appendix cited four validation scripts "included in the source
   archive" that were not present. The four scripts now exist, run with
   `make check`, and implement the checks the appendix describes. The
   scalar-coding counts (66 pairs, 11,264 ranks, 103,423 comparisons,
   114,687 invariants) reproduced the appendix exactly; the other three
   scripts produce different totals (for example 283,199 rectangle entries
   instead of 289,819, and a largest conditional run entropy of 1.306174414
   instead of 1.304510338 bits), so the appendix now reports the numbers
   these scripts actually produce. Several further checks were added and are
   listed in the appendix.
2. Two leftover integration comments ("Standalone insert body ...
   second_pass_menus.md" and "Suggested placement: Section 6 ...") were
   removed.
3. After inequality (18), "In that case" now makes explicit that the
   `s >= 2^{(1/2 - o(1)) N}` conclusion is the one-direction case.
4. In the proof of Lemma 5.6, the sentence bounding the discarded dyadic
   blocks now says "twice the sum of the dyadic sizes below a cutoff h is
   less than 4h", which is what the estimate uses.

### Open items for the author

- The `\author` and `pdfauthor` fields read "Research draft prepared for
  Samuel Schlesinger" and should be replaced before circulation.
- The draft's own caveat stands: historical priority of the formulations is
  not established, and the second-pass results were not checked by a proof
  assistant.
- Coverage gaps in the finite checks: the cutting lemma (Section 8), the
  entropy identities of Theorems 5.2 and 5.5 beyond the run and tail
  estimates, and the running-time claims of Section 9 are not exercised.

## Priority check, October 18, 2026

The question was whether the draft's headline inequality, the depth-free and
weight-free size–energy tradeoff of Theorem 3.2, and the `2^{m/2}` leaf bound
for linear decision trees behind it, already exist. Sources read in full text
or abstract: Uchizawa, Douglas, Maass (Neural Computation 2006), Uchizawa's
survey (Interdisciplinary Information Sciences 2014), Uchizawa and Abe (MFCS
2023, arXiv 2107.00223), Uchizawa and Takimoto (TCS 2008; MFCS 2011; SOFSEM
2015, abstracts), Gröger and Turán (ICALP 1991, abstract), Turán and Vatan
(1997, as cited), Chattopadhyay, Mahajan, Mande, Saurabh (CJTCS 2020), and
Dahiya, Vignesh K, Mahajan, Sreenivasaiah (IPL 2024), plus three web searches
for unrestricted leaf-count bounds.

Findings:

- The pattern-to-tree conversion is Lemma 1 of Uchizawa, Douglas, and Maass:
  a circuit with `2^H` attainable gate-state patterns has a threshold decision
  tree with `2^H` leaves and depth at most the gate count. They use it only
  for upper bounds (energy reduction, depth-three collapse).
- Prior energy lower bounds for inner product all carry a depth factor:
  `s >= 2^{Omega(m / e^d)}` (Uchizawa–Takimoto 2008, as restated by
  Uchizawa–Abe) and `log rk(M_C) <= e d (log s + log w + log m)`
  (Uchizawa–Abe 2023, Theorem 1), the latter also needing integer weights at
  most w. Uchizawa–Abe's introduction says these bounds "break down when we
  consider threshold circuits of larger depth and energy".
- The 2014 survey's open-problems section asks, verbatim, whether the depth
  term "can also be dropped off, which implies a lower bound for threshold
  circuits of bounded energy, but with unbounded weight and depth". Theorem
  3.2 is that statement.
- Linear decision tree size for inner product with unrestricted weights: the
  best recorded bound is `exp(Omega(sqrt m))`, via the Turán–Vatan rank–depth
  tradeoff (Dahiya et al. 2024 record it in this form). Uchizawa–Takimoto
  obtain `2^{Omega(m / log w)}` only for total absolute weight at most w
  (2015, improving their 2011 energy-based argument), which indicates the
  unrestricted exponential bound was not known to them. No later unrestricted
  bound was found.
- The quantile rectangle step is in CMMS20 (Lemma 16) for decision lists;
  the leaf-weighted induction that extends it to trees was not found.

Conclusion: no earlier depth-independent energy bound and no earlier
unrestricted `2^{Omega(m)}` leaf bound were located. The manuscript now states
the comparison precisely in Section 1.2, adds a comparison paragraph and the
general rectangle-density corollary after Theorem 3.2, and cites UT08, UT11,
UT15, GT91, Muroga 1971, and Håstad 1994. The paywalled 2011 paper's exact
exponent was not retrievable; it is cited for its method and weight
restriction only. A search is not a proof of priority, and the ledger says so.

For the other components, two quick searches found nothing on enumerative
coding inside threshold circuits or on zonotope-based counting for threshold
circuits; those searches were shallow and the draft's own caveat stands.

## Check-script coverage

| Script | Statements exercised | Method |
|---|---|---|
| `weighted_pivot_validation.py` | Lemma 2.1, Theorem 2.2, Lemma 2.4, Lemma 3.1 with Theorem 3.2, Section 3.3, Corollary 3.4 constants | Random trees with tied scores; exhaustive rectangles for m <= 3; random circuits; exhaustive block circuits for m <= 5 |
| `interval_entropy_validation.py` | Lemma 5.1, Lemma 5.6, the inequalities in Corollary 5.3 and Theorem 5.5 | All orders for m <= 2 and sampled orders to m = 7; every sign sequence of length <= 8; exact tails for m <= 5 |
| `zonotope_validation.py` | Lemma 9.2, Proposition 9.5, Corollary 6.2's identity, Theorem 9.3, Theorem 6.1, Proposition 6.3, Theorem 11.1 | Exact rational arithmetic; gcd-of-minors lattice indices; dictionary counting against enumeration; all faces of the gate cube |
| `scalar_coding_validation.py` | Theorem 4.1 (construction), Theorem 4.3 (coding and gate level), Corollary 4.4 | Exhaustive over q <= 10 ranks and all inputs for N <= 9 |

## Lean 4 formalization prospects

The question is whether the new results can be checked in Lean 4 in the style
of the author's existing complexity developments. The honest answer is that a
substantial core can, building directly on what already exists, while the
entropy, cutting, and zonotope results would need infrastructure that neither
Mathlib nor CSLib has today.

What exists already:

- `complexitylib` (`Complexitylib/Algebraic/LowerBound/Threshold`) defines
  `Program n s`, the dense real-weight threshold straight-line program this
  paper uses, with `inputGates`, `evalFrom` for runs along one real parameter,
  `TwoSidedRectangleFree`, `BlockDecomposition`, and the capacity theorem
  `two_pow_lt_mul_cost` (Theorem 7.2 here), together with
  `sourceReductionHardFamily` and its asymptotic corollaries (the explicit
  family of Section 7). The capacity proof reproduced in Section 7 is
  therefore already formal.
- `complexitylib/Metacomplexity/Hamming.lean` has Boolean Hamming balls with
  their exact binomial volume, which is the `B(q, e)` of this paper.
- Mathlib has `Real.binEntropy` and its inverse monotonicity facts,
  `Matrix.IsTotallyUnimodular`, Hadamard matrices, and `hammingDist`.
- CSLib's `Circuit σ` is generic over a `Signature` and `Interpretation`; a
  threshold signature (`Op := Σ k, (Fin k → ℝ) × ℝ` with arity k) would give
  threshold circuits with CSLib's composition and size bookkeeping.

A tiered plan, from most to least tractable:

1. **Direct extensions of the capacity API.** Theorem 7.3 (direction menus)
   is a `BlockDecomposition` whose scalar block for direction a is the vector
   of all comparisons at all menu offsets; `ChangesAtMost` for that vector
   is a sum of single-change coordinates, and the corollaries are arithmetic.
   Corollary 4.2 and the lower direction of Theorem 4.1 need one new lemma:
   the fiber of a restricted transcript is an interval of scores, so the
   number of monochromatic runs is at most the number of attainable
   transcripts, which `Hamming.lean` bounds by `B(q, e)`. Proposition 5.4 is
   the same counting on the inner-product matrix.
2. **The leaf-count theorem without trees.** Theorem 2.2 is best stated for
   programs directly: by induction on gates, every rectangle `A × B` has a
   sub-rectangle with sides at least `|A| / P` and `|B| / P` on which the
   input-reading transcript is constant, where P counts transcripts
   attainable on `A × B`. This avoids defining comparison trees. The sorting
   machinery in `Internal/Cut.lean` (`exists_median_internal`,
   `exists_sorted_group`) already handles the quantile step. Lemma 2.4 is a
   finite character-sum computation over `ZMod 2`; `Circuits/Correlation.lean`
   and `BooleanAnalysis.lean` are the places to look for parity
   orthogonality. Theorem 3.2 and Corollary 3.3 then follow from
   `Hamming.lean` and `binEntropy`.
3. **Explicit constructions with gate counts.** Section 3.3's block circuit,
   Theorem 4.3's enumerative-coding normalization, Theorem 6.1's pending-bias
   simulation, Proposition 6.3, and Theorem 11.1's unique-sink potential are
   finite constructions of the kind the iterated-functions Lean development
   already does. Theorem 4.3 is the longest: it needs the `B` recursion, the
   flag gates, and a decoder, with a size bound `S + O(q e + N)`. The Python
   check mirrors the gate list exactly and can serve as the specification.
4. **Blocked or expensive.** Theorems 5.2 and 5.5 use Shannon entropy, which
   Mathlib does not provide (the PFR project has its own); they can be
   restated as bounds on the expected logarithm of fiber sizes, which is how
   the proofs run, but the step to `q h2(mu / q)` needs subadditivity. Lemma
   5.6 is finite but intricate. Lemma 8.2 depends on Chazelle's cutting
   theorem, and Lemma 9.2 on the lattice-zonotope formula; neither is
   realistic to formalize soon, so Sections 8 and 9.1 would remain prose
   except for the dictionary algorithm's correctness (easy) and the
   totally unimodular case, where the independent-set count might be proved
   combinatorially using `IsTotallyUnimodular`.

Model recommendation: keep the lower bounds in `complexitylib`'s
`Algebraic.Threshold.Program`, whose gate count is the index `s` and whose
input-reading set is already defined, and use CSLib's `Circuit` only where
composition matters (the constructions in tier 3), connected by a simulation
lemma that evaluates a CSLib threshold circuit as a `Program`. The uniform
polynomial-time hypothesis in Corollary 4.5 is the one place where CSLib's
Boolean circuit model enters, to import polynomial-size Boolean circuits for
the explicit family.

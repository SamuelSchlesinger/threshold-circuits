# Threshold Circuits as Decision Trees and Finite-State Programs

Research draft, second pass, October 18, 2026.

Which resource of a threshold circuit does a target function force to be
large? The paper reads an arbitrary-depth, arbitrary-real-weight threshold
circuit through four exact representations and attaches a quantitative
consequence to each: the activation transcript of its input-reading gates,
adaptive comparisons along a scalar input direction, the vector of pending
biases across a topological cut, and the additive summary of its input rows.

Main statements, with their hypotheses kept separate:

- **Activity at arbitrary depth.** A leaf-count rectangle theorem for
  separable-comparison trees gives, for inner product on m + m bits,
  `sum_{j <= e} C(q, j) >= 2^{m/2}` where q counts input-reading gates and
  e is their worst-case activity. Linear size forces linear activity; the
  known polynomial-size `O(m / log m)`-energy construction is included.
- **Exact scalar capacity.** For one input direction, q scalar queries of
  activity at most e compute f exactly iff the number of monochromatic score
  runs is at most the Hamming-ball size B(q, e). An enumerative-coding
  construction realizes any size-S circuit with N + 1 collinear input-reading
  gates, activity at most `ceil(N/2)`, and size `S + O(N^2)`.
- **Average activity under one-sided rank one.** If one party's direct input
  vectors are collinear, the transcript has conditional entropy above m - 2;
  with error at most eps it is at least `(1 - 2 eps) m - 4 log m - 14`, and
  the leading coefficient is optimal.
- **Directions, menus, and input rank.** Static expansion of threshold offsets
  and Chazelle's cutting theorem give exponential input-reading gate lower
  bounds for bounded predecessor fan-in with few directions or bounded input
  rank, on a prior explicit all-splits rectangle-free family.
- **Exact counting from additive summaries.** Integer zonotopes bound the
  number of summaries by `(e N W / sqrt r)^r`, giving an exact #SAT algorithm
  for rational circuits whose input rows lie in a supplied integer rank-r span;
  totally unimodular summaries are counted by matroid independent sets, and
  totally unimodular internal cuts give Fibonacci-size linear-query programs.

Polyhedral, integer-programming, and quadratic unique-sink interpretations
complete the map. No unrestricted superlinear total-gate lower bound is
claimed; the final section states what each argument does not prove.

## Read and use

- [Paper (PDF)](threshold_circuit_representations.pdf) and
  [LaTeX source](threshold_circuit_representations.tex). The bibliography is
  embedded; no BibTeX run or external figure is needed.
- [Review and validation record](REVIEW.md) and
  [finite-check output](validation.txt).
- Finite check scripts: [weighted_pivot_validation.py](weighted_pivot_validation.py)
  (leaf-count rectangles, activation patterns, the block construction),
  [interval_entropy_validation.py](interval_entropy_validation.py) (ordered
  runs, dyadic certificates, biased-interval tails),
  [zonotope_validation.py](zonotope_validation.py) (additive images, exact
  counting, pending biases, the quadratic potential), and
  [scalar_coding_validation.py](scalar_coding_validation.py) (comparison
  trees, enumerative coding, the gate-level normalization).
- [Submission metadata](submission_metadata.txt).

## Build and check

```sh
make         # Compile the paper; requires TeX Live or MiKTeX.
make check   # Run the four finite semantic checks; requires Python 3 only.
```

The checks compare constructions with exhaustive enumeration and verify the
exact inequalities used in the proofs on small instances. They do not prove
any statement for unbounded parameters.

## Status and attribution

The all-splits capacity theorem and its explicit family are prior results
from the author's `complexitylib` Lean development and are reproduced with
attribution; the enumerative coding scheme, hyperplane cutting theorem, and
lattice-zonotope formula are classical ingredients cited in the text. The
remaining statements are proved in the manuscript. Historical priority of
the formulations and applications remains unresolved, the author field is a
placeholder pending confirmation, and no proof-assistant verification of the
new results is claimed. The review record lists the open items and a plan
for formalizing parts of the paper in Lean 4.

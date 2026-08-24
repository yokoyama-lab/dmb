# dmb — computing discrete Morse–Bott theory

[![CI](https://github.com/yokoyama-lab/dmb/actions/workflows/ci.yml/badge.svg)](https://github.com/yokoyama-lab/dmb/actions/workflows/ci.yml)

Computes and visualises **discrete Morse–Bott functions** on finite complexes
(simplicial complexes and CW complexes that need not be regular), following

> Y. Nishikawa and T. Yokoyama, *On discrete Morse–Bott theory*,
> [arXiv:2511.07864](https://arxiv.org/abs/2511.07864) v2.

It is a rewrite of `dmf.py`, which drew discrete Morse functions on a torus;
discrete Morse theory is recovered as a special case (Theorem 3.2 of the paper).

| collections of the torus height function | the same on the torus |
|---|---|
| ![collections](docs/img/collections-2d.png) | ![torus](docs/img/collections-3d.png) |

| a discrete Morse function with 4 critical cells | a rotation-invariant one (4·n critical cells) |
|---|---|
| ![dmf](docs/img/dmf-vectorfield-2d.png) | ![invariant](docs/img/invariant-dmf-2d.png) |

Arrows are strictly noncritical pairs (Forman's V-paths when the function is a
discrete Morse function); red circles mark weakly critical cells, i.e. the cells
of the reduced collections.

## Install

The theory code has **no dependencies** and runs on Python 3.9+. Only the
visualiser (`dmb.py`) needs dash, numpy and plotly.

```bash
git clone https://github.com/yokoyama-lab/dmb
cd dmb
pip install -r requirements.txt      # only for the visualiser
                                     # (or: pip install -e ".[app]")
```

## Usage

```bash
python3 examples.py                  # worked examples, then what DMBT can do that DMT cannot
python3 examples.py tutorial         # small examples you can check by hand
python3 examples.py strength         # only the things discrete Morse theory cannot do

python3 dmb_core.py                  # report on 7 functions on T(4,4)
python3 dmb_core.py 5 4              # on the torus with ni=5, nj=4
python3 dmb_core.py --table          # how symmetry degrades DMT but not DMBT
python3 dmb_core.py --field 2        # over F_2 (torsion becomes visible)
python3 dmb_core.py --json           # machine-readable output
python3 dmb_core.py --complex my.json --json     # your own complex

python3 complexes.py                 # catalogue of test complexes with Betti numbers
python3 search.py --ni 4 --nj 3 --values 2       # exhaustive search of invariant DMBFs
python3 export.py --format tikz -o fig.tex       # TikZ figure for the paper
python3 dmb.py                       # visualiser at http://127.0.0.1:8050/
python3 dmf.py                       # the original discrete Morse figure (unchanged)
```

As a library:

```python
import dmb_core as D

K = D.torus(5, 4)                    # triangulated torus
f = D.height_fn(K, 5, 4)             # rotation-invariant discrete Morse-Bott function
M = D.MorseBott(K, f)
M.is_dmb()                           # True: satisfies (M1)(MB2)(M3)(MB4)
M.collections()                      # equivalence classes of the r-path relation
M.reduced_collections()              # their weakly critical parts
M.report()["R_MB"]                   # [] i.e. R(t) = 0 (equality)

D.homology_z(K.cells)                # integral homology, including torsion
D.MorseBott(K, f, p=2).report()      # over F_2
```

Bring your own complex with `D.Complex([...])` (pass the facets), or
`D.CWComplex(...)` for a CW complex that need not be regular; JSON input is
supported through `complexes.load_json`.

```python
K = D.Complex([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)])   # S^2 = boundary of a 3-simplex
K = D.cw_projective_plane_minimal()                            # RP^2 with 3 cells
```

## Layout

```
dmb_core.py        the theory (simplicial and CW complexes, DMB/DMF tests, collections,
                   homology, Theorem 4.12, Morsification, drawing coordinates)  no deps
complexes.py       test complexes (S^1, S^2, disc, Mobius, RP^2, Klein bottle, ...)
                   and JSON input/output                                        no deps
examples.py        worked examples and the strengths of DMBT                     no deps
search.py          exhaustive search of Z_ni-invariant discrete Morse-Bott functions
export.py          figure export (TikZ needs nothing; HTML/PNG/SVG need plotly)
dmb.py             Dash visualiser                                      dash/numpy/plotly
dmf.py             the original discrete Morse figure (unchanged)       dash/numpy/plotly
tests/             165 tests (visualiser tests skip when dash is absent)
docs/results.md    results and proofs (layer A: machine-checked / layer B: on paper)
docs/img/          figures used in the README
```

## Definitions implemented

For a k-cell σ,

- U(σ)     = #{τ : σ ≺ τ, dim τ = k+1, f(σ) ≥ f(τ)}
- D(σ)     = #{ν : ν ≺ σ, dim ν = k−1, f(ν) ≥ f(σ)}
- U^snc(σ) = #{τ : σ ≺ τ, dim τ = k+1, f(σ) > f(τ)}
- D^snc(σ) = #{ν : ν ≺ σ, dim ν = k−1, f(ν) > f(σ)}

- **discrete Morse** (Definition 12): (M1) (M2) U(σ) ≤ 1 (M3) (M4) D(σ) ≤ 1
- **discrete Morse–Bott** (Definition 18): (M1) (MB2) U^snc(σ) ≤ 1 (M3) (MB4) D^snc(σ) ≤ 1

The only difference between (M2) and (MB2) is whether neighbours of equal value
are counted; that is where collections come from. (M1)/(M3) require the value to
increase strictly along **irregular faces of any codimension**. On a simplicial
complex every face is regular, so they hold automatically; on a non-regular CW
complex they bite, and reading them in codimension one only gives the
counterexample to Lemma 3.1 of v1 of the paper (see example 6 of
`python3 examples.py tutorial`).

- **collection**: an equivalence class of the r-path relation (cells of equal
  value joined through codimension-one incidences) — the discrete analogue of a
  critical submanifold
- **weakly critical**: U^snc(σ) = D^snc(σ) = 0
- **reduced collection** C: all weakly critical cells of a collection
- **Theorem 4.12**: Σ_C P_t(C) = P_t(K) + (1+t) R(t) with R(t) ≥ 0, where P_t(C)
  is the Poincaré polynomial of (C_*(C), ∂^C), the boundary restricted to C

## Results

Four things discrete Morse theory cannot do (`python3 examples.py strength`;
proofs in [`docs/results.md`](docs/results.md)).

**1. On spaces with torsion the discrete Morse inequality is never an equality.**
A closed surface has H₂(K; Z/2) ≠ 0, so every discrete Morse function has m₂ ≥ 1,
while RP² and the Klein bottle have b₂ = 0 over the rationals. Hence R(t) ≠ 0
always — even for the minimal CW structure of RP², which has three cells.
A discrete Morse–Bott function reaches R(t) = 0. **The gap depends on the
coefficient field**: over F₂, where the torsion is visible, the same discrete
Morse function becomes sharp. Discrete Morse–Bott functions are sharp over every
field.

| complex | field | P_t(K) | M(t) of a minimal DMF | R_DMF | R_DMBF |
|---|---|---|---|---|---|
| RP² | Q | 1 | 1+t+t² | t | 0 |
| RP² | F₂ | 1+t+t² | 1+t+t² | 0 | 0 |
| Klein bottle | Q | 1+t | 1+2t+t² | t | 0 |
| Klein bottle | F₂ | 1+2t+t² | 1+2t+t² | 0 | 0 |

**2. Imposing symmetry degrades discrete Morse theory but not discrete Morse–Bott
theory.** For the rotation Z_ni acting freely on the cells of T(ni, nj):

| function | invariant | critical cells / circles | Σ_C P_t(C) | R(t) |
|---|---|---|---|---|
| DMF with 4 critical cells (tree–cotree) | no | 4 cells | 1+2t+t² | 0 |
| minimal Z_ni-invariant DMF | yes | 4·ni cells | ni(1+2t+t²) | (ni−1)(1+t) |
| Z_ni-invariant discrete Morse–Bott height | yes | 2 critical circles | 1+2t+t² | **0** |

Exhaustive search (`search.py`) shows that on T(3,3), T(4,3), T(3,4) and T(4,4)
with two values there are **only two kinds of sharp invariant discrete Morse–Bott
function**: the one with two critical circles contributing (1+t) and t(1+t) —
exactly the smooth Morse–Bott picture — and the trivial one whose single
collection carries all of P_t(T²). On T(4,4), for instance, 4410 of the 121842
invariant discrete Morse–Bott functions are sharp: 4408 of the first kind and 2
of the second.

**4. Torsion and symmetry obstruct independently (the Klein bottle).**
On the rotation-invariant triangulation of the Klein bottle (`klein_bottle_sym`)
the translation g has order 2·ni, and **only for odd ni** does g² generate a free
Z_ni action. An invariant discrete Morse function then suffers both penalties:

| field | minimal DMF (not invariant) | Z_ni-invariant DMF (lower bound) | invariant DMBF |
|---|---|---|---|
| Q | t | ≥ (ni−1)(1+t) + t | 0 |
| F₂ | 0 | ≥ (ni−1)(1+t) | 0 |

The difference is exactly t, the torsion part, which disappears over F₂.
Discrete Morse–Bott functions reach R = 0 over both.

**3. It reproduces smooth Morse–Bott functions.** The distance from the axis on a
rotationally symmetric torus has critical circles of index 0 and 1, giving
Σ_i t^λi P_t(S¹) = (1+t) + t(1+t) = 1 + 2t + t² = P_t(T²); the discrete version
produces two reduced collections with P_t = 1+t and t+t².

## Testing

```bash
python3 -m unittest discover -s tests -t . -v          # 165 tests
DMB_TRIALS=60 python3 -m unittest tests.test_properties   # more random trials
DMB_SLOW=1    python3 -m unittest tests.test_search       # compare against brute force
DMB_LATEX=1   python3 -m unittest tests.test_export       # actually compile the TikZ
DMB_RENDER=1  python3 -m unittest tests.test_dmb_app      # actually render the figures
ruff check .
```

Detection power is measured by mutation: lowering the value of a single edge
breaks (MB4) in all 48 cases; dropping the signs of the incidence numbers changes
the Betti numbers (so orientations are really used). The fast sparse elimination
is checked against a naive dense reference implementation.

## Limitations

- Theorem 4.12 is verified **over a field** (Q or F_p). Integral homology with
  torsion is available separately (`homology_z`, via Smith normal form).
- **Incidence numbers of a CW complex must be supplied** by the user
  (`CWComplex(..., incidence=...)`); `check_boundary()` verifies ∂∘∂ = 0.
- The visualiser and the torus-specific tools (`height_fn`, `invariant_dmf`,
  `export.py`, `search.py`) are for the torus. The theory (`MorseBott`) works on
  any finite complex.
- The visualiser allows 3–16 subdivisions (`MAX_GRID`); a report on T(16,16)
  (1536 cells) takes about half a second.

## Citation

```bibtex
@misc{nishikawa2025discrete,
  title  = {On discrete Morse--Bott theory},
  author = {Nishikawa, Y. and Yokoyama, T.},
  year   = {2025},
  eprint = {2511.07864},
  archivePrefix = {arXiv},
  primaryClass  = {math.GT},
  url    = {https://arxiv.org/abs/2511.07864}
}
```

## License

MIT License (see [LICENSE](LICENSE)).

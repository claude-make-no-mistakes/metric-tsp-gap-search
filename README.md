# Metric TSP four-thirds counterexample search

This experiment certified the exact worst metric ratio of each of 24 labeled non-half-integral subtour vertices on 17 cities. The largest was \(15/13\). The only sampled vertex with 26 support edges had ratio \(10/9\). The half-integral positive control reproduced the published \(53/43\). The search found no counterexample.

[Research note](note.pdf) · [LaTeX source](note.tex)

| Search | Seed | Attempts | Certified candidates | Largest ratio |
| --- | ---: | ---: | ---: | ---: |
| [Unrestricted chords](global-ledger.jsonl) | 61125 | 1,600 | 20 | \(23/20\) |
| [Restricted chord endpoints](local-ledger.jsonl) | 61717 | 1,800 | 2 | \(15/13\) |
| [Objectives biased toward the control metric](guided-ledger.jsonl) | 61917 | 2,200 | 2 | \(129/113\) |

There were 15 supports with 24 edges, eight with 25, and one with 26. Counts exclude the repeated positive control. Distinctness uses the labeled rational vector. Candidate generation is a sample of objectives and supports.

## Reproduce

The verifier uses the Python 3.10+ standard library. It checks degrees and cuts with exact flows, confirms extremality by rational row reduction, checks each Eulerian multigraph and its mixture weight, and recomputes the minimum tour with integer subset dynamic programming.

    python verify.py global-ledger.jsonl
    python verify.py local-ledger.jsonl
    python verify.py guided-ledger.jsonl

To regenerate the ledgers, use Python 3.12, NumPy 2.5.3, SciPy 1.18.1, and a GCC-compatible 64-bit C++ compiler. On Windows with MinGW-w64, build the library from the repository directory.

    g++ -std=c++17 -O3 -Wall -Wextra -shared -static exact.cpp -o exact.dll

On Linux, use the shared-object build.

    g++ -std=c++17 -O3 -Wall -Wextra -shared -fPIC exact.cpp -o exact.so

The output path must be new. LP tie-breaking can vary by build. Verify the generated certificates with the same standard-library checker.

    python search.py replay-global.jsonl --attempts 1600 --cases 30 --rounds 320
    python search.py replay-local.jsonl --attempts 1800 --cases 30 --rounds 320 --seed 61717 --local
    python search.py replay-guided.jsonl --attempts 2200 --cases 30 --rounds 320 --seed 61917 --guided

## Scope and prior work

Cook, Hougardy, and Petrich give exhaustive results through 16 cities in general and 18 on half-integral vertices. Villa, Vercesi, Barta, and Mastrolilli give a theorem for support size at most \(n+6\). Further public work claims extensions to \(n+7\) and \(n+8\). These claims overlap with most of this pilot sample. The 26-edge case falls outside those support thresholds. [Exact small-order computation](https://arxiv.org/html/2603.12995v5) · [Published sparse-support theorem](https://arxiv.org/abs/2507.07003v3) · [\(n+7\) artifact](https://github.com/GiladR1979/metric-tsp-nplus7-certificates) · [\(n+8\) preprint](https://ssrn.com/abstract=7395957) · [\(n+8\) certificates](https://zenodo.org/records/21870868)

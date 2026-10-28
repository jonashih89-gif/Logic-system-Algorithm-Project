# Quine-McCluskey Solver (Logic System Algorithm, HW2)

A C++ implementation of the **Quine-McCluskey (tabular) method** for Boolean function minimization. Given the minterms and don't-care terms, it finds all **prime implicants** and **every minimum Sum-of-Products (SOP)** solution, using **Petrick's method** when the essential prime implicants alone are not enough.

Unlike a K-map, the tabular method is not limited to 4 variables; it has been tested on functions with 8 and 10 variables.

This repository contains only the core algorithm (`solveQm()`) and a result-printing helper. The program was originally submitted to a course auto-grading platform that supplied the I/O driver; that platform has since been closed.

## Interface

```cpp
QmSolution* solveQm(int numVar, vector<int> minterms,
                    vector<int> dontcares, bool verbose);
```

The result `QmSolution` contains the number of variables, the prime implicants in literal form (e.g. `a'bc'd`), and a list of all minimum SOP solutions. Variables are named `a, b, c, ...` from the most significant bit. When `verbose` is `true`, the program also prints the column-by-column merging table and the prime implicant chart.

## Algorithm

**1. Initialization.** Every minterm and don't care becomes an implicant `Imp` holding a binary string (`parseString`, MSB first, zero-padded) and a `set<int>` of the original indices it covers. Carrying the covered set along means merged implicants never need to be traced back to their sources.

**2. Column-by-column merging.** Each column, all pairs of implicants are compared. Two implicants merge if they differ in exactly one bit and their `-` positions match. The differing bit is replaced with `-`, and the two covered sets are united. A `map<string, set<int>>` collects the next column, so identical implicants produced by different pairs are automatically deduplicated.

**3. Collecting prime implicants.** Any implicant that could not be merged in its column is prime. Only the actual minterms in its covered set are kept, so implicants made purely of don't cares are discarded. The loop stops when no new merges occur.

**4. Prime implicant chart.** A P × M boolean matrix records which prime implicant covers which minterm.

**5. Essential prime implicants.** The chart is scanned repeatedly: if an uncovered minterm is covered by only one remaining prime implicant, that implicant is essential, and all minterms it covers are removed. This repeats until no new essential implicant is found. If every minterm is now covered, the essentials form the unique minimum SOP.

**6. Petrick's method.** For the remaining minterms, each one contributes a sum (OR) of the primes that cover it. These sums are multiplied out step by step as a Cartesian product; after each step every product term is sorted and deduplicated (absorbing `X·X = X`), and duplicate terms are removed. The terms with the fewest primes are the minimum covers, and each one combined with the essentials is returned as a separate minimum SOP.

## Verbose Output

In verbose mode, each column of the merging table is grouped by the number of 1s, and each entry is marked:

| Mark | Meaning |
|---|---|
| `v` | merged into the next column |
| `d` | a single don't-care term |
| `x` | duplicate of another entry in the same column |
| (blank) | not merged, i.e. a prime implicant |

The prime implicant chart is then printed with an `x` wherever an implicant covers a minterm.

## Techniques Used

Bit manipulation for binary conversion, string-based representation of implicants with `-` for eliminated variables, set union to track coverage, `std::map` for automatic deduplication of merged terms, an iterative essential-implicant search on a 2D boolean chart, and Petrick's method implemented as a Cartesian product with `sort` and `unique` for simplification. The code uses STL containers (`vector`, `set`, `map`), lambdas, and `iomanip` for formatted table output.

## Example

```
numVar    = 10
minterms  = 0, 1, 16, 17, 128, 343, 512, 640, 1023
dontcares = 341

primes = ["abcdefghij", "a'bc'de'fg'hj", "b'd'e'f'g'h'i'j'", "a'b'c'd'e'g'h'i'"]

minimumSops = [
  ["abcdefghij", "a'bc'de'fg'hj", "b'd'e'f'g'h'i'j'", "a'b'c'd'e'g'h'i'"]
]
```

Here the don't care 341 lets minterm 343 merge into the larger implicant `a'bc'de'fg'hj`.

## Limitations

Merging compares every pair of implicants in a column (O(n²) per column) rather than only adjacent 1-count groups. Petrick's method returns exact minimum solutions, but its product expansion grows exponentially in the worst case, so it is intended for moderate problem sizes. "Minimum" is measured by the number of product terms, not by the total number of literals.

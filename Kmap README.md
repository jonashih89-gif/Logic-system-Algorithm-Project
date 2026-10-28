# K-map Solver (Logic System Algorithm, HW1)

A C++ implementation that solves 2- to 4-variable Karnaugh maps. Given the minterms and don't-care terms of a Boolean function, it builds the K-map and finds the **prime implicants**, the **essential prime implicants**, and a **minimum Sum-of-Products (SOP)** expression.

This repository contains only the core algorithm (`solveKmap()`). The program was originally submitted to a course auto-grading platform that supplied the I/O driver; that platform has since been closed.

## Interface

```cpp
KmapSolution *solveKmap(int numVar, const vector<int> minterms,
                        const vector<int> dontcares);
```

The result is returned in the fixed structure `KmapSolution`, which holds the number of variables, the K-map (`'1'` = minterm, `'x'` = don't care, `'0'` = otherwise), and the three lists of product terms. Variables are named `a, b, c, d` from the most significant bit to the least significant bit.

## Algorithm

**1. Gray code layout.** Rows and columns are ordered by Gray code (`{0, 1, 3, 2}` for 2-bit axes), so physically adjacent cells differ in exactly one bit. The low bits (`c` or `cd`) index the rows and the high bits (`ab`) index the columns. Inverse lookup functions (`invertRowGray`, `invertColGray`) map a K-map position back to its minterm index.

**2. Filling the map.** A lambda `fillCell` splits each input value into its row and column bits and writes `'1'` or `'x'` into the corresponding Gray-coded cell.

**3. Exhaustive group enumeration.** Every cell is tried as a starting point, combined with every rectangle size whose height and width are powers of two. Indices are taken modulo the map size, so groups that wrap around the edges are handled naturally. A rectangle is a valid candidate only if all of its cells are `'1'` or `'x'` and at least one is `'1'`. Each group is stored as a sorted list of minterm indices, and a string key in a `std::set` removes duplicates.

**4. Prime implicant filtering.** A candidate is discarded if it is a proper subset of a larger candidate (checked with `std::includes` on the sorted lists). The remaining groups are the prime implicants.

**5. Expression conversion.** `bitsToExpression` checks each bit position across the whole group (including don't cares): a bit that is constant becomes a literal (`a` or `a'`), and a bit that varies is dropped. Only actual minterms are recorded in each implicant's coverage set, so don't cares help form larger groups but never need to be covered.

**6. Essential prime implicants.** A coverage table (`std::map<int, vector<string>>`) lists which implicants cover each minterm. Any minterm covered by exactly one implicant makes that implicant essential.

**7. Minimum SOP.** Starting from the essential implicants, the remaining uncovered minterms are covered greedily: at each step the implicant covering the most uncovered minterms is added, until all minterms are covered.

## Techniques Used

Gray code ordering and inverse lookup, bit manipulation for splitting indices and building literals, modular arithmetic for wrap-around adjacency, brute-force enumeration over power-of-two rectangles, subset elimination with `std::includes`, and a greedy set-cover heuristic. The implementation relies on STL containers (`vector`, `set`, `map`, `ostringstream`) and a C++11 lambda.

## Example

```
numVar = 4
minterms:  0 4 5 10 11
dontcares: 1 13 14 15

K-map:
1 1 0 0
x 1 x 0
0 0 x 1
0 0 x 1

Prime implicants:           a'c', bc'd, ac
Essential prime implicants: a'c', ac
Minimum SOP:                a'c' + ac
```

## Limitations

Only 2 to 4 variables are supported. The final cover step is greedy, so for maps with cyclic prime-implicant charts (where no single choice is forced) the result is a valid cover but is not guaranteed to be the absolute minimum.

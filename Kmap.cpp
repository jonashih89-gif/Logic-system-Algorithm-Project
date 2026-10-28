#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <algorithm>
using namespace std;

// ==============================
// The structure below must NOT be modified. It stores the final K-map solution,
// including the number of variables, the K-map itself, the prime implicants,
// the essential prime implicants, and the minimum SOP expression.
struct KmapSolution {
    int numVar;
    char **kmap;
    vector<string> primes;
    vector<string> essentials;
    vector<string> minimumSop;
};

// ==============================
// Convert a group of cells (given as decimal minterm indices) into a product term.
// If every cell in the group has the same value at a bit position, that variable is
// kept in the term (complemented, e.g. a', if the bit is 0).
// If the bit differs within the group, the variable is a "don't care" ("-") and is
// simply omitted from the term.
string bitsToExpression(int numVar, const vector<int>& group) {
    if (group.empty()) return "";
    string expr = "";
    int first = group[0];
    // Check bits from the most significant to the least significant
    for (int i = 0; i < numVar; i++) {
        bool same = true;
        int bit = (first >> (numVar - 1 - i)) & 1;
        // Check whether all cells agree on bit i
        for (size_t j = 0; j < group.size(); j++) {
            if (((group[j] >> (numVar - 1 - i)) & 1) != bit) {
                same = false;
                break;
            }
        }
        // If they agree, keep this variable in the term
        if (same) {
            char var = 'a' + i; // Variable order: a is the MSB, then b, ...
            expr += (bit ? string(1, var) : string(1, var) + "'");
        }
    }
    return expr;
}

// ==============================
// Gray code tables (used to order K-map rows and columns so that adjacent cells
// differ in exactly one bit)
// 2 variables: simply {0,1}
// 3 variables: rows use {0,1} (variable c), columns use {0,1,3,2} (Gray code of ab)
// 4 variables: rows and columns both use {0,1,3,2}, i.e. low bits (cd) and high bits (ab)
static int rowGray2[2] = {0, 1};
static int colGray2[2] = {0, 1};

static int rowGray3[2] = {0, 1};         // 3 variables: rows represent c
static int colGray3[4] = {0, 1, 3, 2};   // 3 variables: columns represent ab (Gray code)

static int rowGray4[4] = {0, 1, 3, 2};   // 4 variables: rows represent cd (Gray code)
static int colGray4[4] = {0, 1, 3, 2};   // 4 variables: columns represent ab (Gray code)

// ------------------------------
// The functions below perform the inverse Gray code lookup: they convert a K-map
// position (r, c) back to the original binary value, which is needed to compute
// minterm indices.
// e.g. invertRowGray() searches the row Gray code table for the entry matching the
// given row position and returns its original-order value.
int invertRowGray(int nRow, int r) {
    if (nRow == 2) {
        // For 2 rows the mapping is the identity
        for (int i=0; i<2; i++)
            if (rowGray2[i] == r) return i;
    }
    else {
        for (int i=0; i<4; i++)
            if (rowGray4[i] == r) return i;
    }
    return 0; // Default return value
}
int invertColGray(int nCol, int c) {
    if (nCol == 2) {
        for (int j=0; j<2; j++)
            if (colGray2[j] == c) return j;
    }
    else if (nCol == 4) {
        for (int j=0; j<4; j++)
            if (colGray4[j] == c) return j;
    }
    return 0;
}

// ==============================
// solveKmap(): main routine. Solves the K-map and produces the prime implicants,
// essential prime implicants, and minimum SOP expression.
KmapSolution *solveKmap(int numVar, const vector<int> minterms,
                        const vector<int> dontcares)
{
    // Create the result structure
    KmapSolution *sol = new KmapSolution;
    sol->numVar = numVar;
    
    // Determine K-map dimensions (rows and columns) from the number of variables
    int nRow, nCol;
    if (numVar == 2) { nRow = 2; nCol = 2; }
    else if (numVar == 3) { nRow = 2; nCol = 4; }
    else if (numVar == 4) { nRow = 4; nCol = 4; }
    else { nRow = 0; nCol = 0; }
    
    // Allocate and initialize the K-map array (each cell stored as a character)
    sol->kmap = new char*[nRow];
    for (int i = 0; i < nRow; i++) {
        sol->kmap[i] = new char[nCol];
        for (int j = 0; j < nCol; j++) {
            sol->kmap[i][j] = '0';
        }
    }
    
    // Fill in minterms and don't cares using Gray code ordering
    // Select the row/column Gray code tables according to the number of variables
    int *rg, *cg;
    if (numVar == 2) { rg = rowGray2; cg = colGray2; }
    else if (numVar == 3) { rg = rowGray3; cg = colGray3; }
    else { rg = rowGray4; cg = colGray4; }
    
    // Lambda that fills a single cell: splits a value into its row and column parts
    auto fillCell = [&](int val, char ch) {
        if (numVar == 2) {
            int rowBit = val & 1;      // Low bit -> row (b)
            int colBit = val >> 1;     // High bit -> column (a)
            sol->kmap[rg[rowBit]][cg[colBit]] = ch;
        }
        else if (numVar == 3) {
            int rowBit = val & 1;      // c -> row
            int colBit = (val >> 1);   // ab -> column (2 bits)
            sol->kmap[rg[rowBit]][cg[colBit]] = ch;
        }
        else { // numVar==4
            int rowBits = val & 3;     // cd -> row (2 bits)
            int colBits = (val >> 2);  // ab -> column (2 bits)
            sol->kmap[rg[rowBits]][cg[colBits]] = ch;
        }
    };
    
    // Fill minterms with '1' and don't cares with 'x'
    for (int val : minterms) {
        fillCell(val, '1');
    }
    for (int val : dontcares) {
        fillCell(val, 'x');
    }
    
    // ------------------------------
    // Find prime implicants by exhaustive search
    // Both '1' and 'x' cells may be grouped, but every group must contain at least one '1'
    // Rectangle sizes are powers of 2 (1, 2, 4) and cannot exceed the K-map size
    vector<int> possibleHeights, possibleWidths;
    for (int h = 1; h <= nRow; h *= 2)
        possibleHeights.push_back(h);
    for (int w = 1; w <= nCol; w *= 2)
        possibleWidths.push_back(w);
    
    // Store every cell index of each candidate group (as the original minterm value,
    // obtained through the inverse Gray code lookup)
    vector<vector<int>> candidateGroups;
    set<string> groupKeys; // String keys used to avoid duplicate candidate groups
    // Iterate over every starting cell (r, c) in K-map (Gray code) coordinates
    for (int r = 0; r < nRow; r++) {
        for (int c = 0; c < nCol; c++) {
            // Try every possible rectangle size
            for (int h : possibleHeights) {
                for (int w : possibleWidths) {
                    vector<int> groupIndices;
                    bool valid = true;   // Whether every cell in the rectangle is 1 or x
                    bool hasOne = false; // Whether the rectangle contains at least one 1
                    // Visit every cell in the rectangle, using modulo to handle wrap-around
                    for (int dr = 0; dr < h; dr++) {
                        for (int dc = 0; dc < w; dc++) {
                            int rr = (r + dr) % nRow;
                            int cc = (c + dc) % nCol;
                            char cell = sol->kmap[rr][cc];
                            if (cell != '1' && cell != 'x') {
                                valid = false;
                                break;
                            }
                            if (cell == '1') hasOne = true;
                            // Convert Gray code coordinates back to the original value:
                            int realRow = invertRowGray(nRow, rr);
                            int realCol = invertColGray(nCol, cc);
                            int val;
                            if (numVar == 2) {
                                val = (realCol << 1) + realRow; // (col*2 + row)
                            }
                            else if (numVar == 3) {
                                val = (realCol << 1) + realRow; // Same as the 2-variable case
                            }
                            else { // numVar==4
                                val = (realCol << 2) + realRow; // (col*4 + row)
                            }
                            groupIndices.push_back(val);
                        }
                        if (!valid) break;
                    }
                    // Add the rectangle as a candidate if it is valid and contains a 1
                    if (valid && hasOne) {
                        sort(groupIndices.begin(), groupIndices.end());
                        ostringstream oss;
                        for (int idx : groupIndices)
                            oss << idx << ",";
                        string key = oss.str();
                        if (groupKeys.find(key) == groupKeys.end()) {
                            groupKeys.insert(key);
                            candidateGroups.push_back(groupIndices);
                        }
                    }
                }
            }
        }
    }
    
    // ------------------------------
    // Filter: if a candidate group is a subset of another candidate group,
    // discard the smaller one
    int nCand = candidateGroups.size();
    vector<bool> isPrime(nCand, true);
    for (int i = 0; i < nCand; i++) {
        for (int j = 0; j < nCand; j++) {
            if (i == j) continue;
            if (candidateGroups[i].size() < candidateGroups[j].size() &&
                includes(candidateGroups[j].begin(), candidateGroups[j].end(),
                         candidateGroups[i].begin(), candidateGroups[i].end())) {
                isPrime[i] = false;
                break;
            }
        }
    }
    
    // ------------------------------
    // Convert the remaining candidate groups into Boolean expressions
    // [Key point]: call bitsToExpression on the WHOLE group (including 'x' don't cares),
    // but only include actual minterms in the coverage set
    vector<pair<string, set<int>>> implicants;
    for (int i = 0; i < nCand; i++) {
        if (!isPrime[i]) continue;
        const auto &grp = candidateGroups[i]; // grp holds all indices (both 1 and x)
        string expr = bitsToExpression(numVar, grp);
        if (expr.empty()) continue;
        // Only indices that are minterms go into the coverage set
        set<int> covered;
        for (int idx : grp) {
            if (find(minterms.begin(), minterms.end(), idx) != minterms.end())
                covered.insert(idx);
        }
        if (covered.empty()) continue;
        // Skip duplicates
        bool dup = false;
        for (auto &p : implicants) {
            if (p.first == expr) {
                dup = true;
                break;
            }
        }
        if (!dup) {
            implicants.push_back({expr, covered});
        }
    }
    // Store all prime implicant expressions in the solution structure
    for (auto &p : implicants) {
        sol->primes.push_back(p.first);
    }
    
    // ------------------------------
    // Next, compute the essential prime implicants and the minimum SOP
    // Build a coverage table: for each input minterm, list the implicants covering it
    map<int, vector<string>> coverage;
    for (int m : minterms) {
        for (auto &imp : implicants) {
            if (imp.second.find(m) != imp.second.end())
                coverage[m].push_back(imp.first);
        }
    }
    // If a minterm is covered by only one implicant, that implicant is essential
    set<string> essentialSet;
    for (auto &kv : coverage) {
        if (kv.second.size() == 1)
            essentialSet.insert(kv.second[0]);
    }
    // Store the essential implicants and record the minterms they cover
    set<int> coveredMinterms;
    for (auto &e : essentialSet) {
        sol->essentials.push_back(e);
        for (auto &imp : implicants) {
            if (imp.first == e) {
                for (int m : imp.second)
                    coveredMinterms.insert(m);
            }
        }
    }
    // Cover the remaining minterms with a greedy strategy
    set<int> uncovered;
    for (int m : minterms) {
        if (coveredMinterms.find(m) == coveredMinterms.end())
            uncovered.insert(m);
    }
    vector<string> selected(sol->essentials.begin(), sol->essentials.end());
    while (!uncovered.empty()) {
        string best;
        int bestCover = 0;
        for (auto &imp : implicants) {
            if (find(selected.begin(), selected.end(), imp.first) != selected.end())
                continue;
            int count = 0;
            for (int m : uncovered) {
                if (imp.second.find(m) != imp.second.end())
                    count++;
            }
            if (count > bestCover) {
                bestCover = count;
                best = imp.first;
            }
        }
        if (best.empty())
            break;
        selected.push_back(best);
        for (auto &imp : implicants) {
            if (imp.first == best) {
                for (int m : imp.second)
                    uncovered.erase(m);
            }
        }
    }
    sol->minimumSop = selected;
    
    return sol;
}

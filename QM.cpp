#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <algorithm>
#include <functional>
#include <iomanip>   // for setw(), setfill()

using namespace std;

// Structure: stores the result of the Quine-McCluskey algorithm, including the number of variables,
// the literal form of all prime implicants, and the set of minimum SOP solutions.
struct QmSolution {
	int numVar;
	vector<string> primes;
	vector<vector<string>> minimumSops;
};

/**
 * Convert integer n into a binary string of length numVar (MSB first).
 */
static string parseString(int n, int numVar) {
	string s(numVar, '0');
	for (int i = numVar - 1; i >= 0; --i) {
		s[i] = (n & 1) ? '1' : '0';
		n >>= 1;
	}
	return s;
}

/**
 * Convert a bit string into its literal (product-term) representation.
 */
static string toLiteral(const string &bits) {
	string lit;
	for (int i = 0; i < (int)bits.size(); ++i) {
		char var = char('a' + i);
		if (bits[i] == '1') lit.push_back(var);
		else if (bits[i] == '0') { lit.push_back(var); lit.push_back('\''); }
		// '-' is skipped
	}
	return lit;
}

/**
 * Main algorithm: performs Quine-McCluskey minimization.
 */
QmSolution* solveQm(int numVar,
	vector<int> minterms,
	vector<int> dontcares,
	bool verbose)
{
	QmSolution *sol = new QmSolution;
	sol->numVar = numVar;

	struct Imp { string bits; set<int> covered; bool used = false; };
	vector<Imp> curr;
	set<int> mtset(minterms.begin(), minterms.end()), dcset(dontcares.begin(), dontcares.end());
	auto addInit = [&](int x) { Imp ip; ip.bits = parseString(x, numVar); ip.covered.insert(x); curr.push_back(ip); };
	for (int m : minterms) addInit(m);
	for (int d : dontcares) addInit(d);

	struct Entry { int ones; string bits; char mark; vector<int> idx; };
	vector<vector<Entry>> logEntries;  // Grouping log for each column
	vector<pair<string, vector<int>>> primes;
	int column = 1;

	while (!curr.empty()) {
		// 1. Determine which implicants can be merged
		vector<bool> willMerge(curr.size(), false);
		for (int i = 0; i < (int)curr.size(); ++i) for (int j = i + 1; j < (int)curr.size(); ++j) {
			const string &A = curr[i].bits, &B = curr[j].bits;
			int diff = 0;
			for (int k = 0; k < numVar; ++k) {
				if (A[k] != B[k]) {
					if (A[k] == '-' || B[k] == '-') { diff = -1; break; }
					++diff;
				}
			}
			if (diff == 1) willMerge[i] = willMerge[j] = true;
		}

		// 2. Collect grouping information
		vector<Entry> entries;
		map<string, int> freq;
		for (auto &ip : curr) freq[ip.bits]++;
		for (int idx = 0; idx < (int)curr.size(); ++idx) {
			auto &ip = curr[idx]; Entry e;
			e.bits = ip.bits;
			e.ones = count(ip.bits.begin(), ip.bits.end(), '1');
			if (freq[ip.bits] > 1) e.mark = 'x';
			else if (ip.covered.size() == 1) {
				int x = *ip.covered.begin();
				if (!mtset.count(x) && dcset.count(x)) e.mark = 'd';
				else if (mtset.count(x) && willMerge[idx]) e.mark = 'v';
				else e.mark = ' ';
			}
			else e.mark = willMerge[idx] ? 'v' : ' ';
			e.idx.assign(ip.covered.begin(), ip.covered.end());
			sort(e.idx.begin(), e.idx.end());
			entries.push_back(e);
		}
		if (verbose) logEntries.push_back(entries);

		// 3. Merge
		map<string, set<int>> nextMap;
		for (int i = 0; i < (int)curr.size(); ++i) for (int j = i + 1; j < (int)curr.size(); ++j) {
			const string &A = curr[i].bits, &B = curr[j].bits; int diff = 0, pos = -1;
			for (int k = 0; k < numVar; ++k) if (A[k] != B[k]) {
				if (A[k] == '-' || B[k] == '-') { diff = -1; break; } diff++; pos = k;
			}
			if (diff == 1) {
				curr[i].used = curr[j].used = true;
				string C = A; C[pos] = '-';
				auto &cov = nextMap[C];
				cov.insert(curr[i].covered.begin(), curr[i].covered.end());
				cov.insert(curr[j].covered.begin(), curr[j].covered.end());
			}
		}
		// Implicants that were not merged are prime
		for (auto &ip : curr) if (!ip.used) {
			vector<int> cov;
			for (int x : ip.covered) if (mtset.count(x)) cov.push_back(x);
			if (!cov.empty()) primes.emplace_back(ip.bits, cov);
		}
		if (nextMap.empty()) break;
		curr.clear();
		for (auto &kv : nextMap) {
			Imp ip; ip.bits = kv.first; ip.covered = kv.second; curr.push_back(ip);
		}
		++column;
	}

	// 4. Print the merging process of each column
	if (verbose) {
		for (int h = 0; h < (int)logEntries.size(); ++h) {
			auto &ents = logEntries[h];
			if (ents.empty()) continue;
			cout << "\n\n\n";
			cout << "====================\n";
			cout << "Column " << (h + 1) << "\n";
			cout << "====================\n";
			sort(ents.begin(), ents.end(), [&](auto &a, auto &b) {
				if (a.ones != b.ones) return a.ones < b.ones;
				return a.bits < b.bits;
			});
			int lastOnes = ents[0].ones;
			for (auto &e : ents) {
				if (e.ones != lastOnes) { cout << "----------------\n"; lastOnes = e.ones; }
				cout << e.mark << " " << setw(numVar) << setfill('0') << e.bits << ": ";
				for (int v : e.idx) cout << v << ",";
				cout << "\n";
			}
			cout << "----------------\n";
		}
	}

	// 5. Convert prime implicants to literal form
	QmSolution* solptr = sol;
	for (auto &p : primes) solptr->primes.push_back(toLiteral(p.first));

	// 6. Build the prime implicant chart
	int P = primes.size(), M = minterms.size();
	vector<vector<bool>> chart(P, vector<bool>(M, false));
	for (int i = 0; i < P; ++i) for (int j = 0; j < M; ++j)
		if (find(primes[i].second.begin(), primes[i].second.end(), minterms[j]) != primes[i].second.end())
			chart[i][j] = true;

	// 7. Print the chart
	if (verbose) {
		cout << "\n====================\nPrime Implicant Chart\n====================\n";
		cout << setw(numVar * 2 + 1) << setfill(' ') << "|";
		for (int j = 0; j < M; ++j) cout << "|" << setw(3) << setfill('0') << minterms[j];
		cout << "|\n";
		cout << string(numVar * 2 + 1, '-');
		for (int j = 0; j < M; ++j) cout << "----";
		cout << "\n";
		for (int i = 0; i < P; ++i) {
			string lit = toLiteral(primes[i].first);
			cout << setw(numVar * 2) << setfill(' ') << lit;
			for (int j = 0; j < M; ++j) cout << "| " << (chart[i][j] ? 'x' : ' ') << " ";
			cout << "|\n";
		}
		cout << "\n";
	}

	// 8. Essential prime implicants & Petrick's method
	vector<bool> essential(P, false);
	set<int> uncovered; for (int j = 0; j < M; ++j) uncovered.insert(j);
	bool found;
	do {
		found = false;
		for (int j : vector<int>(uncovered.begin(), uncovered.end())) {
			int cnt = 0, who = -1;
			for (int i = 0; i < P; ++i) if (!essential[i] && chart[i][j]) { cnt++; who = i; }
			if (cnt == 1) {
				essential[who] = true;
				for (int k = 0; k < M; ++k) if (chart[who][k]) uncovered.erase(k);
				found = true; break;
			}
		}
	} while (found);
	vector<int> essIdx;
	for (int i = 0; i < P; ++i) if (essential[i]) essIdx.push_back(i);
	if (uncovered.empty()) {
		vector<string> one;
		for (int i : essIdx) one.push_back(solptr->primes[i]);
		solptr->minimumSops.push_back(one);
		return solptr;
	}
	// Petrick's method
	vector<vector<int>> sumOfPrimes;
	for (int j : uncovered) {
		vector<int> covers;
		for (int i = 0; i < P; ++i) if (chart[i][j]) covers.push_back(i);
		if (sumOfPrimes.empty()) for (int pi : covers) sumOfPrimes.push_back({ pi });
		else {
			vector<vector<int>> nextSum;
			for (auto &term : sumOfPrimes) for (int pi : covers) {
				auto nt = term; nt.push_back(pi);
				sort(nt.begin(), nt.end()); nt.erase(unique(nt.begin(), nt.end()), nt.end());
				nextSum.push_back(nt);
			}
			sort(nextSum.begin(), nextSum.end()); nextSum.erase(unique(nextSum.begin(), nextSum.end()), nextSum.end());
			sumOfPrimes.swap(nextSum);
		}
	}
	size_t bestSize = sumOfPrimes[0].size();
	for (auto &t : sumOfPrimes) bestSize = min(bestSize, t.size());
	for (auto &t : sumOfPrimes) if (t.size() == bestSize) {
		vector<string> one;
		for (int i : essIdx) one.push_back(solptr->primes[i]);
		for (int pi : t) one.push_back(solptr->primes[pi]);
		solptr->minimumSops.push_back(one);
	}
	return solptr;
}

// Helper that prints the result (same format as the original version)

static void printSolution(const QmSolution *sol) {
	cout << "numVar = " << sol->numVar << "\n";
	cout << "primes = [";
	for (size_t i = 0; i < sol->primes.size(); ++i) { cout << "\"" << sol->primes[i] << "\"" << (i + 1 < sol->primes.size() ? ", " : ""); }
	cout << "]\n\nminimumSops = [\n";
	for (auto&sop : sol->minimumSops) { cout << "  ["; for (size_t j = 0; j < sop.size(); ++j)cout << "\"" << sop[j] << "\"" << (j + 1 < sop.size() ? ", " : ""); cout << "]\n"; }cout << "]\n";
}





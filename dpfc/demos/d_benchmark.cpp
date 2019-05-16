/*
 * A program comparing the accuracy and speed of the histogram algorithms in the dpfc package.
 */

#include "dpfc.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>


using namespace dpfc;
using std::accumulate;
using std::cout;
using std::endl;
using std::get;
using std::ifstream;
using std::make_unique;
using std::ofstream;
using std::pair;
using std::setw;
using std::sort;
using std::string;
using std::stringstream;
using std::tuple;
using std::unique_ptr;
using std::vector;

using clk = std::chrono::high_resolution_clock;
using hist_pairs_t = vector<pair<unsigned long, unsigned long>>;
using meta = tuple<mpz_class, unsigned long, string>;
using vmpz = vector<mpz_class>;

// privacy parameters
mpq_class eps("1/10");
mpq_class delta("1/340282366920938463463374607431768211456");

// scale the delta parameter to get (eps, delta)-differential privacy from statistical distance
mpq_class tau = delta / 211;

// the number of times to run each algorithm on each dataset
int trials = 10;

// algorithms used in benchmarking (tau <= delta / (100 * (1 + e^eps)))
vector<string> algorithms {
	"basic-g", // Basic() Geo()
	"stab-g",  // Stability() Geo()
	"stab-f",  // Stability() Fast(gamma=delta/2)
	"pure-g",  // PureSparse(beta=1/80, tvd=100*tau) Fast(gamma=1/(160*|X|))
	"sp-g",    // Sparse(delta=99*tau, zeta=0, tvd=tau) Geo()
	"sp-f",    // Sparse(delta=99*tau, zeta=0, tvd=tau) Fast(gamma=1/(40*|X|))
	"sp-f-z",  // Sparse(delta=9*tau, zeta=90*tau, tvd=tau) Fast(gamma=zeta/|X|)
	"sp-f-k",  // Sparse(delta=50*tau, zeta=49*tau, tvd=tau) Fast(gamma=zeta/|X|)
	"sp-f-kp", // Sparse(delta=50*tau, zeta=49*tau, tvd=tau) Fast(gamma=zeta/|X|) Min(1.5/eps*log|X|)
	"ff-f",    // Compact(eps=9*eps/10, eps2=eps/10, tvd=100*tau) Fast(gamma=1/(40*|X|))
};

// the list of datasets to use for benchmarking
mpz_class X2_32("4294967296");
mpz_class X2_64("18446744073709551616");
vector<meta> dataset_meta {
	{X2_32,    100, "all-in-one"},
	{X2_64,    100, "all-in-one"},
	{ 4096,   1000, "all-in-one"},
	{X2_32,   1000, "all-in-one"},
	{X2_64,   1000, "all-in-one"},
	{X2_32,  10000, "all-in-one"},
	{X2_64,  10000, "all-in-one"},
	{X2_32, 100000, "all-in-one"},
	{X2_64, 100000, "all-in-one"},
	{X2_32,   1000, "exponential"},
	{X2_64,   1000, "exponential"},
	{X2_32,  10000, "exponential"},
	{X2_64,  10000, "exponential"},
	{X2_32, 100000, "exponential"},
	{X2_64, 100000, "exponential"},
	{X2_32,  10000, "uniform-600"},
	{X2_32, 100000, "uniform-600"},
	{X2_32,    500, "unique"},
	{X2_32,   5000, "unique"},
	{X2_32, 100000, "uniform-2000"},
	{X2_64, 100000, "uniform-2000"},
	{X2_32,  50000, "uniform-2500"},
	{X2_32, 100000, "uniform-3500"},
	{X2_64, 100000, "uniform-3500"},
};

// Create a histogram of n elements based on the given type string.
histogram_t datasetGenerator(unsigned long n, string type) {
	histogram_t hist;
	if (type == "all-in-one") {
		hist.emplace_back(1, n);
	}
	else if (type == "exponential") {
		unsigned long c;
		for (mpz_class label(1); n > 0; n -= c, label++) {
			hist.emplace_back(label, c = (n > 1) ? n / 2 : 1);
		}
	}
	else if (!type.rfind("uniform-", 0)) {
		unsigned long count = std::stoi(type.substr(8)), c;
		for (mpz_class label(1); n > 0; n -= c, label++) {
			hist.emplace_back(label, c = (count < n) ? count : n);
		}
	}
	else if (type == "unique") {
		for (mpz_class label(1); n > 0; n--, label++) {
			hist.emplace_back(label, 1);
		}
	}
	return hist;
}

// An iterator for our datasets where the internal storage representation is a histogram.
class HistogramIterator {
	public:
		const histogram_t* const hist;
		unsigned long count = 0;
		unsigned long index = 0;

	public:
		HistogramIterator(histogram_t* hist) : hist(hist) { }

		void operator ++(int) {
			if (++this->count == (*hist)[this->index].count) {
				this->count = 0;
				this->index++;
			}
		}

		mpz_class operator *() {
			return (*hist)[this->index].label;
		}

		// only returns true after having iterated over the entire histogram
		bool operator !=(const HistogramIterator& that) {
			return this->index != this->hist->size();
		}
};

// Sort a histogram by bin label.
void histogramSort(histogram_t& hist) {
	sort(hist.begin(), hist.end(), [](const bin_t& a, const bin_t& b) {
			return a.label < b.label;
		});
}

// Given two histograms sorted by label, this function returns a vector of pairs for each label with
// either count nonzero holding both counts for that label.
hist_pairs_t histogramZip(histogram_t* a, histogram_t* b) {
	hist_pairs_t pairs;

	mpz_class la(0), lb(0);
	size_t ia = 0, ib = 0;
	unsigned long ca, cb;
	int k = 0;

	while (ia < a->size() || ib < b->size()) {
		if (ia == a->size()) {
			pairs.push_back({0, (*b)[ib].count});
			ib++;
			continue;
		}
		if (ib == b->size()) {
			pairs.push_back({(*a)[ia].count, 0});
			ia++;
			continue;
		}

		la = (*a)[ia].label;
		ca = (*a)[ia].count;
		lb = (*b)[ib].label;
		cb = (*b)[ib].count;

		if (lb < la) {
			pairs.push_back({0, cb});
			ib++;
		}
		else if (la < lb) {
			pairs.push_back({ca, 0});
			ia++;
		}
		else if (lb == la) {
			pairs.push_back({ca, cb});
			ia++;
			ib++;
		}
	}
	return pairs;
}

// Returns a triple of the l-infinity error between the pairs. The first entry contains the error
// over pairs where the second value is zero, the second entry is over pairs where the second value
// is larger than zero and the final entry is over all pairs.
vmpz linferror(const hist_pairs_t& pairs) {
	vmpz dist = {0, 0, 0};
	mpz_class curr;
	for (const auto c : pairs) {
		curr = (c.first > c.second) ? c.first - c.second : c.second - c.first;
		if (curr > dist[c.second > 0]) {
			dist[c.second > 0] = curr;
		}
	}
	dist[2] = dist[0] > dist[1] ? dist[0] : dist[1];
	return dist;
}

// Returns a triple of the l1 error between the pairs. The structure of the returned value is the
// same as linferror.
vmpz l1error(const hist_pairs_t& pairs) {
	vmpz dist = {0, 0, 0};
	for (const auto c : pairs) {
		dist[c.second != 0] += (c.first > c.second) ? c.first - c.second : c.second - c.first;
	}
	dist[2] = dist[0] + dist[1];
	return dist;
}

// Sum over of the specified index in the given vector.
mpz_class sum(const vector<vmpz>& vec, int index) {
	mpz_class acc(0);
	for(auto it = vec.cbegin(); it != vec.cend(); it++) {
		acc += (*it)[index];
	}
	return acc;
}

// Benchmark all algorithms on the given dataset.
void runConfiguration(mpz_class X, unsigned long n, string type) {
	stringstream ssname;
	ssname << "n" << n << " " << type << " X2^" << mpz_clog2(X.get_mpz_t()) << ".txt";

	// skip if the file already exists
	if (ifstream(ssname.str())) {
		return;
	}
	cout << ssname.str() << endl;

	stringstream outStream;

	// create the dataset
	histogram_t trueHist = datasetGenerator(n, type);
	HistogramIterator ds(&trueHist);

	// print out header information
	outStream << "X=" << X << endl << "n=" << n << endl << "type=" << type << endl << endl;
	outStream << setw(7) << "alg";
	outStream << setw(12) << "L1";
	outStream << setw(12) << "L1 (>0)";
	outStream << setw(12) << "L1 (=0)";
	outStream << setw(12) << "Linf (>0)";
	outStream << setw(12) << "Linf (=0)";
	outStream << setw(12) << "time (ms)" << endl;

	histogram_t hist;
	unique_ptr<CountSample> cs;

	for (auto alg : algorithms) {

		// drop algorithms based on parameters
		if (alg == "basic-g" && (n > 1000 || X > 65536)) {
			continue;
		}
		else if (alg == "pure-g" && (n > 100)) {
			continue;
		}
		else if (alg == "sp-g" && (n > 1000)) {
			continue;
		}
		else if (alg == "sp-f" && (n > 50000)) {
			continue;
		}
		else if (alg == "ff-f" && (n > 500)) {
			continue;
		}
		else if (alg == "sp-f-z" && (n > 10000)) {
			continue;
		}

		vector<vmpz> l1s, linfs;
		vector<unsigned long> speeds;

		for (int i = 0; i < trials; i++) {
			auto time_start = clk::now();

			if (alg == "basic-g") {
				cs = make_unique<GeoSample>(n, eps);
				hist = BasicHistogram(X, ds, ds, *cs);
			}
			else if (alg == "stab-g") {
				cs = make_unique<GeoSample>(n, eps);
				hist = StabilityHistogram(X, ds, ds, *cs, delta);
				for (auto bin : hist) {
					cout << bin.count << endl;
				}
			}
			else if (alg == "stab-f") {
				cs = make_unique<FastSample>(n, eps, delta / 2);
				hist = StabilityHistogram(X, ds, ds, *cs, delta);
			}
			else if (alg == "pure-g") {
				cs = make_unique<FastSample>(n, eps, mpq_class(1, 160 * X));
				hist = PureSparseHistogram(X, ds, ds, mpq_class(1, 80), *cs, 100 * tau);
				histogramSort(hist);
			}
			else if (alg == "sp-g") {
				cs = make_unique<GeoSample>(n, eps);
				hist = SparseHistogram(X, ds, ds, 99 * tau, 0, *cs, tau);
				histogramSort(hist);
			}
			else if (alg == "sp-f") {
				cs = make_unique<FastSample>(n, eps, mpq_class(1, 40 * X));
				hist = SparseHistogram(X, ds, ds, 99 * tau, 0, *cs, tau);
				histogramSort(hist);
			}
			else if (alg == "sp-f-z") {
				mpq_class gamma(90 * tau / X);
				mpq_zinv_round(gamma.get_mpq_t(), gamma.get_mpq_t());
				cs = make_unique<FastSample>(n, eps, gamma);
				hist = SparseHistogram(X, ds, ds, 9 * tau, 90 * tau, *cs, tau);
				histogramSort(hist);
			}
			else if (alg == "sp-f-k") {
				mpq_class gamma(49 * tau / X);
				mpq_zinv_round(gamma.get_mpq_t(), gamma.get_mpq_t());
				cs = make_unique<FastSample>(n, eps, gamma);
				hist = SparseHistogram(X, ds, ds, 50 * tau, 49 * tau, *cs, tau);
				histogramSort(hist);
			}
			else if (alg == "sp-f-kp") {
				mpq_class gamma(49 * tau / X);
				mpq_zinv_round(gamma.get_mpq_t(), gamma.get_mpq_t());
				cs = make_unique<FastSample>(n, eps, gamma);
				hist = SparseHistogram(X, ds, ds, 50 * tau, 49 * tau, *cs, tau);
				histogramSort(hist);

				mpz_class cap_mpz(3 * (unsigned long) mpz_clog2(X.get_mpz_t()) * eps.get_den());
				cap_mpz /= eps.get_num();
				unsigned long cap = cap_mpz.get_ui();
				for (auto& bin : hist) {
					if (bin.count < cap) {
						bin.count = 0;
					}
				}
			}
			else if (alg == "ff-f") {
				cs = make_unique<FastSample>(n, 9 * eps / 10, mpq_class(1, 40 * X));
				auto chist = CompactHistogram(X, ds, ds, *cs, eps / 10, 100 * tau);
			}
			else {
				continue;
			}
			auto time_end = clk::now();

			auto pairs = histogramZip(&hist, &trueHist);
			l1s.push_back(l1error(pairs));
			linfs.push_back(linferror(pairs));
			auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start);
			speeds.push_back(dur.count());
		}

		outStream << setw(7) << alg;
		if (alg == "ff-f") {
			outStream << setw(12) << "-";
			outStream << setw(12) << "-";
			outStream << setw(12) << "-";
			outStream << setw(12) << "-";
			outStream << setw(12) << "-";
		}
		else {
			outStream << setw(12) << sum(l1s, 2) / trials;
			outStream << setw(12) << sum(l1s, 1) / trials;
			outStream << setw(12) << sum(l1s, 0) / trials;
			outStream << setw(12) << sum(linfs, 1) / trials;
			outStream << setw(12) << sum(linfs, 0) / trials;
		}
		outStream << setw(12) << accumulate(speeds.begin(), speeds.end(), 0) / trials << endl;
	}
	outStream << endl;

	ofstream outFile;
	outFile.open(ssname.str());
	outFile << outStream.str();
	outFile.close();
}


int main(int argc, char **argv) {
	for (auto ds : dataset_meta) {
		runConfiguration(get<0>(ds), get<1>(ds), get<2>(ds));
	}
}


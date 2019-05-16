/*
 * A demo showing how to make a noisy histogram for datasets with string labels.
 */

#include "dpfc.hpp"

#include <chrono>
#include <iomanip>
#include <ostream>
#include <random>
#include <string>
#include <vector>


using dpfc::histogram_t;
using std::chrono::system_clock;
using std::cout;
using std::endl;
using std::default_random_engine;
using std::setw;
using std::string;
using std::vector;

// length of string labels
static const int LEN = 4;

// map a string to a number
mpz_class to_num(string s) {
	mpz_class num;

	for (int i = 0; i < LEN; i++) {
		num *= 26;
		num += (int) s[i] - 97;
	}

	return num + 1;
}

// map a number to a lowercase alphabetic string (the inverse of to_num)
string to_label(mpz_class num) {
	num -= 1;
	mpz_class mod;
	string s(LEN, '\0');

	for (int i = LEN - 1; i >= 0; i--) {
		mod = num % 26;
		s[i] = 97 + mod.get_ui();
		num /= 26;
	}
	return s;
}

// Defines a class for creating an iterator over a dataset of strings.
// To work with the histogram algorithms in the dpfc package, the iterator must support
// "postfix increment", "indirection" and "not equal to".
class StringDatasetIterator {
	private:
		vector<string>* dataset;
		unsigned long index = 0;

	public:
		StringDatasetIterator(vector<string>* dataset) {
			this->dataset = dataset;
		}

		void operator ++(int) {
			this->index++;
		}

		// get the current element and convert it to a number
		mpz_class operator *() {
			return to_num((*dataset)[this->index]);
		}

		// returns true until the index reaches the end of the dataset
		bool operator !=(const StringDatasetIterator& that) {
			return this->index < this->dataset->size();
		}
};

int main(int argc, char **argv) {
	// determine the number of labels
	mpz_class X;
	mpz_ui_pow_ui(X.get_mpz_t(), 26, LEN);

	// construct a dataset of random LEN character strings
	unsigned int seed = system_clock::now().time_since_epoch().count();
	default_random_engine gen(seed);

	vector<string> dataset;
	string s;
	for (int i = 0; i < 20; i++) {
		// add fifty copies of each string to the dataset
		s = to_label(gen());
		for (int j = 0; j < 50; j++) {
			dataset.push_back(s);
		}
	}

	// create an iterator that returns a number representation for each string upon indirection
	StringDatasetIterator it(&dataset);

	// run the noisy histogram algorithm
	unsigned long n = 1000;
	mpq_class eps(1, 4), delta(1, 64);
	GeoSample gs(n, eps);

	histogram_t hist = StabilityHistogram(X, it, it, gs, delta);

	// display noisy counts
	cout << "Noisy Histogram" << endl << endl;
	cout << string(LEN > 5 ? LEN - 3 : 2, ' ') << "label  count" << endl;
	for (const auto& bin : hist) {
		cout << string(LEN <= 5 ? 7 - LEN : 2, ' ') << to_label(bin.label);
		cout << setw(7) << bin.count << endl;
	}
}


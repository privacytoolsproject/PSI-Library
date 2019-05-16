/*
 * Implementation of the histogram algorithms from "Differential Privacy on Finite Computers"
 * Section 5.
 */

#ifndef DPFC_HIST_SIMPLE_HPP
#define DPFC_HIST_SIMPLE_HPP

#include "dpfc/count.hpp"

#include "gmpxx.h"
#include "sodium.h"

#include <map>
#include <stdexcept>
#include <utility>
#include <vector>


// Datasets are passed as iterators of positive integers (mpz_class) with a maximum value of X. The
// histogram algorithms in this source file (BasicHistogram, StabilityHistogram and
// PartialBasicHistogram) take as input X defining the data universe, a dataset of elements in
// {1, ..., X} and a dpfc::CountSample object to generate noisy counts.

namespace dpfc {
	// data type definitions
	struct bin_t {
		mpz_class label;
		unsigned long count;

		bin_t() : label(0), count(0) {}
		bin_t(const mpz_class& label, unsigned long count) : label(label), count(count) {}
		bin_t(const std::pair<mpz_class, unsigned long>& p) : label(p.first), count(p.second) {}
	};
	using histogram_t = std::vector<bin_t>;
	using histogram_map_t = std::map<mpz_class, unsigned long>;

	// (Algorithm 5.2) This is an implementation of BasicHistogram that constructs a histogram only
	// for the bins with nonzero true count. The counts are perturbed in the output histogram by
	// applying the given counting query algorithm (CountSample::sample). Note that in addition to
	// having only bins for nonzero counts, privacy is possibly violated in structure of the returned
	// map.
	template <class Iterator>
	histogram_map_t PartialBasicHistogram(const mpz_class& X, Iterator dataset_begin,
			Iterator dataset_end, const CountSample& cs);
}

// (Algorithm 5.2) BasicHistogram returns a noisy histogram with a count for each bin that is
// sampled from the given counting query algorithm using the true count as the center in
// CountSample::sample.
// Warning: This algorithm should not be applied to datasets drawn from a large universe as it
// returns a count for each bin. No checks are done prior to instantiating the histogram and the
// process will fail if the universe is too large.
template <class Iterator>
dpfc::histogram_t BasicHistogram(const mpz_class& X, Iterator dataset_begin, Iterator dataset_end,
		const dpfc::CountSample& cs);

// (Algorithm 5.4) StabilityHistogram adds noise (using CountSample::sample) to each bin with
// nonzero true count. Noisy counts below a threshold (derived from the privacy parameter delta and
// the counting query algorithm used) are filtered out before releasing the noisy histogram.
template <class Iterator>
dpfc::histogram_t StabilityHistogram(const mpz_class& X, Iterator dataset_begin,
		Iterator dataset_end, const dpfc::CountSample& cs, const mpq_class& delta);


// --- Implementation ---

namespace dpfc {
	template <class Iterator>
	histogram_map_t PartialBasicHistogram(const mpz_class& X, Iterator dataset_begin,
			Iterator dataset_end, const CountSample& cs) {
		if (X <= 0) {
			throw std::logic_error("The data universe parameter X must be a positive integer.");
		}

		histogram_map_t histogram;
		mpz_class label;
		unsigned long count = 0;

		// iterate over the dataset and maintain counts for the valid labels
		for (; dataset_begin != dataset_end; dataset_begin++) {
			label = *dataset_begin;
			if (label >= 1 && label <= X) {
				histogram[label] = histogram.count(label) ? histogram[label] + 1 : 1;
			}
			count++;
		}

		if (count != cs.get_n()) {
			throw std::logic_error("The size of the dataset must match n from the given CountSample "
					"object.");
		}

		// replace the true counts with noisy counts
		for (auto& bin : histogram) {
			bin.second = cs.sample(bin.second);
		}

		return histogram;
	}
}

template <class Iterator>
dpfc::histogram_t BasicHistogram(const mpz_class& X, Iterator dataset_begin, Iterator dataset_end,
		const dpfc::CountSample& cs) {
	if (!X.fits_uint_p()) {
		throw std::runtime_error("The given data universe is too large for BasicHistogram.");
	}

	// add noise to nonzero counts
	auto histTree = PartialBasicHistogram(X, dataset_begin, dataset_end, cs);

	// add noise to all zero counts and convert to type histogram_t
	dpfc::histogram_t histogram(X.get_ui());
	size_t i = 0;
	for (mpz_class label = 1; label <= X; label++, i++) {
		histogram[i].label = label;
		histogram[i].count = histTree.count(label) ? histTree[label] : cs.sample(0);
	}

	return histogram;
}

template <class Iterator>
dpfc::histogram_t StabilityHistogram(const mpz_class& X, Iterator dataset_begin,
		Iterator dataset_end, const dpfc::CountSample& cs, const mpq_class& delta) {
	if (delta <= 0 || delta >= 1) {
		throw std::logic_error("The privacy parameter delta must be in the range (0,1).");
	}

	// add noise to nonzero counts
	auto histTree = PartialBasicHistogram(X, dataset_begin, dataset_end, cs);

	// filter out small counts and convert to type histogram_t
	dpfc::histogram_t hist;
	unsigned long b = cs.threshold(delta);
	for (auto it = histTree.cbegin(); it != histTree.cend(); it++) {
		if (it->second > b) {
			hist.emplace_back(*it);
		}
	}
	return hist;
}

#endif


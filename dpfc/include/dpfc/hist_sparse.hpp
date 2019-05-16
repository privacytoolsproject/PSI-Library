/*
 * Header file for the algorithms from "Differential Privacy on Finite Computers" Section 6.
 */

#ifndef DPFC_HIST_SPARSE_HPP
#define DPFC_HIST_SPARSE_HPP

#include "dpfc/count.hpp"
#include "dpfc/hist_simple.hpp"
#include "dpfc/utils.hpp"

#include "gmpxx.h"
#include "NTL/ZZ.h"
#include "NTL/ZZX.h"
#include "sodium.h"

#include <algorithm>
#include <stdexcept>
#include <vector>


// Some functions in this module will generate random samples than can only be generated
// approximately. These functions will have an extra parameter, error_tvd of type mpq_class, to
// control the statistical distance of their output distributions from their true distributions.

namespace dpfc {
	// (Algorithm 6.16) ApproxConvExp returns an approximation to the first t terms of the i-fold
	// convolution of the polynomial a (intermediate steps are rounded to multiples of 1/s).
	NTL::ZZX ApproxConvExp(const NTL::ZZX& a, const mpz_class& i, const NTL::ZZ& s, unsigned long t);

	// (Algorithm 6.18) ApproxBinSample returns a sample from the min of t and an approximate
	// (error bounded by a function of s and the other parameters) binomial distribution with
	// parameters m and p. u is the randomness used by this function and assumed to be uniform over
	// {1, ..., s}.
	unsigned long ApproxBinSample(const mpz_class& m, const mpq_class& p, const mpz_class& s,
			const unsigned long t, const mpz_class& u);

	// (Algorithm 6.20) ApproxOrdSample returns an approximation (error bounded by a function of s
	// and the other parameters) of the counts of the n+1 heaviest noisy bins from a set of m bins
	// with a true count of 0.
	std::vector<unsigned long> ApproxOrdSample(const mpz_class& m, const mpz_class& s,
			const CountSample& cs, const mpq_class& error_tvd);

	// ApproxOrdSampleClamp returns an approximation of the algorithm that invokes Algorithm 6.20
	// and then prior to releasing counts clamps each to the range [0, k].
	std::vector<unsigned long> ApproxOrdSampleClamp(const mpz_class& m, const mpz_class& s,
			const CountSample& cs, unsigned long k, const mpq_class& error_tvd);

	// (Algorithm A.3) DistinctSample returns a list of r unique labels drawn almost uniformly at
	// random from the data universe X such that none are a label in the given histogram A.
	// (Note: You must first initialize Libsodium.)
	// (Note: Behavior is undefined if there are not enough labels to sample from.)
	std::vector<mpz_class> DistinctSample(const mpz_class& X, const histogram_map_t& A,
			unsigned long r, const mpq_class& error_tvd);

	// Return a random histogram drawn from data universe X with at most n bins as defined in
	// PureSparseHistogram (Algorithm 6.9).
	// (Note: You must first initialize Libsodium.)
	histogram_t RandomHistogram(const mpz_class& X, unsigned long n, const mpq_class& error_tvd);
}

// SparseHistogram returns a noisy sparse histogram. It creates a noisy histogram on the bins with
// nonzero true count (dpfc::PartialBasicHistogram) and for the remaining bins it simulates their
// noisy heavy hitters. The algorithm returns the heaviest (at most n) of the merged set of bins.
//
// By taking zeta equal to 0, this algorithm is Algorithm 6.21. For zeta > 0, the algorithm uses
// ApproxOrdSampleClamp with k set such that the true order statistics when compared to the true
// order statistics clamped to k have statistical distance at most delta from each other.
//
// Overall, we have (epsilon, (1 + e^epsilon) * (delta + zeta + error_tvd))-differential privacy.
template <class Iterator>
dpfc::histogram_t SparseHistogram(const mpz_class& X, Iterator dataset_begin, Iterator dataset_end,
		const mpq_class& delta, const mpq_class& zeta, const dpfc::CountSample& cs,
		const mpq_class& error_tvd);

// (Algorithm 6.9) PureSparseHistogram mixes between SparseHistogram and dpfc::RandomHistogram with
// probability beta1.
template <class Iterator>
dpfc::histogram_t PureSparseHistogram(const mpz_class& X, Iterator dataset_begin,
		Iterator dataset_end, const mpq_class& beta, const dpfc::CountSample& cs,
		const mpq_class& error_tvd);


// --- Implementation ---

template <class Iterator>
dpfc::histogram_t SparseHistogram(const mpz_class& X, Iterator dataset_begin, Iterator dataset_end,
		const mpq_class& delta, const mpq_class& zeta, const dpfc::CountSample& cs,
		const mpq_class& error_tvd) {
	if (sodium_init() == -1) {
		throw std::runtime_error("Unable to initialize Libsodium.");
	}

	if (delta <= 0 || delta >= 1) {
		throw std::logic_error("The privacy parameter delta must be in the range (0,1).");
	}
	if (zeta < 0 || zeta >= 1) {
		throw std::logic_error("The clamping parameter zeta must be in the range [0,1).");
	}
	if (error_tvd <= 0 || error_tvd >= 1) {
		throw std::logic_error("The randomness approximation parameter error_tvd must be in the range"
				"(0,1).");
	}

	unsigned long n = cs.get_n();
	if (X < 2*n + 1) {
		throw std::logic_error("The data universe size is too small. Use BasicHistogram instead.");
	}

	// set up parameters
	unsigned long k = n;
	if (zeta > 0) {
		mpq_class u = cs.get_d() * (1 - zeta / X);
		mpz_class threshold;
		mpz_cdiv_q(threshold.get_mpz_t(), u.get_num().get_mpz_t(), u.get_den().get_mpz_t());
		k = cs.sample(0, threshold);
	}
	mpq_class deltaZInv;
	dpfc::mpq_zinv_round(deltaZInv.get_mpq_t(), delta.get_mpq_t());
	mpz_class s = k * X * (n + 2) * deltaZInv.get_den();

	// main algorithm
	auto histTree = dpfc::PartialBasicHistogram(X, dataset_begin, dataset_end, cs);

	auto labels = dpfc::DistinctSample(X, histTree, n + 1, error_tvd / 2);

	mpz_class m;
	dpfc::mpz_set_ull(m.get_mpz_t(), histTree.size());
	m = X - m;
	auto counts = dpfc::ApproxOrdSampleClamp(m, s, cs, k, error_tvd / 2);

	dpfc::histogram_t hist(histTree.begin(), histTree.end());

	for(auto i = 0; i < labels.size(); i++) {
		hist.emplace_back(labels[i], counts[i]);
	}

	std::sort(hist.begin(), hist.end(), [](const dpfc::bin_t& a, const dpfc::bin_t& b) {
			return (a.count == b.count) ? (a.label < b.label) : (a.count > b.count);
		});
	unsigned long b = hist[n].count;

	dpfc::histogram_t rop;
	for(const auto& bin : hist) {
		if (bin.count == b) {
			break;
		}
		rop.emplace_back(bin.label, bin.count);
	}
	return rop;
}

template <class Iterator>
dpfc::histogram_t PureSparseHistogram(const mpz_class& X, Iterator dataset_begin,
		Iterator dataset_end, const mpq_class& beta, const dpfc::CountSample& cs,
		const mpq_class& error_tvd) {
	if (beta <= 0 || beta >= 1) {
		throw std::logic_error("The mixing parameter beta must be in the range (0,1).");
	}

	mpz_class op = 3 * X;
	mpz_pow_ui(op.get_mpz_t(), op.get_mpz_t(), cs.get_n());
	mpq_class delta = beta * (cs.get_eps() / 9) / op;

	dpfc::histogram_t hist = SparseHistogram(X, dataset_begin, dataset_end, delta, 0, cs, error_tvd);

	mpz_class u;
	dpfc::mpz_randint_reject(u.get_mpz_t(), beta.get_den().get_mpz_t());

	if (u <= beta.get_num()) {
		hist = dpfc::RandomHistogram(X, cs.get_n(), error_tvd);
	}
	return hist;
}

#endif


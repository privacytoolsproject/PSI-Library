/*
 * Header file for the algorithms from "Differential Privacy on Finite Computers" Section 8.
 */

#ifndef DPFC_HIST_COMPACT_HPP
#define DPFC_HIST_COMPACT_HPP

#include "dpfc/count.hpp"
#include "dpfc/hist_simple.hpp"
#include "dpfc/utils.hpp"

#include "gmpxx.h"
#include "NTL/GF2E.h"
#include "NTL/GF2EX.h"
#include "NTL/GF2X.h"
#include "NTL/vector.h"

#include <cstdint>
#include <stdexcept>
#include <vector>


namespace dpfc {
	// Convert a byte array to a finite field element. Bits after the first field's degree number of
	// bits are ignored.
	NTL::GF2E GF2E_from_bytes(std::uint8_t const* const data, const size_t nBytes);

	// Convert an integer into a finite field element (bits are interpreted as coefficients).
	NTL::GF2E GF2E_from_mpz(const mpz_class& op);

	// Return a random finite field element.
	// (Note: You must first initialize Libsodium.)
	NTL::GF2E GF2E_random();

	// Convert a finite field element into an integer (coefficients are interpreted as bits).
	mpz_class mpz_from_GF2E(const NTL::GF2E& op);

	// Return the inverse of the function f in the empty-bin sampler (Algorithm 8.5).
	mpz_class f_inv_max(const mpz_class& u, const mpz_class& d0, const mpz_class& d);

	// Return a sample from S_v (defined in Lemma 8.4) where S_v is the uniform distribution over
	// numbers u_0 such that cs.sample(0, f(u_0)) = v.
	// (Note: You must first initialize Libsodium.)
	mpz_class S(unsigned long v, const mpz_class& d0, const dpfc::CountSample& cs);

	// histogram datatype where the histogram is stored as a polynomial over a finite field.
	class compact_histogram_t {
		private:
			mpz_class d0, q, q_a1, r;

		public:
			// the counting query object used to generate noisy counts
			const CountSample* const cs;

			// data universe size
			const mpz_class X;

			// a description of the underlying field
			const NTL::GF2X modulus;

			// the polynomial determining the seeds for the count at each bin
			NTL::GF2EX hash;

			compact_histogram_t(const mpz_class& X, const NTL::GF2X& modulus,
					const CountSample* const cs);

			// Return the noisy count for the given label.
			unsigned long operator [](const mpz_class& label) const;
	};
}

// (Algorithm 8.2) CompactHistogram equipped with the empty-bin sampler Algorithm 8.5 returns a
// (n + 1)-wise independent version of a BasicHistogram represented compactly by the coefficients of
// a polynomial over a finite field. The resulting algorithm is (cs.eps + eps2, 0)-differentially
// private (ignoring privacy loss from error_tvd).
template <class Iterator>
dpfc::compact_histogram_t CompactHistogram(const mpz_class& X, Iterator dataset_begin,
		Iterator dataset_end, const dpfc::CountSample& cs, const mpq_class& eps2,
		const mpq_class& error_tvd);


// --- Implementation ---

template <class Iterator>
dpfc::compact_histogram_t CompactHistogram(const mpz_class& X, Iterator dataset_begin,
		Iterator dataset_end, const dpfc::CountSample& cs, const mpq_class& eps2,
		const mpq_class& error_tvd) {
	if (sodium_init() == -1) {
		throw std::runtime_error("Unable to initialize Libsodium.");
	}
	if (eps2 <= 0) {
		throw std::logic_error("The privacy parameter eps2 must be nonnegative.");
	}
	if (error_tvd <= 0 || error_tvd >= 1) {
		throw std::logic_error("The randomness approximation parameter error_tvd must be in the range"
				"(0,1).");
	}
	if (X <= cs.get_n() || X < 16) {
		throw std::logic_error("The data universe size is too small. Use BasicHistogram instead.");
	}

	auto histTree = dpfc::PartialBasicHistogram(X, dataset_begin, dataset_end, cs);

	// initialize the finite field modulus
	mpq_class d0_q_low = (1 + 1 / eps2) * cs.get_d();
	mpz_class d0_z_low;
	mpz_cdiv_q(d0_z_low.get_mpz_t(), d0_q_low.get_num().get_mpz_t(), d0_q_low.get_den().get_mpz_t());
	if (X > d0_z_low) {
		d0_z_low = X;
	}

	unsigned long ell = dpfc::mpz_clog2(d0_z_low.get_mpz_t());
	d0_z_low = ell;
	ell = dpfc::mpz_clog2(d0_z_low.get_mpz_t()) - 1;

	mpz_class pow;
	unsigned long pow_ui;
	mpz_ui_pow_ui(pow.get_mpz_t(), 3, ell);
	pow_ui = mpz_get_ui(pow.get_mpz_t());

	mpz_class d0;
	mpz_ui_pow_ui(d0.get_mpz_t(), 2, 2 * pow_ui);

	NTL::GF2X modulus;
	NTL::SetCoeff(modulus, 0, 1);
	NTL::SetCoeff(modulus, pow_ui, 1);
	NTL::SetCoeff(modulus, 2 * pow_ui, 1);
	NTL::GF2E::init(modulus);

	NTL::Vec<NTL::GF2E> vecX, vecY;
	vecX.SetLength(cs.get_n() + 1);
	vecY.SetLength(cs.get_n() + 1);

	// pick a random field element that will map to the count given by histTree
	int i = 0;
	mpz_ui_pow_ui(pow.get_mpz_t(), 2, 2 * pow_ui);
	for (const auto& bin : histTree) {
		vecX[i] = dpfc::GF2E_from_mpz(bin.first);
		vecY[i] = dpfc::GF2E_from_mpz(S(bin.second, d0, cs));
		i++;
	}

	// fill up the vector to a total n + 1 bins and set count randomly
	mpz_class label(1);
	for (; i < cs.get_n() + 1; i++) {
		while (histTree.count(label)) {
			label++;
		}
		vecX[i] = dpfc::GF2E_from_mpz(label++);
		vecY[i] = dpfc::GF2E_random();
	}

	dpfc::compact_histogram_t hist(X, modulus, &cs);
	hist.hash = interpolate(vecX, vecY);
	return hist;
}

#endif


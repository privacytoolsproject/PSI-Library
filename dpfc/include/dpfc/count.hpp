/*
 * Header file for the counting query algorithms from "Differential Privacy on Finite Computers"
 * Section 4.
 */

#ifndef DPFC_COUNT_HPP
#define DPFC_COUNT_HPP

#include "gmpxx.h"


namespace dpfc {
	// An abstract class meant to be used for generating noisy counting queries based strictly on
	// using integer arithmetic. This class will handle logic shared by GeoSample and FastSample.
	class CountSample {
		protected:
			// The maximum value achieved by the scaled cdf (i.e. the value returned by
			// CountSample::cdf(center, this->n) should be d for all valid values of center).
			mpz_class d;

			// Privacy parameter epsilon.
			const mpq_class eps;

			// Maximum count value (i.e. the number of rows in the dataset).
			const unsigned long n;

			// Holds the value "2^{ceil(log2(2/eps))}".
			mpz_class pow_2k;

			// Holds the value "2^{ceil(log2(2/eps))} + 1".
			mpz_class pow_2k_a1;

		public:
			// Constructor to initialize n, the privacy parameter and intermediate values. (Keeping
			// with the notation in the paper, passing the privacy parameter eps leads to a single
			// count being eps/2 differentially private.)
			CountSample(const unsigned long n, const mpq_class& eps);

			// Return data member d.
			mpz_class get_d() const;

			// Return data member eps.
			mpq_class get_eps() const;

			// Return data member n.
			unsigned long get_n() const;

			// Perform inverse transform sampling to generate a single sample from the underlying cdf
			// with the given center c. Default randomness is generated using rejection sampling.
			unsigned long sample(const unsigned long c) const;
			unsigned long sample(const unsigned long c, const mpz_class& random) const;

			// Return the smallest b in [0,n] such that Pr[sample(1) > b] <= delta.
			unsigned long threshold(const mpq_class& delta) const;

			// Return an additive bound such that except with probability at most beta a sample will be
			// within this bound of the given center (less than or equal to).
			virtual unsigned long accuracy(const mpq_class& beta) const =0;

			// The underlying cdf distribution centered at c evaluated at z scaled by this->d used to
			// generate samples.
			virtual mpz_class cdf(const unsigned long c, const unsigned long z) const =0;
	};
}

// (Algorithm 4.8) Class to generate random samples using GeoSample.
class GeoSample : public dpfc::CountSample {
	public:
		GeoSample(const unsigned long n, const mpq_class& eps);

		unsigned long accuracy(const mpq_class& beta) const override;
		mpz_class cdf(const unsigned long c, const unsigned long z) const override;
};

// (Algorithm 4.9) Class to generate random samples using FastSample.
class FastSample : public dpfc::CountSample {
	private:
		// Probability of outputting a uniformly random count.
		mpq_class gamma;

		// The maximum value achieved by the scaled cdf before mixing with uniform.
		mpz_class d_prime;

		// The threshold on which to cut off counts from the center.
		unsigned long t;

	public:
		FastSample(const unsigned long n, const mpq_class& eps, const mpq_class& gamma);

		// Return data member gamma.
		mpq_class get_gamma() const;

		unsigned long accuracy(const mpq_class& beta) const override;
		mpz_class cdf(const unsigned long c, const unsigned long z) const override;
};

#endif


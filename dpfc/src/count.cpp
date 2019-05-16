/*
 * Implementation of the counting query algorithms from "Differential Privacy on Finite Computers"
 * Section 4.
 */

#include "dpfc/count.hpp"

#include "dpfc/utils.hpp"

#include "sodium.h"

#include <stdexcept>


using dpfc::CountSample;
using dpfc::mpz_clog2;
using dpfc::mpq_zinv_round;
using std::logic_error;
using std::runtime_error;

namespace dpfc {
	CountSample::CountSample(const unsigned long n, const mpq_class& eps) : n(n), eps(eps) {
		// error handling of parameters
		if (this->n == 0) {
			throw logic_error("n cannot be zero.");
		}

		if (eps <= 0 || eps > 1) {
			throw logic_error("The privacy parameter eps must be in the range (0,1].");
		}

		if (sodium_init() == -1) {
			throw runtime_error("Unable to initialize Libsodium.");
		}

		// eps is being replaced with 1/ceil(1/eps)
		// initialize the intermediary data members holding 2^ceil(log2(2/eps)) and
		// 2^ceil(log2(2/eps)) + 1
		mpq_class eps_zinv;
		mpq_zinv_round(eps_zinv.get_mpq_t(), eps.get_mpq_t());
		unsigned long k = mpz_clog2(eps_zinv.get_den().get_mpz_t());

		mpz_ui_pow_ui(this->pow_2k.get_mpz_t(), 2, k);
		this->pow_2k *= 2;
		this->pow_2k_a1 = this->pow_2k + 1;
	}

	mpz_class CountSample::get_d() const {
		return this->d;
	}

	mpq_class CountSample::get_eps() const {
		return this->eps;
	}

	unsigned long CountSample::get_n() const {
		return this->n;
	}

	unsigned long CountSample::sample(unsigned long c) const {
		mpz_class random;
		mpz_randint_reject(random.get_mpz_t(), this->d.get_mpz_t());

		return this->sample(c, random);
	}

	// Perform inverse transform sampling using binary search on the output space to generate a
	// single sample from the underlying cdf function with the given center center and randomness
	// threshold random. For this function to work properly random must be drawn uniformly from the
	// range [1, ..., this->d].
	unsigned long CountSample::sample(const unsigned long c, const mpz_class& random) const {
		unsigned long low = 0, mid = 0, high = this->n;
		while (low < high) {
			mid = low + (high - low) / 2;
			if (this->cdf(c, mid) >= random) {
				high = mid;
			}
			else {
				low = mid + 1;
			}
		}
		return low;
	}

	// Perform binary search over [0, ..., n] to find the desired threshold value. This value always
	// exists as we have the invariant cdf(1, n) = d.
	unsigned long CountSample::threshold(const mpq_class& delta) const {
		if (delta <= 0 || delta > 1) {
			throw logic_error("The privacy parameter delta must be in the range (0, 1].");
		}

		mpq_class breakpoint = this->d * (1 -  delta);
		unsigned long low = 0, mid = 0, high = this->n;
		while (low < high) {
			mid = low + (high - low) / 2;
			if (this->cdf(1, mid) >= breakpoint) {
				high = mid;
			}
			else {
				low = mid + 1;
			}
		}
		return low;
	}
}

GeoSample::GeoSample(const unsigned long n, const mpq_class& eps) : CountSample(n, eps) {
	// set d to (2*pow_2k + 1) * (pow_2k + 1)^(n-1)
	mpz_pow_ui(this->d.get_mpz_t(), this->pow_2k_a1.get_mpz_t(), n - 1);
	this->d *= 2 * this->pow_2k + 1;
}

// Perform binary search over [0, ..., n] to find an additive bound such that except with
// probability at most beta a sample will be within this bound of the given center.
unsigned long GeoSample::accuracy(const mpq_class& beta) const {
	// returns the smallest a in [0, ..., n] such that
	// 1 - beta <= F(c + a) - F(c - a - 1) which is equivalent to
	// (2^(k+1) + 1) * (2^k + 1)^a * beta >= 2^(k+1) * 2^(ka)

	if (beta <= 0 || beta > 1) {
		throw logic_error("The accuracy parameter beta must be in the range (0, 1].");
	}

	mpq_class lhs_base = beta * (2 * this->pow_2k + 1) / 2 / this->pow_2k;
	mpz_class lhs_var, rhs_var;

	unsigned long low = 0, mid = 0, high = this->n;
	while (low < high) {
		mid = low + (high - low) / 2;

		mpz_pow_ui(lhs_var.get_mpz_t(), this->pow_2k_a1.get_mpz_t(), mid);
		mpz_pow_ui(rhs_var.get_mpz_t(), this->pow_2k.get_mpz_t(), mid);
		if (lhs_base * lhs_var >= rhs_var) {
			high = mid;
		}
		else {
			low = mid + 1;
		}
	}
	return low;
}

mpz_class GeoSample::cdf(unsigned long c, unsigned long z) const {
	if (z >= this->n) {
		return this->d;
	}

	mpz_class rop, temp;
	if (z < c) {
		// set rop to pow_2k^(c-z) * (pow_2k + 1)^(n-(c-z))
		mpz_pow_ui(temp.get_mpz_t(), this->pow_2k.get_mpz_t(), c - z);
		mpz_pow_ui(rop.get_mpz_t(), this->pow_2k_a1.get_mpz_t(), this->n - c + z);
		rop *= temp;
	}
	else {
		// set rop to d - pow_2k^(z-c+1) * (pow_2k + 1)^(n-1-(z-c))
		mpz_pow_ui(temp.get_mpz_t(), this->pow_2k_a1.get_mpz_t(), (this->n - z - 1) + c);
		mpz_pow_ui(rop.get_mpz_t(), this->pow_2k.get_mpz_t(), (z - c) + 1);
		temp *= rop;
		rop = this->d - temp;
	}
	return rop;
}

FastSample::FastSample(const unsigned long n, const mpq_class& eps, const mpq_class& gamma) :
		CountSample(n, eps), gamma(gamma) {
	if (gamma <= 0 || gamma >= 1 || gamma.get_num() != 1) {
		throw logic_error("The parameter gamma must be in the range (0,1) and the reciprocal of "
				"an integer.");
	}

	mpq_class epsRound;
	mpq_zinv_round(epsRound.get_mpq_t(), eps.get_mpq_t());

	// set t to min(n, ceil(9/(2*eps) * ceil(log(8*(n+1)*(1-gamma)/(eps*gamma)))) - 1)
	mpz_class op;
	op = 8 * (this->gamma.get_den() - 1) * (mpz_class(n) + 1) * epsRound.get_den();
	mpz_set_ui(op.get_mpz_t(), mpz_clog2(op.get_mpz_t()));
	op *= 9 * epsRound.get_den();
	mpz_cdiv_q_ui(op.get_mpz_t(), op.get_mpz_t(), 2);

	if (!op.fits_ulong_p() || (this->t = mpz_get_ui(op.get_mpz_t()) - 1) > n) {
		this->t = n;
	}

	// set d' to (2*pow_2k + 1) * (pow_2k + 1)^t
	mpz_pow_ui(this->d_prime.get_mpz_t(), this->pow_2k_a1.get_mpz_t(), this->t);
	this->d_prime *= 2 * this->pow_2k + 1;

	//set d to (n+1) * d' / gamma
	this->d = (mpz_class(n) + 1) * gamma.get_den() * d_prime;
}

mpq_class FastSample::get_gamma() const {
	return this->gamma;
}

// Perform binary search over [0, ..., n] to find an additive bound such that except with
// probability at most beta a sample will be within this bound of the given center.
unsigned long FastSample::accuracy(const mpq_class& beta) const {
	// returns the smallest a in [0, ..., n] such that
	// (2^(k+1) + 1) * (2^k + 1)^a * (beta - gamma) >= 2^(k+1) * 2^(ka)

	if (beta <= this->gamma || beta > 1) {
		throw logic_error("The accuracy parameter beta must be in the range (gamma, 1].");
	}

	mpq_class lhs_base = (beta - this->gamma) * (2 * this->pow_2k + 1) / 2 / this->pow_2k;
	mpz_class lhs_var, rhs_var;

	unsigned long low = 0, mid = 0, high = this->n;
	while (low < high) {
		mid = low + (high - low) / 2;

		mpz_pow_ui(lhs_var.get_mpz_t(), this->pow_2k_a1.get_mpz_t(), mid);
		mpz_pow_ui(rhs_var.get_mpz_t(), this->pow_2k.get_mpz_t(), mid);
		if (lhs_base * lhs_var >= rhs_var) {
			high = mid;
		}
		else {
			low = mid + 1;
		}
	}
	return low;
}

mpz_class FastSample::cdf(unsigned long c, unsigned long z) const {
	if (z >= this->n) {
		return this->d;
	}

	mpz_class rop(0), temp;
	if (z < c) {
		// if z < max{0, c-t} set rop to 0
		// if z in [max{0, c-t}, c) set rop to 2^(k*(c-z)) * (2^k + 1)^(t+1-(c-z)) - 2^(k*(t+1))
		if  (c <= this->t || z >= c - this->t) {
			mpz_pow_ui(rop.get_mpz_t(), this->pow_2k.get_mpz_t(), c - z);
			mpz_pow_ui(temp.get_mpz_t(), this->pow_2k_a1.get_mpz_t(), this->t + 1 - (c - z));
			rop *= temp;

			mpz_pow_ui(temp.get_mpz_t(), this->pow_2k.get_mpz_t(), this->t + 1);
			rop -= temp;
		}
	}
	else if (c >= this->n - this->t || z < c + this->t) {
		// if z in [c, min{c+t, n}) set rop to d' - 2^(k*(z-c+1)) * (2^k + 1)^(t-(z-c)) + 2^(k*(t+1))
		mpz_pow_ui(temp.get_mpz_t(), this->pow_2k.get_mpz_t(), z - c + 1);
		mpz_pow_ui(rop.get_mpz_t(), this->pow_2k_a1.get_mpz_t(), this->t - (z - c));
		temp *= rop;

		rop = this->d_prime - temp;
		mpz_pow_ui(temp.get_mpz_t(), this->pow_2k.get_mpz_t(), this->t + 1);
		rop += temp;
	}
	else {
		// if z in [min{c+t, n}, n) set rop to d'
		rop = this->d_prime;
	}

	// scale the cdf value F'(z) (stored in rop) to account for mixing with uniform
	// (z+1)*d' + (1/gamma - 1)*(n+1)*F'(z)
	temp = this->n;
	temp += 1;
	rop *= temp;

	temp = gamma.get_den();
	temp -= 1;
	rop *= temp;

	temp = z + 1;
	temp *= this->d_prime;
	rop += temp;

	return rop;
}


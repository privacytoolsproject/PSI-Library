/*
 * Miscellaneous functions used in the dpfc project.
 */

#ifndef DPFC_UTILS_HPP
#define DPFC_UTILS_HPP

#include "gmp.h"


namespace dpfc {
	// Set rop to the largest rational less than or equal to op whose reciprocal is an integer
	// (behavior is undefined if op <= 0).
	void mpq_zinv_round(mpq_t rop, const mpq_t op);

	// Return the ceiling of the base 2 logarithm of op (behavior is undefined if op <= 0).
	size_t mpz_clog2(const mpz_t op);

	// Set rop to an integer drawn almost uniformly at random from 1 to op (inclusive). The
	// distribution of rop is guaranteed to be at most error > 0 away from uniform in statistical
	// distance. (Note: You must first initialize Libsodium.)
	void mpz_randint_bias(mpz_t rop, const mpz_t op, const mpq_t error);

	// Set rop to an integer drawn uniformly at random from 1 to op (inclusive) using rejection
	// sampling. When op is a power of 2, this algorithm is guaranteed to only take one sample.
	// (Note: You must first initialize Libsodium.)
	void mpz_randint_reject(mpz_t rop, const mpz_t op);

	// Set rop to op.
	void mpz_set_ull(mpz_t rop, const unsigned long long op);
}

#endif


/*
 * Miscellaneous functions used in the dpfc project.
 */

#include "dpfc/utils.hpp"

#include "sodium.h"

#include <cstdint>


using std::uint8_t;

namespace dpfc {
	void mpq_zinv_round(mpq_t rop, const mpq_t op) {
		mpz_cdiv_q(mpq_denref(rop), mpq_denref(op), mpq_numref(op));
		mpz_set_ui(mpq_numref(rop), 1);
	}

	size_t mpz_clog2(const mpz_t op) {
		const size_t opSize = mpz_sizeinbase(op, 2);
		const size_t opBit1 = mpz_scan1(op, 0);

		return (opSize == opBit1 + 1) ? opBit1 : opSize;
	}

	void mpz_randint_bias(mpz_t rop, const mpz_t op, const mpq_t error) {
		// number of bits to generate: ceil(log_2(ceil(op/error))) - 2
		mpz_mul(rop, op, mpq_denref(error));
		mpz_cdiv_q(rop, rop, mpq_numref(error));

		const size_t bitCount = mpz_clog2(rop);
		const size_t byteCount = (bitCount + 5) / 8;
		uint8_t* random = (uint8_t*)_alloca(byteCount);

		// generate randomness
		randombytes_buf(random, byteCount);
		mpz_import(rop, byteCount, 1, 1, 0, 0, random);

		mpz_tdiv_r(rop, rop, op);
		mpz_add_ui(rop, rop, 1);
	}

	void mpz_randint_reject(mpz_t rop, const mpz_t op) {
		// compute the number of bits to generate
		const size_t bitCount = mpz_clog2(op);
		const size_t byteCount = (bitCount + 7) / 8;
		const size_t bitGenCount = 8 * byteCount;
		uint8_t* random = (uint8_t*)_alloca(byteCount);

		// repeat until we get a good sample
		do {
			// generate randomness
			randombytes_buf(random, byteCount);
			mpz_import(rop, byteCount, 1, 1, 0, 0, random);

			// truncate extra bits
			for (int i = bitCount ; i < bitGenCount; i++) {
				mpz_clrbit(rop, i);
			}
		} while (mpz_cmp(rop, op) >= 0);

		mpz_add_ui(rop, rop, 1);
	}

	void mpz_set_ull(mpz_t rop, const unsigned long long op) {
		mpz_import(rop, 1, 1, sizeof(op), 0, 0, &op);
	}
}


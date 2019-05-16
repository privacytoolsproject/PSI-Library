/*
 * Implementation of algorithms from "Differential Privacy on Finite Computers" Section 8.
 */

#include "dpfc/hist_compact.hpp"

#include <algorithm>
#include <sstream>
#include <string>


using NTL::deg;
using NTL::GF2;
using NTL::GF2E;
using NTL::GF2EX;
using NTL::GF2X;
using NTL::SetCoeff;
using std::hex;
using std::logic_error;
using std::reverse;
using std::string;
using std::stringstream;
using std::uint8_t;
using std::vector;

namespace dpfc {
	GF2E GF2E_from_bytes(uint8_t const* const data, const size_t nBytes) {
		size_t byteCap = (GF2E::degree() + 7) / 8;
		size_t minCap = (nBytes < byteCap ? nBytes : byteCap - 1);

		stringstream ss;
		ss << "0x" << hex;
		size_t i = 0;
		for (; i < minCap; i++) {
			ss << (data[i] & 15) << ((data[i] >> 4) & 15);
		}

		if (i == byteCap - 1 && i < nBytes) {
			uint8_t b = data[i] % (1 << (GF2E::degree() % 8));
			ss << (b & 15) << ((b >> 4) & 15);
		}

		GF2E rop;
		ss >> rop;
		return rop;
	}

	GF2E GF2E_from_mpz(const mpz_class& op) {
		stringstream ss;
		ss << "0x" << hex << (op - 1);
		string s = ss.str();
		reverse(s.begin()+2, s.end());
		ss.str(s);
		ss.clear();

		GF2E rop;
		ss >> rop;
		return rop;
	}

	GF2E GF2E_random() {
		size_t nbytes = (GF2E::degree() + 7) / 8;
		uint8_t* rand = (uint8_t*)_alloca(nbytes);
		randombytes_buf(rand, nbytes);

		return GF2E_from_bytes(rand, nbytes);
	}

	mpz_class mpz_from_GF2E(const GF2E& op) {
		stringstream ss;
		int temp = GF2X::HexOutput;
		GF2X::HexOutput = 1;
		ss << op;
		GF2X::HexOutput = temp;

		string s = ss.str();
		s.erase(0, 2);
		reverse(s.begin(), s.end());
		ss.str(s);
		ss.clear();

		mpz_class rop;
		ss >> hex >> rop;
		return ++rop;
	}

	mpz_class f_inv_max(const mpz_class& u, const mpz_class& d0, const mpz_class& d) {
		mpz_class r = d0 % d;
		return (d0 / d) * u + (((r != 0) && (u < r)) ? u : r);
	}

	mpz_class S(unsigned long v, const mpz_class& d0, const CountSample& cs) {
		mpz_class FMin(0);
		if (v > 0) {
			FMin = cs.cdf(0, v - 1);
		}
		mpz_class RMin = f_inv_max(FMin, d0, cs.get_d());
		mpz_class RMax = f_inv_max(cs.cdf(0, v), d0, cs.get_d()) - RMin;

		mpz_class rop;
		mpz_randint_reject(rop.get_mpz_t(), RMax.get_mpz_t());
		rop += RMin;
		return rop;
	}

	compact_histogram_t::compact_histogram_t(const mpz_class& X, const GF2X& modulus,
			const CountSample* const cs) : X(X), modulus(modulus), cs(cs) {
		mpz_ui_pow_ui(this->d0.get_mpz_t(), 2, deg(modulus));

		if (this->d0 < X) {
			throw logic_error("The cardinality of the field must be at least the data universe size.");
		}
		if (this->d0 < cs->get_d()) {
			throw logic_error("The cardinality of the field is too small for the given CountSample.");
		}

		this->q = d0 / cs->get_d();
		this->q_a1 = this->q + 1;
		this->r = d0 % cs->get_d();
	}

	unsigned long compact_histogram_t::operator [](const mpz_class& label) const {
		if ((label < 1) || (label > this->X)) {
			throw logic_error("The given label does not correspond to an element of the data "
					"universe.");
		}
		if (this->modulus != GF2E::modulus()) {
			GF2E::init(this->modulus);
		}

		NTL::GF2E field = GF2E_from_mpz(label);
		NTL::GF2E u0;
		NTL::eval(u0, this->hash, field);

		// set u to f(u0)
		mpz_class u = mpz_from_GF2E(u0);
		if ((this->r != 0) && (u <= this->r * (this->q + 1))) {
			mpz_cdiv_q(u.get_mpz_t(), u.get_mpz_t(), this->q_a1.get_mpz_t());
		}
		else {
			u -= this->r;
			mpz_cdiv_q(u.get_mpz_t(), u.get_mpz_t(), this->q.get_mpz_t());
		}
		return this->cs->sample(0, u);
	}
}


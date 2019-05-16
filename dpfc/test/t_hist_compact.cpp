
#include "dpfc/hist_compact.hpp"

#include "gtest/gtest.h"

#include <sstream>
#include <string>


using dpfc::compact_histogram_t;
using dpfc::CountSample;
using dpfc::f_inv_max;
using dpfc::GF2E_from_bytes;
using dpfc::GF2E_from_mpz;
using dpfc::GF2E_random;
using dpfc::mpz_from_GF2E;
using NTL::deg;
using NTL::GF2E;
using NTL::GF2X;
using NTL::SetCoeff;
using std::logic_error;
using std::runtime_error;
using std::string;
using std::stringstream;
using std::uint8_t;
using std::vector;

TEST(GF2E_from_bytes, Standard) {
	GF2X modulus;
	SetCoeff(modulus, 0, 1);
	SetCoeff(modulus, 20, 1);
	SetCoeff(modulus, 75, 1);
	GF2E::init(modulus);

	uint8_t bytes[10] = {1, 128, 255, 20, 5, 247, 0, 0, 202, 255};

	stringstream ss;
	ss << GF2E_from_bytes(bytes, 9);
	EXPECT_EQ(ss.str(), "[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 1 0 1 0 0 0 "
								"1 0 1 0 0 0 0 0 1 1 1 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
								"0 1 0 1 0 0 1 1]");
	ss.str("");
	ss.clear();

	ss << GF2E_from_bytes(bytes, 10);
	EXPECT_EQ(ss.str(), "[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 1 0 1 0 0 0 "
								"1 0 1 0 0 0 0 0 1 1 1 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
								"0 1 0 1 0 0 1 1 1 1 1]");
}

TEST(GF2E_from_bytes, Small) {
	GF2X modulus;
	SetCoeff(modulus, 5, 1);
	GF2E::init(modulus);

	uint8_t byte = 132;
	stringstream ss;
	ss << GF2E_from_bytes(&byte, 1);
	EXPECT_EQ(ss.str(), "[0 0 1]");
}

TEST(GF2E_from_mpz, Standard) {
	GF2X modulus;
	SetCoeff(modulus, 0, 1);
	SetCoeff(modulus, 12, 1);
	SetCoeff(modulus, 42, 1);
	GF2E::init(modulus);

	EXPECT_EQ(GF2E_from_mpz(mpz_class(1)), GF2E::zero());

	stringstream ss;
	ss << GF2E_from_mpz(mpz_class("3401302760720"));
	EXPECT_EQ(ss.str(), "[1 1 1 1 0 0 0 0 1 0 1 1 1 0 1 0 1 0 0 0 1 1 1 0 1 0 1 1 0 1 1 1 "
								"1 1 1 0 1 0 0 0 1 1]");
}

TEST(GF2E_random, Nonzero) {
	if (sodium_init() == -1) {
		throw runtime_error("Unable to initialize Libsodium.");
	}

	GF2X modulus;
	SetCoeff(modulus, 0, 1);
	SetCoeff(modulus, 1000, 1);
	GF2E::init(modulus);

	EXPECT_NE(GF2E_random(), GF2E::zero()); // whp
}

TEST(mpz_from_GF2E, Standard) {
	GF2X modulus;
	SetCoeff(modulus, 0, 1);
	SetCoeff(modulus, 52, 1);
	SetCoeff(modulus, 73, 1);
	GF2E::init(modulus);

	EXPECT_EQ(mpz_from_GF2E(GF2E::zero()), 1);

	mpz_class num("291150252670408097792");
	EXPECT_EQ(mpz_from_GF2E(GF2E_from_mpz(num)), num);

	num = "9444732965739290427392";
	EXPECT_EQ(mpz_from_GF2E(GF2E_from_mpz(num)), num);
}

TEST(f_inv_max, Standard) {
	mpz_class d0(5000), d(250);
	EXPECT_EQ(f_inv_max(1, d0, d), 20);
	EXPECT_EQ(f_inv_max(149, d0, d), 2980);
	EXPECT_EQ(f_inv_max(250, d0, d), 5000);

	d = 247;
	EXPECT_EQ(f_inv_max(1, d0, d), 21);
	EXPECT_EQ(f_inv_max(60, d0, d), 1260);
	EXPECT_EQ(f_inv_max(61, d0, d), 1280);
	EXPECT_EQ(f_inv_max(247, d0, d), 5000);
}

unsigned long Mf(const mpz_class& u0, const mpz_class& d0, const CountSample& cs) {
	mpz_class q(d0 / cs.get_d()), r(d0 % cs.get_d()), u;

	if ((r != 0) && (u0 <= r * (q + 1))) {
		q++;
		mpz_cdiv_q(u.get_mpz_t(), u0.get_mpz_t(), q.get_mpz_t());
	}
	else {
		u = u0 - r;
		mpz_cdiv_q(u.get_mpz_t(), u.get_mpz_t(), q.get_mpz_t());
	}
	return cs.sample(0, u);
}

TEST(S, Standard) {
	if (sodium_init() == -1) {
		throw runtime_error("Unable to initialize Libsodium.");
	}

	mpz_class d0;
	mpz_ui_pow_ui(d0.get_mpz_t(), 2 , 486);

	GeoSample gs(100, mpq_class("1/2"));
	mpz_class u = S(35, d0, gs);
	EXPECT_EQ(Mf(S(0, d0, gs), d0, gs), 0);
	EXPECT_EQ(Mf(S(35, d0, gs), d0, gs), 35);
	EXPECT_EQ(Mf(S(78, d0, gs), d0, gs), 78);
	EXPECT_EQ(Mf(S(100, d0, gs), d0, gs), 100);
}

TEST(compact_histogram_t, Init) {
	mpq_class eps("1/2");
	mpz_class X;
	mpz_ui_pow_ui(X.get_mpz_t(), 2, 100);
	GeoSample gs(10, eps);
	GF2X modulus;
	SetCoeff(modulus, 100, 1);
	EXPECT_THROW(compact_histogram_t(X + 1, modulus, &gs), logic_error);

	GeoSample gs2(100, eps);
	EXPECT_THROW(compact_histogram_t(1000, modulus, &gs2), logic_error);

	compact_histogram_t hist(X, modulus, &gs);
	EXPECT_EQ(hist.X, X);
	EXPECT_EQ(hist.modulus, modulus);
	EXPECT_EQ(hist.cs->get_n(), 10);
	EXPECT_EQ(hist.cs->get_eps(), eps);
}

// testing counting query class that returns the randomness minus 1 as the noisy count
class TestSample : public CountSample {
	public:
		TestSample(const unsigned long X) : CountSample(X - 1, 1) {
			this->d = X;
		}

		mpz_class cdf(const unsigned long c, const unsigned long z) const {
			return z + 1;
		}

		unsigned long accuracy(const mpq_class& beta) const {
			return 0;
		}
};

TEST(compact_histogram_t, Subscript) {
	GF2X modulus;
	SetCoeff(modulus, 0, 1);
	SetCoeff(modulus, 3, 1);
	SetCoeff(modulus, 6, 1);
	GF2E::init(modulus);
	TestSample ts(31);

	compact_histogram_t hist(31, modulus, &ts);

	uint8_t bytes[2] = {20, 3};
	SetCoeff(hist.hash, 0, GF2E_from_bytes(bytes, 1));
	SetCoeff(hist.hash, 1, GF2E_from_bytes(bytes + 1, 1));

	EXPECT_THROW(hist[0], logic_error);
	EXPECT_THROW(hist[40], logic_error);

	EXPECT_EQ(hist[1], 9);
	EXPECT_EQ(hist[20], 15);

	SetCoeff(modulus, 6, 0);
	GF2E::init(modulus);

	EXPECT_EQ(hist[5], 11);
}

TEST(CompactHistogram, Error) {
	mpz_class X("500");
	vector<mpz_class> ds(99, 1);
	mpq_class eps("1/2"), error_tvd("1/100");
	GeoSample gs(100, eps);

	EXPECT_THROW(CompactHistogram(X, ds.cbegin(), ds.cend(), gs, eps, error_tvd), logic_error);

	ds.push_back(2);
	EXPECT_THROW(CompactHistogram(X, ds.cbegin(), ds.cend(), gs, 0, error_tvd), logic_error);
	EXPECT_THROW(CompactHistogram(X, ds.cbegin(), ds.cend(), gs, eps, 0), logic_error);
	EXPECT_THROW(CompactHistogram(X, ds.cbegin(), ds.cend(), gs, eps, 1), logic_error);

	EXPECT_THROW(CompactHistogram(100, ds.cbegin(), ds.cend(), gs, eps, error_tvd), logic_error);
}

TEST(CompactHistogram, Standard) {
	mpz_class X("500");
	vector<mpz_class> ds(200, 1);
	mpq_class eps("1/2"), error_tvd("1/100");
	GeoSample gs(200, eps);

	auto hist = CompactHistogram(X, ds.cbegin(), ds.cend(), gs, eps / 2, error_tvd);

	EXPECT_NE(hist.hash, GF2E::zero()); // whp

	mpz_class d0;
	mpz_ui_pow_ui(d0.get_mpz_t(), 2, deg(hist.modulus));

	EXPECT_GE(d0, X);
	EXPECT_GE(d0, mpq_class("5/4") * gs.get_d());

	EXPECT_EQ(hist.X, X);
	EXPECT_GE(hist[1], 100); // whp
	EXPECT_LE(hist[2], 100); // whp
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


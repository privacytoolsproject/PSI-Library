
#include "dpfc/utils.hpp"

#include "gmpxx.h"
#include "gtest/gtest.h"
#include "sodium.h"

#include <stdexcept>


using dpfc::mpq_zinv_round;
using dpfc::mpz_clog2;
using dpfc::mpz_randint_bias;
using dpfc::mpz_randint_reject;
using dpfc::mpz_set_ull;
using std::runtime_error;

TEST(mpq_zinv_round, Standard) {
	mpq_class op(13, 50), eop(1, 4), rop(0,1);
	mpq_zinv_round(rop.get_mpq_t(), op.get_mpq_t());
	EXPECT_EQ(rop, eop);

	op = "3/20";
	eop = "1/7";
	mpq_zinv_round(rop.get_mpq_t(), op.get_mpq_t());
	EXPECT_EQ(rop, eop);

	op = "343/15";
	eop = "1/1";
	mpq_zinv_round(rop.get_mpq_t(), op.get_mpq_t());
	EXPECT_EQ(rop, eop);
}

TEST(mpq_zinv_round, InPlace) {
	mpq_class op(13, 50), eop(1, 4);
	mpq_zinv_round(op.get_mpq_t(), op.get_mpq_t());
	EXPECT_EQ(op, eop);
}

TEST(mpz_clog2, Standard) {
	mpz_class op(1000);
	size_t rop = mpz_clog2(op.get_mpz_t());
	EXPECT_EQ(rop, 10);

	op = 1;
	rop = mpz_clog2(op.get_mpz_t());
	EXPECT_EQ(rop, 0);

	op = 2048;
	rop = mpz_clog2(op.get_mpz_t());
	EXPECT_EQ(rop, 11);

	// op is set to "2^(2^31) + 1"
	mpz_ui_pow_ui(op.get_mpz_t(), 2, 2147483648);
	op += 1;
	rop = mpz_clog2(op.get_mpz_t());
	EXPECT_EQ(rop, 2147483649);
}

TEST(mpz_randint_bias, Range) {
	if (sodium_init() == -1) {
		throw std::runtime_error("Unable to initialize Libsodium.");
	}

	mpq_class error("7/1000000000000");
	mpz_class op(1), rop(0);
	mpz_randint_bias(rop.get_mpz_t(), op.get_mpz_t(), error.get_mpq_t());
	EXPECT_EQ(rop, 1);

	op = 1;
	rop = 0;
	mpz_setbit(op.get_mpz_t(), 100);
	mpz_randint_bias(rop.get_mpz_t(), op.get_mpz_t(), error.get_mpq_t());
	EXPECT_GE(rop, 1);
	EXPECT_LE(rop, op);
}

TEST(mpz_randint_reject, Range) {
	if (sodium_init() == -1) {
		throw runtime_error("Unable to initialize Libsodium.");
	}

	mpz_class op(1), rop(0);
	mpz_randint_reject(rop.get_mpz_t(), op.get_mpz_t());
	EXPECT_EQ(rop, 1);

	op = 1;
	rop = 0;
	mpz_setbit(op.get_mpz_t(), 100);
	mpz_randint_reject(rop.get_mpz_t(), op.get_mpz_t());
	EXPECT_GE(rop, 1);
	EXPECT_LE(rop, op);
}

TEST(mpz_set_ull, Standard) {
	mpz_class rop, eop("1099511627777");
	unsigned long long op = 1099511627777;

	mpz_set_ull(rop.get_mpz_t(), op);
	EXPECT_NE(rop, 1);
	EXPECT_EQ(rop, eop);
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


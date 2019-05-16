
#include "dpfc/hist_sparse.hpp"

#include "gtest/gtest.h"


using dpfc::ApproxBinSample;
using dpfc::ApproxConvExp;
using dpfc::DistinctSample;
using dpfc::histogram_map_t;
using dpfc::histogram_t;
using dpfc::mpz_randint_bias;
using dpfc::mpz_randint_reject;
using dpfc::RandomHistogram;
using NTL::conv;
using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZX;
using std::logic_error;
using std::runtime_error;
using std::vector;

mpz_class bigLabel("19204812490834545"), sX("400"), bX("292502840823049823432442");
vector<mpz_class> ds_sX, ds_bX;
unsigned long sn = 200, bn = 63453;
GeoSample scs(sn, mpq_class("1/10"));
FastSample bcs(bn, mpq_class("1/100"), mpq_class(1, 8 * bX));

TEST(ApproxConvExp, Exact) {
	// compute exactly the distribution of flipping 30 fair coins
	ZZ p = conv<ZZ>("536870912");
	ZZ s = conv<ZZ>("1073741824");
	ZZX a;
	SetCoeff(a, 0, s - p);
	SetCoeff(a, 1, p);

	ZZX rop = ApproxConvExp(a, 30, s, 31);

	EXPECT_EQ(deg(rop), 30);
	for (int k = 0; k <= 30; k++) {
		ZZ binom(1);
		for (int i = 0; i < k; i++) {
			binom *= (30 - i);
			binom /= (i + 1);
		}
		EXPECT_EQ(rop[k], binom);
	}

	// approximation does not change as t varies
	ZZX ropTrunc = ApproxConvExp(a, 30, s, 15);

	EXPECT_EQ(deg(ropTrunc), 14);
	for (int i = 0; i <= 14; i++) {
		EXPECT_EQ(ropTrunc[i], rop[i]);
	}
}

TEST(ApproxConvExp, Standard) {
	// a case with rounding
	ZZ s = conv<ZZ>("20");
	ZZX a, rop;
	SetCoeff(a, 0, s - 7);
	SetCoeff(a, 1, 7);

	rop = ApproxConvExp(a, 6, s, 10);

	EXPECT_EQ(deg(rop), 3);
	EXPECT_EQ(rop[0], 1);
	EXPECT_EQ(rop[1], 4);
	EXPECT_EQ(rop[2], 5);
	EXPECT_EQ(rop[3], 3);
}

void checkRange(const vector<unsigned long> res, unsigned long len, unsigned long top) {
	EXPECT_EQ(res.size(), len);
	unsigned long p = top;
	for (auto count : res) {
		EXPECT_LE(count, p);
		p = count;
	}
}

TEST(ApproxOrdSample, Range) {
	if (sodium_init() == -1) {
		throw runtime_error("Unable to initialize Libsodium.");
	}

	unsigned long n = 100;
	FastSample cs(n, mpq_class("1/10"), mpq_class("1/13123193"));
	mpz_class m("1231028"), s("2934823904");
	mpq_class error_tvd("1/1000000000000");

	auto res = ApproxOrdSample(m, s, cs, error_tvd);
	checkRange(res, n + 1, n);

	res = ApproxOrdSample(m, s, cs, error_tvd);
	checkRange(res, n + 1, n);
}

TEST(ApproxOrdSampleClamp, Range) {
	if (sodium_init() == -1) {
		throw runtime_error("Unable to initialize Libsodium.");
	}

	unsigned long n = 1000, k = 20;
	GeoSample cs(n, mpq_class("1/2"));
	mpz_class m("43083404"), s("7326814253");
	mpq_class error_tvd("1/1000000000000");

	auto res = ApproxOrdSampleClamp(m, s, cs, k, error_tvd);
	checkRange(res, n + 1, k);

	k = 24;
	res = ApproxOrdSampleClamp(m, s, cs, k, error_tvd);
	checkRange(res, n + 1, k);
}


TEST(DistinctSample, Standard) {
	if (sodium_init() == -1) {
		throw runtime_error("Unable to initialize Libsodium.");
	}

	histogram_map_t A;
	for (int i = 3; i <= 100; i += 3) {
		A[i] = 1;
	}

	auto S = DistinctSample(100, A, 67, mpq_class(1, 3000));
	EXPECT_EQ(S.size(), 67);

	// S should contain all number less than or equal to 100 not divisible by 3
	for (int i = 1; i <= 100; i++) {
		if (i % 3 != 0) {
			EXPECT_TRUE(find(S.begin(), S.end(), i) != S.end());
		}
	}

	// larger random test
	A.empty();
	mpz_class sample, m(500000);
	for (int i = 0; i < 20000; i++) {
		mpz_randint_bias(sample.get_mpz_t(), m.get_mpz_t(), mpq_class(1, 10000000).get_mpq_t());
		A[sample] = 1;
	}

	// sort results to check distinctness and check that each is not in A
	S = DistinctSample(m, A, 10000, mpq_class(1, 3000));
	EXPECT_EQ(S.size(), 10000);

	sort(S.begin(), S.end());
	EXPECT_EQ(A.count(S[0]), 0);
	for (size_t i = 1; i < 10000; i++) {
		EXPECT_GT(S[i], S[i - 1]);
		EXPECT_EQ(A.count(S[i]), 0);
	}
}

// ensure histogram is sorted by count then by label
void isSortedCountLabel(const histogram_t& res, unsigned long n, const mpz_class& X) {
	EXPECT_LE(res.size(), sn);

	mpz_class pc = n, pl = 0;
	for (auto bin : res) {
		EXPECT_LE(bin.count, pc);
		if (bin.count == pc) {
			EXPECT_GT(bin.label, pl);
		}
		pc = bin.count;
		pl = bin.label;
	}
	EXPECT_GT(pc, 0);
	EXPECT_LE(pl, X);
}

TEST(RandomHistogram, Range) {
	mpq_class error_tvd("1/1000000000");
	mpz_class X("10000000");
	unsigned long n = 1000;

	for (int k = 0; k < 2; k++) {
		auto hist = RandomHistogram(X, n, error_tvd);
		isSortedCountLabel(hist, n, X);
	}
}

TEST(SparseHistogram, Error) {
	mpq_class delta(1, bX), zeta(1, bX), error_tvd("1/12301239473290384234");

	EXPECT_THROW(SparseHistogram(bX, ds_bX.cbegin(), ds_bX.cend(),
				0, zeta, bcs, error_tvd), logic_error);
	EXPECT_THROW(SparseHistogram(bX, ds_bX.cbegin(), ds_bX.cend(),
				delta, zeta, bcs, 1), logic_error);
	EXPECT_THROW(SparseHistogram(bX, ds_bX.cbegin(), ds_bX.cend(),
				delta, -1, bcs, error_tvd), logic_error);
	EXPECT_THROW(SparseHistogram(bX, ds_bX.cbegin(), ds_bX.cend(),
				4, zeta, bcs, error_tvd), logic_error);
	EXPECT_THROW(SparseHistogram(bX, ds_bX.cbegin(), ds_bX.cend(),
				delta, zeta, bcs, 0), logic_error);
	EXPECT_THROW(SparseHistogram(sX, ds_bX.cbegin(), ds_bX.cend(),
				delta, zeta, bcs, error_tvd), logic_error);
	EXPECT_THROW(SparseHistogram(sX, ds_sX.cbegin(), ds_sX.cend(),
				delta, zeta, scs, error_tvd), logic_error);
}

TEST(SparseHistogram, Standard) {
	mpq_class delta(1, bn * bX), zeta(1, bn * bX), error_tvd("1/8976576675645465475685");

	auto res = SparseHistogram(bX, ds_sX.cbegin(), ds_sX.cend(), delta, zeta, scs, error_tvd);
	isSortedCountLabel(res, sn, bX);
}

TEST(PurseSparseHistogram, Error) {
	mpq_class beta("1/50"), error_tvd("1/12301239473290384234");

	EXPECT_THROW(PureSparseHistogram(bX, ds_bX.cbegin(), ds_bX.cend(),
				0, bcs, error_tvd), logic_error);
	EXPECT_THROW(PureSparseHistogram(bX, ds_bX.cbegin(), ds_bX.cend(),
				2, bcs, error_tvd), logic_error);
	EXPECT_THROW(PureSparseHistogram(bX, ds_bX.cbegin(), ds_bX.cend(),
				beta, bcs, 5), logic_error);
}

TEST(PurseSparseHistogram, Standard) {
	mpq_class beta("1/100"), error_tvd("1/75456486565831256");

	auto res = PureSparseHistogram(bX, ds_sX.cbegin(), ds_sX.cend(), beta, scs, error_tvd);
	isSortedCountLabel(res, sn, bX);
}

int main(int argc, char **argv) {
	for(int i = 0; i < sn; i++) {
		ds_sX.push_back(1 + (i % 20));
	}

	for(int i = 1; i <= bn; i++) {
		ds_bX.push_back(100 * i);
	}

	random_shuffle(ds_sX.begin(), ds_sX.end());
	random_shuffle(ds_bX.begin(), ds_bX.end());

	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


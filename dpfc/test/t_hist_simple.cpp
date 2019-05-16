
#include "dpfc/hist_simple.hpp"

#include "gtest/gtest.h"


using dpfc::bin_t;
using dpfc::CountSample;
using std::logic_error;
using std::runtime_error;
using std::vector;

mpz_class bigLabel("23492340983248549"), sX("10"), bX("100000000000000000000000000000");
vector<mpz_class> ds_sX(5000), ds_bX(45000);

// if the true count is 1, then the "noisy" count is b and otherwise the noisy count is one more
// than the true count (with wraparound for n)
class TestSample : public CountSample {
	public:
		unsigned long b = 2;

		TestSample(const unsigned long n, const mpq_class& eps) : CountSample(n, eps) {
			this->d = 1;
		}

		mpz_class cdf(const unsigned long c, const unsigned long z) const override {
			if (c == 1) {
				return (z >= this->b);
			}
			else if (c == n) {
				return 1;
			}
			return z >= c + 1;
		}

		unsigned long accuracy(const mpq_class& beta) const override {
			return 1;
		}
};

TEST(PartialBasicHistogram, Error) {
	GeoSample gs(ds_bX.size(), mpq_class(1, 2));

	// Invalid universe label
	EXPECT_THROW(PartialBasicHistogram(0, ds_bX.cbegin(), ds_bX.cend(), gs), logic_error);

	// size mismatch
	EXPECT_THROW(PartialBasicHistogram(sX, ds_sX.cbegin(), ds_sX.cend(), gs), logic_error);
}

TEST(PartialBasicHistogram, Standard) {
	mpq_class eps(1, 1);
	TestSample ts_s(ds_sX.size(), eps);

	auto hist = PartialBasicHistogram(sX, ds_sX.cbegin(), ds_sX.cend(), ts_s);

	// check it has entries for bins 1 to 9 and not 10
	for(int i = 1; i < 10; i++) {
		EXPECT_EQ(hist.count(i), 1);
	}
	EXPECT_EQ(hist.count(10), 0);

	// test on larger dataset
	TestSample ts_b(ds_bX.size(), eps);
	hist = PartialBasicHistogram(bX, ds_bX.cbegin(), ds_bX.cend(), ts_b);

	EXPECT_EQ(hist.size(), 20001);
	for(int i = 1; i <= 20000; i++) {
		EXPECT_EQ(hist.count(5 * i), 1);
		EXPECT_EQ(hist[5 * i], 3);
	}
	EXPECT_EQ(hist.count(bigLabel), 1);
	EXPECT_EQ(hist[bigLabel], 5001);
}

TEST(BasicHistogram, Error) {
	GeoSample gs(ds_bX.size(), mpq_class(1, 2));
	EXPECT_THROW(BasicHistogram(bX, ds_bX.cbegin(), ds_bX.cend(), gs), runtime_error);
}

TEST(BasicHistogram, Standard) {
	TestSample ts(ds_sX.size(), mpq_class(1, 1));
	auto hist = BasicHistogram(sX, ds_sX.cbegin(), ds_sX.cend(), ts);
	for(int i = 1; i <= 10; i++) {
		EXPECT_EQ(hist[i - 1].label, i);
		EXPECT_EQ(hist[i - 1].count, (i <= 5) ? 557 : ((i < 10) ? 556 : 1));
	}
}

TEST(StabilityHistogram, Error) {
	GeoSample gs(ds_bX.size(), mpq_class(1, 2));
	EXPECT_THROW(StabilityHistogram(bX, ds_bX.cbegin(), ds_bX.cend(), gs, mpq_class(0, 1)),
			logic_error);
	EXPECT_THROW(StabilityHistogram(bX, ds_bX.cbegin(), ds_bX.cend(), gs, mpq_class(1, 1)),
			logic_error);
}

TEST(StabilityHistogram, Standard) {
	mpq_class delta(1, 2);

	// threshold at 2
	TestSample ts(ds_bX.size(), mpq_class(1, 2));
	auto hist = StabilityHistogram(bX, ds_bX.cbegin(), ds_bX.cend(), ts, delta);
	EXPECT_EQ(hist.size(), 20001);
	for(int i = 1; i <= 20000; i++) {
		EXPECT_EQ(hist[i - 1].label, 5 * i);
		EXPECT_EQ(hist[i - 1].count, 3);
	}
	EXPECT_EQ(hist[20000].label, bigLabel);
	EXPECT_EQ(hist[20000].count, 5001);

	// threshold at 3
	ts.b = 3;
	hist = StabilityHistogram(bX, ds_bX.cbegin(), ds_bX.cend(), ts, delta);
	EXPECT_EQ(hist.size(), 1);
	EXPECT_EQ(hist[0].label, bigLabel);
	EXPECT_EQ(hist[0].count, 5001);

}

int main(int argc, char **argv) {
	for(int i = 0; i < 5000; i++) {
		ds_sX.push_back(1 + (i % 9));
		ds_bX.push_back(bigLabel);
	}

	for(int i = 1; i <= 20000; i++) {
		ds_bX.push_back(5 * i);
		ds_bX.push_back(5 * i);
	}

	random_shuffle(ds_sX.begin(), ds_sX.end());
	random_shuffle(ds_bX.begin(), ds_bX.end());

	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


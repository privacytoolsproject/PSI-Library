/*
 * Implementation of helper algorithms from "Differential Privacy on Finite Computers" Section 6.
 */

#include "dpfc/hist_sparse.hpp"

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <functional>


using __gnu_pbds::null_type;
using __gnu_pbds::rb_tree_tag;
using __gnu_pbds::tree;
using __gnu_pbds::tree_order_statistics_node_update;
using dpfc::CountSample;
using dpfc::mpz_randint_bias;
using dpfc::mpz_set_ull;
using NTL::coeff;
using NTL::conv;
using NTL::div;
using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZX;
using std::less;
using std::logic_error;
using std::vector;

namespace dpfc {
	ZZX ApproxConvExp(const ZZX& a, const mpz_class& i, const ZZ& s, unsigned long t) {
		ZZX rop;
		if(i == 0) {
			SetCoeff(rop, 0, s);
			return rop;
		}
		else if(i == 1) {
			return a;
		}

		rop = a;
		ZZ s2 = s * s;

		// iterate over i from most to least significant bit after skipping the first one
		size_t k = mpz_sizeinbase(i.get_mpz_t(), 2) - 1;
		do {
			k--;
			SqrTrunc(rop, rop, t);
			// multiply by an extra a/s factor if the bit is set
			if (mpz_tstbit(i.get_mpz_t(), k)) {
				MulTrunc(rop, rop, a, t);

				// set each coefficient to its floor divided by s^2
				for(long c = 0; c <= deg(rop); c++) {
					div(rop[c], rop[c], s2);
				}
				rop.normalize();
			}
			else {
				// set each coefficient to its floor divided by s
				for(long c = 0; c <= deg(rop); c++) {
					div(rop[c], rop[c], s);
				}
				rop.normalize();
			}
		} while (k > 0);

		return rop;
	}

	unsigned long ApproxBinSample(const mpz_class& m, const mpq_class& p, const mpz_class& s,
			const unsigned long t, const mpz_class& u) {
		if (t == 0) {
			return 0;
		}

		mpq_class sp = s * p;
		mpz_class p_prime = sp.get_num() / sp.get_den();

		ZZ p_zz = conv<ZZ>(p_prime.get_str().c_str());
		ZZ s_zz = conv<ZZ>(s.get_str().c_str());
		ZZ u_zz = conv<ZZ>(u.get_str().c_str());

		ZZX a, d;
		SetCoeff(a, 0, s_zz - p_zz);
		SetCoeff(a, 1, p_zz);

		for (unsigned long k = 1; k < 2 * t; k *= 2) {
			d = ApproxConvExp(a, m, s_zz, k);

			for (long ell = 1; ell < k; ell++) {
				SetCoeff(d, ell, coeff(d, ell) + coeff(d, ell - 1));
			}

			for (long ell = k / 2; ell < k; ell++) {
				if (coeff(d, ell) >= u_zz) {
					return (ell < t) ? ell : t;
				}
			}
		}
		return t;
	}

	vector<unsigned long> ApproxOrdSample(const mpz_class& m, const mpz_class& s,
			const CountSample& cs, const mpq_class& error_tvd) {
		return ApproxOrdSampleClamp(m, s, cs, cs.get_n(), error_tvd);
	}

	vector<unsigned long> ApproxOrdSampleClamp(const mpz_class& m, const mpz_class& s,
			const CountSample& cs, unsigned long k, const mpq_class& error_tvd) {
		unsigned long assigned = 0, n = cs.get_n();
		vector<unsigned long> rop;
		rop.reserve(n + 1);
		mpz_class u;
		mpq_class round_error = error_tvd / k;
		mpq_class p;
		unsigned long t = n + 1, ell;

		for (unsigned long v = k; v > 0; v--) {
			// generate the number of counts taking on the value v
			mpz_randint_bias(u.get_mpz_t(), s.get_mpz_t(), round_error.get_mpq_t());

			if ((v == k) && (k < n)) {
				p = mpq_class(cs.get_d() - cs.cdf(0, k - 1), cs.get_d());
			}
			else {
				p = mpq_class(cs.cdf(0, v) - cs.cdf(0, v - 1), cs.cdf(0, v));
			}

			ell = ApproxBinSample(m - assigned, p, s, t, u);

			// add that number of values to the return vector
			for(unsigned long i = 0; i < ell; i++) {
				rop.push_back(v);
			}
			assigned += ell;
			t -= ell;
		}

		// add zero counts
		for (unsigned long i = 0; i < t; i++) {
			rop.push_back(0);
		}
		return rop;
	}

	vector<mpz_class> DistinctSample(const mpz_class& X, const histogram_map_t& A, unsigned long r,
			const mpq_class& error_tvd) {
		// create a order statistic tree with the labels in A
		tree<mpz_class, null_type, less<mpz_class>, rb_tree_tag, tree_order_statistics_node_update> T;
		mpz_class remaining(X);
		for(const auto x : A) {
			T.insert(x.first);
			remaining -= 1;
		}

		vector<mpz_class> S;
		mpq_class round_error = error_tvd / r;
		mpz_class low, mid, high, z, breakpoint;

		for (; r > 0; r--, remaining--) {
			mpz_randint_bias(z.get_mpz_t(), remaining.get_mpz_t(), round_error.get_mpq_t());

			// map the uniform random to an unpicked label
			low = 1;
			high = X;
			while (low < high) {
				mid = (low + high) / 2;
				mpz_set_ull(breakpoint.get_mpz_t(), T.order_of_key(mid));
				breakpoint += z + (T.find(mid) != T.end());
				if (mid >= breakpoint) {
					high = mid;
				}
				else {
					low = mid + 1;
				}
			}

			// add the new label to the tree and output
			T.insert(low);
			S.push_back(low);
		}

		return S;
	}

	histogram_t RandomHistogram(const mpz_class& X, unsigned long n, const mpq_class& error_tvd) {
		histogram_map_t hist;
		mpq_class round_error = error_tvd / n / 2;
		mpz_class sample, count, n_mpz(n);

		for (unsigned long i = 0; i < n; i++) {
			mpz_randint_bias(sample.get_mpz_t(), X.get_mpz_t(), round_error.get_mpq_t());
			if (hist.count(sample) == 0) {
				mpz_randint_bias(count.get_mpz_t(), n_mpz.get_mpz_t(), round_error.get_mpq_t());
				hist[sample] = count.get_ui();
			}
		}

		// sort the resulting histogram to use the same sort order as SparseHistogram
		histogram_t hist_flat(hist.begin(), hist.end());
		sort(hist_flat.begin(), hist_flat.end(), [](const dpfc::bin_t& a, const dpfc::bin_t& b) {
			return (a.count == b.count) ? (a.label < b.label) : (a.count > b.count);
		});

		dpfc::histogram_t hist_ret;
		for(const auto& bin : hist_ret) {
			hist_ret.emplace_back(std::pair<mpz_class, unsigned long>(bin.label, bin.count));
		}
		return hist_ret;
	}
}


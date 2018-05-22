#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <functional>
#include <vector>
#include <cfenv>
#include <limits>

#include "biquintic.hpp"
#include "elements.hpp"
#include "physcon.hpp"
#include "real.hpp"
#include "fermi_dirac.hpp"

struct thermodynamic_data_t {
	real p;
	real dp_drho;
	real dp_dT;
	real e;
	real de_drho;
	real de_dT;
	real ne;
};

void electron_eta(real ne, real T, real& eta, real& deta_dne, real& deta_dT) {
	static const real c0 = two * std::pow(two * M_PI * me * kb / (h * h), 1.5);
	eta = -std::log(c0) - three / two * std::log(T) + std::log(ne);
	deta_dne = one / ne;
	deta_dT = -three / two / T;
}

thermodynamic_data_t fast_saha2(const std::vector<real>& fracs, real rho,
		real T) {
	thermodynamic_data_t results;

	real n = zero;
	real abar = zero;
	std::vector<real> N(NELE);
	std::vector<std::vector<real> > W(NELE);

	// Compute coefficient matrix
	for (int i = 0; i < NELE; i++) {
		const auto& ele = elements[i + 1];
		const real this_rho = fracs[i] * rho;
		N[i] = this_rho / (amu * ele.A);
		n += N[i];
	}
	abar = rho / (n * amu);

	// Find maximum ne possible
	real ne_max = real(0);
	for (int i = 0; i < NELE; i++) {
		ne_max += N[i] * real(i + 1);
	}

	real eta_0, dummy;
	electron_eta(ne_max * 1.0e-7, T, eta_0, dummy, dummy);
	// Compute coefficient matrix
	for (int i = 0; i < NELE; i++) {
		const auto& ele = elements[i + 1];
		W[i].resize(i + 2);
		W[i][0] = real(0);
		real max_w = real(0);
		for (int j = 0; j < i + 1; j++) {
			W[i][j + 1] = W[i][j] + ele.log_saha(j, T) - eta_0;
			max_w = std::max(max_w, W[i][j + 1]);
		}
		for (int j = 0; j <= i + 1; j++) {
			W[i][j] = std::exp(W[i][j] - max_w);
		}

	}

	real ne = ne_max / real(2);

	real f, df_dne;
	real s0, s1, s2;

	// N-R solve for ne
	int iters = 0;
	real ne_last;
	do {
		f = ne;
		df_dne = real(1);
		real eta, deta_dne, deta_dT;
		electron_eta(ne, T, eta, deta_dne, deta_dT);
		const real exp_neta = std::exp(-eta + eta_0);
		for (int i = 0; i < NELE; i++) {
			s1 = s2 = real(0);
			s0 = W[i][0];
			real exp_n_neta = one;
			for (int j = 1; j <= i + 1; j++) {
				exp_n_neta *= exp_neta;
				const auto tmp = W[i][j] * exp_n_neta;
				s0 += tmp;
				s1 += real(j) * tmp;
				s2 += real(j * j) * tmp;
			}

			f -= N[i] * s1 / s0;
			df_dne += (s2 / s0 - (s1 / s0) * (s1 / s0)) * N[i] * deta_dne;
		}
		printf("%14e %14e %14e %14e %14e %14e %14e\n", real(ne), real(f),
				real(df_dne), real(ne / ne_max), s0, s1, s2);
		ne_last = ne;
		ne = std::min(std::max(ne / real(100), ne - f / df_dne),
				(ne_max * 0.99 + 0.01 * ne));
		iters++;
		if (iters > 400) {
			abort();
		}
	} while (std::abs(ne - ne_last) / ne_max > 1.0e-14);

	real eta, deta_dne, deta_dT;
	electron_eta(ne, T, eta, deta_dne, deta_dT);
	real dne_dn = one;
	real dne_dT = zero;
	const real exp_neta = std::exp(-eta + eta_0);
	real e0, e1, e_i, e2;
	e_i = zero;
	real de_dne, de_dT;
	de_dne = zero;
	de_dT = zero;
	for (int i = 0; i < NELE; i++) {
		s1 = s2 = real(0);
		e0 = e1 = e2 = zero;
		s0 = W[i][0];
		real exp_n_neta = one;
		real ei_sum = zero;
		for (int j = 1; j <= i + 1; j++) {
			exp_n_neta *= exp_neta;
			const auto tmp = W[i][j] * exp_n_neta;
			ei_sum += elements[i + 1].e_i[j];

			e1 += real(j) * tmp * ei_sum;
			s0 += tmp;

			e0 += tmp * ei_sum;
			s1 += real(j) * tmp;

			e2 += tmp * ei_sum * ei_sum;
			s2 += real(j * j) * tmp;
		}
		e1 /= s0;
		e0 /= s0;
		s1 /= s0;
		s2 /= s0;
		e2 /= s0;
		e_i += e0 * fracs[i];
		de_dne -= deta_dne * fracs[i] * (e1 - e0 * s1);
		de_dT -= deta_dT * fracs[i] * (e1 - e0 * s1);
		de_dT += fracs[i] * (e2 - e0 * e0) / (kb * T * T);
		dne_dn += (s2 - s1 * s1) * N[i] * deta_dne;
		real tmp1 = (s2 - s1 * s1) * N[i] * deta_dT;
		real tmp2 = (e1 - s1 * e0) * N[i] / (kb * T * T);
		printf("%e %e %e %e %e %e %e\n", tmp1, tmp2, e0, e1, s0, s1, s2);
		dne_dT -= tmp1 - tmp2;
	}
	de_dT *= n / rho;
	dne_dn = ne / (n * dne_dn);
	dne_dT *= (n / ne) * dne_dn;
	de_dT *= (n / ne) * dne_dn;

	printf("%e\n", (n / ne) * dne_dn);
	printf("%e %e %e %e %e %e\n", n, ne, T, dne_dn, dne_dT, e_i);
	results.ne = ne;

	results.p = kb * (n + ne) * T;
	results.e = (three / two * results.p / rho) + e_i / (abar * amu);
	results.dp_dT = kb * (n + ne + dne_dT * T);
	results.de_dT = three / two * results.dp_dT / rho + de_dT;
	const real dp_dn = kb * T + kb * dne_dn * T;
	results.dp_drho = dp_dn / (amu * abar);
	results.de_drho = de_dne * dne_dn / (amu * amu * abar * abar);
	results.de_drho += (three / two) * kb * T * (dne_dn - ne / n)
			/ (rho * abar * amu);
	return results;
}


int main() {
	feenableexcept(FE_OVERFLOW);
	feenableexcept(FE_DIVBYZERO);

	return 0;


}


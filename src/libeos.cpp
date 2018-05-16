#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <functional>
#include <vector>
#include <cfenv>
#include <limits>

#include "bicubic.hpp"
#include "elements.hpp"
#include "physcon.hpp"
#include "real.hpp"
#include "fermi_dirac.hpp"

real N_ele(real eta, real beta) {
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real f1 = FermiDirac(one / two, eta, beta);
	const real f2 = FermiDirac(three / two, eta, beta);
	return c0 * std::pow(beta, 1.5) * (f1 + beta * f2);
}

real dN_ele_deta(real eta, real beta) {
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real df1_deta = dFermiDirac_deta(one / two, eta, beta);
	const real df2_deta = dFermiDirac_deta(three / two, eta, beta);
	return c0 * std::pow(beta, 1.5) * (df1_deta + beta * df2_deta);
}

real N_pos(real eta, real beta) {
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real eta2 = -eta - two / beta;
	if (eta2 > -200.0) {
		const real f1 = FermiDirac(one / two, eta2, beta);
		const real f2 = FermiDirac(three / two, eta2, beta);
		return c0 * std::pow(beta, 1.5) * (f1 + beta * f2);
	} else {
		return 0.0;
	}
}

real dN_pos_deta(real eta, real beta) {
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real eta2 = -eta - two / beta;
	if (eta2 > -200.0) {
		const real df1_deta = -dFermiDirac_deta(one / two, eta2, beta);
		const real df2_deta = -dFermiDirac_deta(three / two, eta2, beta);
		return c0 * std::pow(beta, 1.5) * (df1_deta + beta * df2_deta);
	} else {
		return 0.0;
	}
}


real electron_chemical_potential(real ne, real beta) {
	const real T = beta * me * c * c / kb;
	const auto F = [ne,beta]( real eta ) {
		return ne - (N_ele(eta,beta)-N_pos(eta,beta));
	};

	const auto dF_deta = [ne,beta]( real eta ) {
		return -dN_ele_deta(eta,beta)+dN_pos_deta(eta,beta);
	};

	const real lambda = std::sqrt(h * h / (two * M_PI * me * kb * T));
	const real eta0 = std::log(half * std::pow(lambda, 3) * ne);
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three) * std::pow(beta, 1.5);
	const real eta1 = std::pow(ne * 1.5 / c0, two / three);
	const real eta2 = std::pow(18.0, 1.0 / 6.0) * std::pow(ne, 1.0 / 3.0)
			* std::pow(beta, -.5) * std::pow(c0, -1.0 / 3.0);
	real eta_guess;
	real min_eta = -1.0 / beta;
	if (eta0 < zero) {
		eta_guess = eta0;
	} else {
		double l = 1.5;
		eta_guess = std::pow(std::pow(eta1, -l) + std::pow(eta2, -l), -1.0 / l);
	}
	real eta = eta_guess;
	real f, deta, err;
	int dir = 0;
	int last_dir = 0;
	double w = 1.0;
	do {
		last_dir = dir;
		f = F(eta);
		real df_deta = dF_deta(eta);
		deta = -f / df_deta;
		if (eta == 0.0) {
			err = 1.0;
		} else {
			err = std::abs(deta / eta);
		}
		dir = deta > zero ? +1 : -1;
		if (dir * last_dir == -1) {
			w *= 0.99;
		}
		real eta0 = eta + w * deta;
		eta = std::max(eta0, (eta + min_eta) / 2.0);
	} while (err > 1.0e-12 && std::abs(deta) > 1.0e-14);
//	printf( "%e %e %e\n", ne, T, eta);
	return eta;
}

real electron_free_energy( real ne, real T ) {
	const real beta =  kb * T / me / c / c;
	const real eta = electron_chemical_potential(ne,beta);


		static const real c0 = real(64) * std::atan(1) * std::sqrt(2)
				* std::pow(me * c / h, three) / real(3) * me * c * c;
		const real eta2 = -eta - two / beta;
		real f3, f4;
		const real f1 = FermiDirac(three / two, eta, beta);
		const real f2 = FermiDirac(five / two, eta, beta);
		if (eta2 > -200.0) {
			f3 = FermiDirac(three / two, eta2, beta);
			f4 = FermiDirac(five / two, eta2, beta);
		} else {
			f3 = f4 = 0.0;
		}
		const real pele = c0 * std::pow(beta, five / two) * (f1 + beta / two * f2);
		const real ppos = c0 * std::pow(beta, five / two) * (f3 + beta / two * f4);
		const real pressure = pele + ppos;
		const real free_energy = eta * kb * ne * T - pressure;

}

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
	const real huge = 1e+100;
	const real large = 1e+20;
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

real fast_saha(const std::vector<real>& n, std::vector<std::vector<real>>& ni,
		real T) {
	static const real c0 = two * std::pow(two * M_PI * me * kb / (h * h), 1.5);
	constexpr
	real saha_zero = real(0);
	constexpr
	real saha_one = real(1);

	constexpr
	real huge = 1.0e+200;
	constexpr
	real large = 1.0e+20;
	constexpr
	real small = 1.0e-20;
	constexpr
	real tiny = 1.0e-200;

	std::vector<std::vector<real>> W;
	ni.resize(NSPECIES);
	W.resize(NSPECIES);
	for (int i = 0; i < NSPECIES; i++) {
		ni[i].resize(i + 1);
		W[i].resize(i + 1);
	}
	real& ne = ni[0][0];
	real ne_max = saha_zero;
	for (int i = 1; i < NSPECIES; i++) {
		ne_max += n[i] * i;
	}
	ne = ne_max * half;
	real ne_last, ne_next;
	real K = c0 * std::pow(T, 1.5);
	real v;
	for (int j = 1; j < NSPECIES; j++) {
		W[j][0] = v = saha_one;
		for (int i = 0; i < j; i++) {
			v *= elements[j].saha(i, T);
			if (v < huge) {
				NULL;
			} else {
				for (int k = 0; k < i; k++) {
					W[j][k] = saha_zero;
				}
				v = saha_one;
			}
			W[j][i + 1] = v;
		}
	}
	std::vector<real> ne_inv_pwr(NSPECIES);
	do {
		ne_last = ne;
		real nsum = saha_zero;
		real ne_inv = saha_one / ne;
		const real kone = K * ne_inv;

		ne_inv_pwr[0] = saha_one;
		for (int i = 1; i < NSPECIES; i++) {
			ne_inv_pwr[i] = kone * ne_inv_pwr[i - 1];
			if (ne_inv_pwr[i] < huge) {
				NULL;
			} else {
				ne_inv_pwr[i] = huge;
			}
		}
		for (int j = 1; j < NSPECIES; j++) {
#pragma ivdep
			for (int i = 0; i <= j; i++) {
				ni[j][i] = W[j][i] * ne_inv_pwr[i];
			}
			nsum = ni[j][0];
			for (int i = 1; i <= j; i++) {
				nsum += ni[j][i];
			}
			const real factor = n[j] / nsum;
			for (int i = 0; i <= j; i++) {
				ni[j][i] *= factor;
			}
		}
		ne_next = saha_zero;
		for (int j = 1; j < NSPECIES; j++) {
			for (int i = 1; i <= j; i++) {
				ne_next += i * ni[j][i];
			}
		}
		printf("%e %e\n", ne, ne_last);
		const real w = real(1) / real(10);
		ne = ne_next * (saha_one - w) + w * ne;
	} while (std::abs(ne - ne_last) / ne_max > 1.0e-10);
	return ne / ne_max;
}

#include <fenv.h>

int main() {
	feenableexcept(FE_OVERFLOW);
	feenableexcept(FE_DIVBYZERO);

	std::vector<real> fracs(NSPECIES, 0.0);
	const real n = 1.0;
	const real np = n * 1.00001;
	const real T = 1.0;
	const real Tp1 = T * 1.00001;
	fracs[0] = 1.0;
	auto th1 = fast_saha2(fracs, n * amu * 1.0079, T);
	auto th2 = fast_saha2(fracs, np * amu * 1.0079, T);
	auto th3 = fast_saha2(fracs, n * amu * 1.0079, Tp1);
	printf("%e %e %e\n", th1.p, th1.dp_drho,
			(th1.p - th2.p) / (n - np) / (amu * 1.0079));
	printf("%e %e %e\n", th1.p, th1.dp_dT, (th1.p - th3.p) / (T - Tp1));

	printf("%e %e %e %e\n", th1.e, th2.e, (th1.de_drho +  th2.de_drho)/2.0* n * (amu * 1.0079),
			(th1.e - th2.e) / (n - np) / (amu * 1.0079) * n * (amu * 1.0079));
	printf("%e %e %e\n", th1.e, th1.de_dT, (th1.e - th3.e) / (T - Tp1));
	return 0;
	fracs[0] = 0.0;
	fracs[1] = 1.0;

	std::vector<std::vector<real>> Ni;
	std::vector<real> N(NSPECIES);
	for (int i = 0; i < NSPECIES; i++) {
		N[i] = n * fracs[i];
	}
	real frac = fast_saha(N, Ni, real(T));
	fprintf(stdout, "%e %e %e\n", n, T, double(frac));

	return 0;

	static constexpr real n_min = 1.0e-4;
	static constexpr real n_max = 1.0e+35;
	static constexpr real T_min = 1.0e+1;
	static constexpr real T_max = 1.0e+12;
	int count = 0;
	for (real n = n_min; n < n_max; n *= 2) {
		for (real T = T_min; T < T_max; T *= 2) {
			std::vector<real> N(NSPECIES);
			std::vector<std::vector<real>> Ni;
			count++;
			for (int i = 0; i < NSPECIES; i++) {
				N[i] = n * fracs[i];
			}
			real frac = fast_saha(N, Ni, real(T));
			fprintf(stdout, "%e %e %e\n", n, T, double(frac));
		}
	}
	fprintf(stderr, "%i\n", count);
	//	printf( "%e\n", electron_chemical_potential(1.000000e+35, kb * 1.000000e+01/me/c/c));
//	electron_table();
	return 0;
}


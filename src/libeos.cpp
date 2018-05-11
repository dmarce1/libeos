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

void electron_state(real eta, real beta, real& pressure, real& energy) {
	static const real c0 = real(64) * std::atan(1) * std::sqrt(2)
			* std::pow(me * c / h, three) / real(3) * me * c * c;
	const real eta2 = -eta - two / beta;

	real pele, ppos, eele, epos;
	real f3, f4;
	const real f1 = FermiDirac(three / two, eta, beta);
	const real f2 = FermiDirac(five / two, eta, beta);
	if (eta2 > -200.0) {
		f3 = FermiDirac(three / two, eta2, beta);
		f4 = FermiDirac(five / two, eta2, beta);
	} else {
		f3 = f4 = 0.0;
	}

	pele = c0 * std::pow(beta, five / two) * (f1 + beta / two * f2);
	ppos = c0 * std::pow(beta, five / two) * (f3 + beta / two * f4);
	eele = 1.5 * c0 * std::pow(beta, five / two) * (f1 + beta * f2);
	epos = 1.5 * c0 * std::pow(beta, five / two) * (f3 + beta * f4);
	pressure = pele + ppos;
	energy = eele + epos;
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

std::vector<real> saha_ratios(int Z, real ne, real T) {
	const real mu1 = electron_chemical_potential(ne, kb * T / me / c / c);
	static const real c0 = two * std::pow(two * M_PI * me * kb / (h * h), 1.5);
	std::vector<real> n(Z + 1);
	std::vector<real> r(Z);
	bool all_le1 = true;
	bool all_ge1 = true;
	const real c1 = c0 * std::pow(T, 1.5) / ne;
	real mu2 = -std::log(c1);
	printf("%e %e %e %e \n", ne, T, mu1, mu2);
	for (int i = 0; i < Z; i++) {
		r[i] = std::exp(-mu1) * elements[Z].saha(i, T);
		all_le1 = all_le1 && (r[i] <= one);
		all_ge1 = all_le1 && (r[i] >= one);
	}
	if (all_le1) {
		n[0] = one;
		for (int i = 0; i < Z; i++) {
			n[i + 1] = r[i] * n[i];
		}
	} else if (all_ge1) {
		n[Z] = one;
		for (int i = Z; i > 0; i--) {
			n[i - 1] = n[i] / r[i - 1];
		}
	} else {
		n[0] = one;
		for (int i = 0; i < Z; i++) {
			n[i + 1] = r[i] * n[i];
			if (n[i + 1] > one) {
				for (int j = 0; j <= i; j++) {
					n[j] /= n[i + 1];
				}
				n[i + 1] = one;
			}
		}
	}
	real nsum = zero;
	for (int i = 0; i <= Z; i++) {
		nsum += n[i];
	}
	for (int i = 0; i <= Z; i++) {
		n[i] /= nsum;
	}
	return n;
}

class electron_table {
private:
	static constexpr real ne_min = 1.0e-4;
	static constexpr real ne_max = 1.0e+35;
	static constexpr real T_min = 1.0e+1;
	static constexpr real T_max = 1.0e+12;
	std::shared_ptr<bicubic_table> mu_table;
	std::shared_ptr<bicubic_table> p_table;
public:
	electron_table() {
		/*
		 mu_table = std::make_shared<bicubic_table>([](real log_ne, real log_T) {
		 const real T = std::exp(log_T);
		 const real ne = std::exp(log_ne);
		 const real beta = kb * T / me / c / c;
		 return electron_chemical_potential(ne,beta);
		 }, std::log(ne_min), std::log(ne_max), std::log(T_min), std::log(T_max),
		 1.0e-2);

		 */
		p_table = std::make_shared<bicubic_table>([](real log_ne, real log_T) {
			const real T = std::exp(log_T);
			const real ne = std::exp(log_ne);
			const real beta = kb * T / me / c / c;
			const real eta= electron_chemical_potential(ne,beta);
			real p, e;
			electron_state(eta,beta,p,e);
			return std::log(p);
		}, std::log(ne_min), std::log(ne_max), std::log(T_min), std::log(T_max),
				1.0e-3);

		FILE* fp = fopen("electron.dat", "wt");
		for (real ne = ne_min; ne < ne_max; ne *= 2.0) {
			for (real T = T_min; T < T_max; T *= 2.0) {
				real p, e, eta;
				eta = electron_chemical_potential(ne, kb * T / me / c / c);
				const real beta = kb * T / me / c / c;
				electron_state(eta, beta, p, e);
				fprintf(fp, "%e %e %e %e\n", ne, T,
						std::exp((*p_table)(std::log(ne), std::log(T))), p);
			}
		}
		fclose(fp);
	}
};

class saha_table {
private:
	int Z;
	std::shared_ptr<bicubic_table> Z_table;
	std::shared_ptr<bicubic_table> e_table;
	static constexpr real ne_min = 1.0e-4;
	static constexpr real ne_max = 1.0e+35;
	static constexpr real T_min = 1.0e+1;
	static constexpr real T_max = 1.0e+12;
public:
	saha_table(int z) {
		Z = z;
		std::function<real(real, real)> Z_func([this](real x, real y) {
			const real ne = std::exp(x);
			const real T = std::exp(y);
			auto n = saha_ratios(Z,ne,T);
			real zbar = zero;
			for( int i = 1; i <=Z; i++ ) {
				zbar += n[i] * real(i);
			}
			return zbar;
		});
		std::function<real(real, real)> e_func([this](real x, real y) {
			const real ne = std::exp(x);
			const real T = std::exp(y);
			auto n = saha_ratios(Z,ne,T);
			real etot = zero;
			for( int i = 0; i <=Z; i++ ) {
				etot += n[i] * elements[Z].e_i[i];
			}
			return etot;
		});
		Z_table = std::make_shared<bicubic_table>(Z_func, std::log(ne_min),
				std::log(ne_max), std::log(T_min), std::log(T_max), 1.0e-6);
		e_table = std::make_shared<bicubic_table>(e_func, std::log(ne_min),
				std::log(ne_max), std::log(T_min), std::log(T_max), 1.0e-6);
		FILE* fp = fopen("eos.dat", "wt");
		for (real ne = ne_min; ne < ne_max; ne *= 1.1) {
			for (real T = T_min; T < T_max; T *= 1.1) {
				fprintf(fp, "%e %e %e %e\n", ne, T,
						(*Z_table)(std::log(ne), std::log(T)),
						(*e_table)(std::log(ne), std::log(T)));
			}
		}
		fclose(fp);
	}
	void save(FILE* fp) const {
		Z_table->save(fp);
	}
	void load(FILE* fp) {
		Z_table = std::make_shared<bicubic_table>(fp);
	}

};

using saha_real = double;
saha_real fast_saha(const std::vector<saha_real>& n,
		std::vector<std::vector<saha_real>>& ni, saha_real T) {
	static const real c0 = two * std::pow(two * M_PI * me * kb / (h * h), 1.5);
	constexpr saha_real saha_zero = saha_real(0);
	constexpr saha_real saha_one = saha_real(1);

	constexpr saha_real huge = 1.0e+200;
	constexpr saha_real large = 1.0e+20;
	constexpr saha_real small = 1.0e-20;
	constexpr saha_real tiny = 1.0e-200;

	std::vector<std::vector<saha_real>> W;
	ni.resize(NSPECIES);
	W.resize(NSPECIES);
	for (int i = 0; i < NSPECIES; i++) {
		ni[i].resize(i + 1);
		W[i].resize(i + 1);
	}
	saha_real& ne = ni[0][0];
	saha_real ne_max = saha_zero;
	for (int i = 1; i < NSPECIES; i++) {
		ne_max += n[i] * i;
	}
	ne = ne_max * half;
	saha_real ne_last, ne_next;
	saha_real K = c0 * std::pow(T, 1.5);
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
	std::vector<saha_real> ne_inv_pwr(NSPECIES);
	do {
		ne_last = ne;
		saha_real nsum = saha_zero;
		saha_real ne_inv = saha_one / ne;
		const saha_real kone = K * ne_inv;

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
			const saha_real factor = n[j] / nsum;
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
		const saha_real w = saha_real(1) / saha_real(10);
		ne = ne_next * (saha_one - w) + w * ne;
	} while (std::abs(ne - ne_last) / ne_max > 1.0e-12);
	return ne / ne_max;
}

#include <fenv.h>

int main() {
	feenableexcept(FE_OVERFLOW);
	feenableexcept(FE_DIVBYZERO);

	std::vector<saha_real> fracs(NSPECIES, 0.0);

	fracs[2] = 1.0;

	static constexpr real n_min = 1.0e-4;
	static constexpr real n_max = 1.0e+35;
	static constexpr real T_min = 1.0e+1;
	static constexpr real T_max = 1.0e+12;
	int count = 0;
	for (real n = n_min; n < n_max; n *= 1.1) {
		for (real T = T_min; T < T_max; T *= 1.1) {
			std::vector<saha_real> N(NSPECIES);
			std::vector<std::vector<saha_real>> Ni;
			count++;
			for (int i = 0; i < NSPECIES; i++) {
				N[i] = n * fracs[i];
			}
			saha_real frac = fast_saha(N, Ni, saha_real(T));
			fprintf( stdout, "%e %e %e\n", n, T, double(frac));
		}
	}
	fprintf( stderr, "%i\n", count);
	//	printf( "%e\n", electron_chemical_potential(1.000000e+35, kb * 1.000000e+01/me/c/c));
//	electron_table();
	return 0;
}


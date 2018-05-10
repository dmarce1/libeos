#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <vector>

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
	const real f1 = FermiDirac(one / two, eta2, beta);
	const real f2 = FermiDirac(three / two, eta2, beta);
	return c0 * std::pow(beta, 1.5) * (f1 + beta * f2);
}

real dN_pos_deta(real eta, real beta) {
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real eta2 = -eta - two / beta;
	const real df1_deta = -dFermiDirac_deta(one / two, eta2, beta);
	const real df2_deta = -dFermiDirac_deta(three / two, eta2, beta);
	return c0 * std::pow(beta, 1.5) * (df1_deta + beta * df2_deta);
}

real electron_state(real eta, real beta, real& pressure, real& energy) {
	static const real c0 = real(64) * std::atan(1) * std::sqrt(2)
			* std::pow(me * c / h, three) / real(3) * me * c * c;
	const real eta2 = -eta - two / beta;

	real pele, ppos, eele, epos;

	const real f1 = FermiDirac(three / two, eta, beta);
	const real f2 = FermiDirac(five / two, eta, beta);
	const real f3 = FermiDirac(three / two, eta2, beta);
	const real f4 = FermiDirac(five / two, eta2, beta);

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
	if (eta0 < zero) {
		eta_guess = eta0;
	} else {
		double l = 1.5;
		eta_guess = std::pow(std::pow(eta1, -l) + std::pow(eta2, -l), -1.0 / l);
	}
	real eta = eta_guess;
	real f, deta, err;
	do {
		f = F(eta);
		real df_deta = dF_deta(eta);
		deta = -f / df_deta;
		if (eta == 0.0) {
			err = 1.0;
		} else {
			err = std::abs(deta / eta);
		}
		eta += deta;
	} while (err > 1.0e-12 && std::abs(deta) > 1.0e-14);
	return eta;
}

std::vector<real> saha_ratios(int Z, real ne, real T) {
	static const real c0 = two * std::pow(two * M_PI * me * kb / (h * h), 1.5);
	std::vector<real> n(Z + 1);
	std::vector<real> r(Z);
	bool all_le1 = true;
	bool all_ge1 = true;
	const real c1 = c0 * std::pow(T, 1.5) / ne;
	for (int i = 0; i < Z; i++) {
		r[i] = c1 * elements[Z].saha(i, T);
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

struct single_state {
	int Z;
	real ne;
	real ni;
	real T;
	real p;
	real e;
	real eta;
	single_state(int _Z, real _ne, real _T, bool fully_ionized = false) {
		Z = _Z;
		ne = _ne;
		T = _T;
		real p_ele;
		real e_ele;
		real p_ion;
		real e_ion;
		real e_exc;
		const real beta = kb * T / (me * c * c);
		eta = electron_chemical_potential(ne, beta);
		electron_state(eta, beta, p_ele, e_ele);
		auto n = saha_ratios(Z, ne, T);
		real ne_per_n = zero;
		e_exc = zero;
		if (!fully_ionized) {
			for (int i = 0; i < Z; i++) {
				ne_per_n += real(i) * n[i];
				e_exc += elements[Z].e_i[i] * n[i];
			}
			ni = ne / ne_per_n;
			e_exc *= ni;
		} else {
			ni = ne / Z;
		}
		p_ion = kb * ni * T;
		e_ion = 1.5 * p_ion;
		e = e_ele + e_ion + e_exc;
		p = p_ele + p_ion;
	}
};

class ionization_fraction_table {
private:
	int Z;
	const std::function<real(real, real)> gen_func;
	const bicubic_table table;
	static constexpr real ne_min = 1.0;
	static constexpr real ne_max = 1.0e+40;
	static constexpr real T_min = 1.0e+0;
	static constexpr real T_max = 1.0e+12;
	static constexpr int NT = 100;
	static constexpr int NE = 100;
public:
	ionization_fraction_table(int z) :
			Z(z), gen_func([this](real x, real y) {
				const real ne = std::exp(x);
				const real T = std::exp(y);

				auto n = saha_ratios(Z,ne,T);
				real zbar = zero;
				for( int i = 1; i <=Z; i++ ) {
					zbar += n[i] * real(i);
				}
				const real ion_frac = zbar / Z;
				const real rc = ion_frac;;
		//		printf( "%e %e %e\n", ne, T, ion_frac);
				return rc;
			}), table(gen_func, std::log(ne_min), std::log(ne_max), std::log(T_min), std::log(T_max), 5.0e-2) {
		FILE* fp = fopen( "eos.dat", "wt" );
		for( real T = T_min; T <= T_max; T *= 1.1 ) {
			for( real ne = ne_min; ne < ne_max; ne *= 1.1 ) {
				fprintf( fp, "%e %e %e\n", ne, T, table(std::log(ne),std::log(T)));
			}
		}
		fclose(fp);

	}
};

int main() {
	ionization_fraction_table(28);
	return 0;
}


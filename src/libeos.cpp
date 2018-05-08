#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <vector>

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

real electron_pressure(real eta, real beta) {
	static const real c0 = real(64) * std::atan(1) * std::sqrt(2)
			* std::pow(me * c / h, three) / real(3) * me * c * c;
	const real eta2 = -eta - two / beta;

	real pele, ppos;

	const real f1 = FermiDirac(three / two, eta, beta);
	const real f2 = FermiDirac(five / two, eta, beta);
	const real f3 = FermiDirac(three / two, eta2, beta);
	const real f4 = FermiDirac(five / two, eta2, beta);

	pele = c0 * std::pow(beta, five / two) * (f1 + beta / two * f2);
	ppos = c0 * std::pow(beta, five / two) * (f3 + beta / two * f4);
	return pele + ppos;
}

real electron_chemical_potential(real ne, real T) {
	const real beta = (kb * T) / (me * c * c);
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
	const real eta2 = std::pow(18.0,1.0/6.0) * std::pow(ne,1.0/3.0) * std::pow(beta,-.5) * std::pow(c0,-1.0/3.0);
	real eta_guess;
	if (eta0 < zero) {
		eta_guess = eta0;
	} else {
		double l = 1.5;
		eta_guess = std::pow(std::pow(eta1,-l) + std::pow(eta2,-l),-1.0/l);
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
	} while (err > 1.0e-12);
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
}

void partial_state(int Z, real ne, real T, real& rho, real& ene, real& P) {
	static const real c0 = two * std::pow(two * M_PI * me * kb / (h * h), 1.5);
	real eta = electron_chemical_potential(ne, T);
	real ne_spec = zero, nsum = one, n = one;
	ene = zero;
	const real c1 = c0 * std::pow(T, 1.5) / ne;
	for (int i = 0; i < Z; i++) {
		n *= c1 * elements[Z].saha(i, T);
		nsum += n;
		ne_spec += real(i + 1) * n;
		ene += elements[Z].e_i[i + 1] * n;
	}
	ne_spec /= nsum;
	ene /= nsum;
	n = ne / ne_spec;
	ene *= n;
	P = kb * n * T;
	P += electron_pressure(eta, T);
	ene += P * three / two + two * me * c * c * N_pos(eta, T);
	rho = n * elements[Z].A * amu;
	ene /= rho;
}

int main() {
	for (real x = 1; x < 1.0e+39; x *= 10.0) {
		const real T = 1e+3;
		const real c0 = two * std::pow(two * M_PI * me * kb / (h * h), 1.5);
		const real beta = kb * T / (me * c * c);
		const real eta = double(electron_chemical_potential(x, T));
		const real mu = kb * T * eta;
		const real p = electron_pressure(eta, beta);
		const real e1 = eta;
		const real rho = amu * 0.75 * x;
		const real e2 = 10.0 * rho;
		const real e3 = 2.17e-11 / (kb * T);
		std::printf("%e %e %16.12e\n", rho, x, eta);
	}
	return 0;
}


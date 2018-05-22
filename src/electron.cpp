/*
 * electron.cpp
 *
 *  Created on: May 21, 2018
 *      Author: dmarce1
 */

#include "electron.hpp"
#include "fermi_dirac.hpp"
#include "physcon.hpp"
#include <cmath>

real electron_eos::P(real ne, real T) const {
	return ne * helmholtz(std::log(ne), std::log(T), 1, 0);
}

real electron_eos::eps(real ne, real T) const {
	const real log_ne = std::log(ne);
	const real log_T = std::log(T);
	const real ST = -helmholtz(log_ne, log_T, 0, 1);
	return (helmholtz(log_ne, log_T) + ST) / ne;
}

real electron_eos::dP_dne(real ne, real T) const {
	const real log_ne = std::log(ne);
	const real log_T = std::log(T);
	return helmholtz(log_ne, log_T, 1, 0) + helmholtz(log_ne, log_T, 2, 0);
}

real electron_eos::dP_T(real ne, real T) const {
	return ne * helmholtz(std::log(ne), std::log(T), 1, 1) / T;
}

real electron_eos::deps_dne(real ne, real T) const {
	const real log_ne = std::log(ne);
	const real log_T = std::log(T);
	return (helmholtz(log_ne, log_T, 1, 0) - helmholtz(log_ne, log_T, 1, 1)
			- helmholtz(log_ne, log_T) + helmholtz(log_ne, log_T, 0, 1))
			/ (ne * ne);
}

real electron_eos::deps_dT(real ne, real T) const {

}

electron_eos::electron_eos() :
		helmholtz([](real log_n, real log_T ) {
			const real n = std::exp(log_n);
			const real T = std::exp(log_T);
			return electron_free_energy(n,T);
		}, std::log(n_min), std::log(n_max), std::log(T_min), std::log(T_max),
				1.0e-3) {
	FILE* fp = fopen("electron_eos.dat", "wb");
	helmholtz.save(fp);
	fclose(fp);
	fp = fopen("eos.dat", "wt");
	for (real n = n_min; n < n_max; n *= 2) {
		for (real T = T_min; T < T_max; T *= 2) {
			fprintf(fp, "%e %e %e %e\n", n, T,
					helmholtz(std::log(n), std::log(T)),
					electron_free_energy(n, T));
		}
	}
	fclose(fp);
}

real electron_eos::N_ele(real eta, real beta) {
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real f1 = FermiDirac(one / two, eta, beta);
	const real f2 = FermiDirac(three / two, eta, beta);
	return c0 * std::pow(beta, 1.5) * (f1 + beta * f2);
}

real electron_eos::dN_ele_deta(real eta, real beta) {
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real df1_deta = dFermiDirac_deta(one / two, eta, beta);
	const real df2_deta = dFermiDirac_deta(three / two, eta, beta);
	return c0 * std::pow(beta, 1.5) * (df1_deta + beta * df2_deta);
}

real electron_eos::N_pos(real eta, real beta) {
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

real electron_eos::dN_pos_deta(real eta, real beta) {
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

real electron_eos::electron_chemical_potential(real ne, real beta) {
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
		//	printf( "    %e %e %e\n", ne, T, eta);
	} while (err > 1.0e-12 && std::abs(deta) > 1.0e-14);
//	printf( "%e %e %e\n", ne, T, eta);
	return eta;
}

real electron_eos::electron_free_energy(real ne, real T) {
	const real beta = kb * T / me / c / c;
	const real etae = electron_chemical_potential(ne, beta);
	static const real c0 = real(64) * std::atan(1) * std::sqrt(2)
			* std::pow(me * c / h, three) / real(3) * me * c * c;
	static const real c1 = real(32) * std::atan(1) * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real etap = -etae - two / beta;
	real f12p, f32p, f52p;
	const real f12e = FermiDirac(one / two, etae, beta);
	const real f32e = FermiDirac(three / two, etae, beta);
	const real f52e = FermiDirac(five / two, etae, beta);
	if (etap > -200.0) {
		f12p = FermiDirac(one / two, etap, beta);
		f32p = FermiDirac(three / two, etap, beta);
		f52p = FermiDirac(five / two, etap, beta);
	} else {
		f12p = f32p = f52p = 0.0;
	}
	const real pele = c0 * std::pow(beta, five / two)
			* (f32e + beta / two * f52e);
	const real ppos = c0 * std::pow(beta, five / two)
			* (f32p + beta / two * f52p);
	const real nele = c1 * std::pow(beta, three / two) * (f12e + beta * f32e);
	const real npos = c1 * std::pow(beta, three / two) * (f12p + beta * f32p);
	const real free_energy = (kb * T * (etae * nele + etap * npos)
			- (pele + ppos)) / ne;
	return free_energy;
}

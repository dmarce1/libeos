#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <vector>

#include "elements.hpp"
#include "physcon.hpp"
#include "real.hpp"

real adaptive_integration(const std::function<real(real)>& func, real a, real b,
		const real toler = 1.0e-6, int depth = 0) {
	const real dx_c = (b - a) / two;
	const real dx_f = dx_c / two;
	const real x1 = a;
	const real x2 = three * a / four + b / four;
	const real x3 = a / two + b / two;
	const real x4 = a / four + three * b / four;
	const real x5 = b;
	const real f1 = func(x1);
	const real f2 = func(x2);
	const real f3 = func(x3);
	const real f4 = func(x4);
	const real f5 = func(x5);
	const real sum_c = (f1 + four * f3 + f5) * three_inv * dx_c;
	const real sum_f1 = (f1 + four * f2 + f3) * three_inv * dx_f;
	const real sum_f2 = (f3 + four * f4 + f5) * three_inv * dx_f;
	const real sum_f = sum_f1 + sum_f2;
	const real err = std::abs(std::log(std::abs(sum_c / sum_f)));
	if ((depth < 4)
			|| (((err > toler) || (sum_f * sum_c < zero)) && (dx_c > 1.0e-16))) {
		return adaptive_integration(func, x1, x3, toler, depth + 1)
				+ adaptive_integration(func, x3, x5, toler, depth + 1);
	} else {
		return sum_f;
	}

}

real fd(const real k, const real eta, const real beta) {
	const auto func =
			[k,eta,beta](real x) {
				if( x != zero && x != real(1)) {
					const real num1 = std::pow(x, k) * std::sqrt(one + half * beta * x);
					const real num2 = std::pow(x, -two - k) * std::sqrt(one + half * beta / x);
					const real den1 = one + std::exp(x - eta);
					const real den2 = one + std::exp(one / x - eta);
					return num1 / den1 + num2 / den2;
				} else {
					return zero;
				}
			};
	if (eta > zero && eta != one) {
		const real a = zero;
		const real b = std::min(eta, one / eta);
		const real c = one;
		return adaptive_integration(func, a, b)
				+ adaptive_integration(func, b, c);
	} else {
		return adaptive_integration(func, zero, one);
	}
}

real dfd_deta(const real k, const real eta, const real beta) {
	const auto func =
			[k,eta,beta](real x) {
				if( x != zero && x != real(1)) {
					const real num1 = std::pow(x, k) * std::sqrt(one + half * beta * x);
					const real num2 = std::pow(x, -two - k) * std::sqrt(one + half * beta / x);
					const real den1 = two + std::cosh(x-eta);
					const real den2 = two + std::cosh(one/x-eta);
					real sum = num1 / den1 + num2 / den2;
					return sum;
				} else {
					return zero;
				}
			};
	if (eta > zero && eta != one) {
		const real a = zero;
		const real b = std::min(eta, one / eta);
		const real c = one;
		return adaptive_integration(func, a, b)
				+ adaptive_integration(func, b, c);
	} else {
		return adaptive_integration(func, zero, one);
	}
}

real N_ele(real eta, real beta) {
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real f1 = fd(one / two, eta, beta);
	const real f2 = fd(three / two, eta, beta);
	return c0 * std::pow(beta, 1.5) * (f1 + f2);
}

real dN_ele_deta(real eta, real beta) {
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real df1_deta = dfd_deta(one / two, eta, beta);
	const real df2_deta = dfd_deta(three / two, eta, beta);
	return c0 * std::pow(beta, 1.5) * (df1_deta + df2_deta);
}

real N_pos(real eta, real beta) {
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real eta2 = -eta - two / beta;
	const real f1 = fd(one / two, eta2, beta);
	const real f2 = fd(three / two, eta2, beta);
	return c0 * std::pow(beta, 1.5) * (f1 + beta * f2);
}

real dN_pos_deta(real eta, real beta) {
	static const real c0 = real(8) * M_PI * std::sqrt(2)
			* std::pow(me * c / h, three);
	const real eta2 = -eta - two / beta;
	const real df1_deta = -dfd_deta(one / two, eta2, beta);
	const real df2_deta = -dfd_deta(three / two, eta2, beta);
	return c0 * std::pow(beta, 1.5) * (df1_deta + beta * df2_deta);
}

real electron_pressure(real eta, real beta) {
	static const real c0 = real(64) * std::atan(1) * std::sqrt(2)
			* std::pow(me * c / h, three) / real(3) * me * c * c;
	const real eta2 = -eta - two / beta;

	real pele, ppos;

	const real f1 = fd(three / two, eta, beta);
	const real f2 = fd(five / two, eta, beta);
	const real f3 = fd(three / two, eta2, beta);
	const real f4 = fd(five / two, eta2, beta);

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

	real eta = zero;
	real f, deta, err;
	do {
		f = F(eta);
		real df_deta = dF_deta(eta);
		deta = -f / df_deta;
		eta += deta;
		if (eta == 0.0) {
			err = 1.0;
		} else {
			err = std::abs(deta / eta);
		}
		//printf("%16.10e %16.10e %16.10e %16.10e %16.10e \n", double(err),
		//		double(eta), double(f), double(deta), double(df_deta));
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
	for (real x = 1.0; x < 1.0e+40; x *= 10.0) {
		const real T = 1000.0;
		const real c0 = two * std::pow(two * M_PI * me * kb / (h * h), 1.5);
		const real beta = kb * T / (me * c * c);
		const real eta = double(electron_chemical_potential(x, 1000.0));
		const real mu = kb * T * eta;
		const real p = electron_pressure(eta, beta);
		const real e1 = std::exp(-eta);
		const real e2 = c0 * std::pow(T, 1.5) / x;
		std::printf("%e %e %e %e %e\n", x, e1, e2, e1 / e2, p / kb / T / x);
	}
	return 0;
}


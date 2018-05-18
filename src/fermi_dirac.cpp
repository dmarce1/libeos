#include "real.hpp"

#include <numeric>
#include <cmath>
#include <functional>
#include "physcon.hpp"
#include <vector>
#include <omp.h>

constexpr int N = 100;

constexpr real quadrature_weights[N] = { 7.3463449050567173E-4,
		0.0017093926535181052395, 0.0026839253715534824194,
		0.0036559612013263751823, 0.0046244500634221193511,
		0.0055884280038655151572, 0.0065469484508453227642,
		0.0074990732554647115788, 0.0084438714696689714026,
		0.0093804196536944579514, 0.0103078025748689695858,
		0.0112251140231859771172, 0.0121314576629794974077,
		0.013025947892971542286, 0.013907710703718772688,
		0.014775884527441301769, 0.015629621077546002724,
		0.0164680861761452126431, 0.017290460568323582439,
		0.0180959407221281166644, 0.018883739613374904553,
		0.0196530874944353058654, 0.0204032326462094327668,
		0.0211334421125276415427, 0.021843002416247386314,
		0.0225312202563362727018, 0.02319742318525412162249,
		0.0238409602659682059626, 0.02446120270795705272,
		0.02505754448157958970376, 0.025629402910208116076,
		0.026176219239545676342, 0.02669745918357096266,
		0.0271926134465768801365, 0.0276611982207923882942,
		0.0281027556591011733176, 0.0285168543223950979909,
		0.0289030896011252031349, 0.0292610841106382766201,
		0.029590488059912642512, 0.029890979593332830917,
		0.030162265105169144919, 0.03040407952645482001651,
		0.0306161865839804484965, 0.030798379031152590428,
		0.030950478850490988234, 0.031072337427566516588,
		0.0311638356962099067838, 0.0312248842548493577324,
		0.031255423453863356948, 0.0312554234538633569476,
		0.0312248842548493577324, 0.031163835696209906784,
		0.0310723374275665165878, 0.0309504788504909882341,
		0.0307983790311525904277, 0.0306161865839804484965,
		0.030404079526454820017, 0.030162265105169144919,
		0.029890979593332830917, 0.0295904880599126425118,
		0.0292610841106382766201, 0.0289030896011252031349,
		0.0285168543223950979909, 0.0281027556591011733176,
		0.027661198220792388294, 0.0271926134465768801365,
		0.0266974591835709626604, 0.0261762192395456763423,
		0.025629402910208116076, 0.025057544481579589704,
		0.02446120270795705272, 0.023840960265968205963,
		0.023197423185254121622, 0.0225312202563362727018,
		0.021843002416247386314, 0.021133442112527641543,
		0.0204032326462094327668, 0.0196530874944353058654,
		0.0188837396133749045529, 0.018095940722128116664,
		0.0172904605683235824393, 0.0164680861761452126431,
		0.0156296210775460027239, 0.0147758845274413017689,
		0.013907710703718772688, 0.0130259478929715422856,
		0.012131457662979497408, 0.0112251140231859771172,
		0.0103078025748689695858, 0.0093804196536944579514,
		0.008443871469668971403, 0.0074990732554647115788,
		0.0065469484508453227642, 0.005588428003865515157,
		0.004624450063422119351, 0.0036559612013263751823,
		0.0026839253715534824194, 0.0017093926535181052395,
		7.346344905056717304E-4

};
constexpr real legendre_roots[N] = { -0.9997137267734412336782,
		-0.9984919506395958184002, -0.996295134733125149186,
		-0.993124937037443459652, -0.9889843952429917480044,
		-0.983877540706057015496, -0.977809358486918288554,
		-0.970785775763706331931, -0.962813654255815527294,
		-0.9539007829254917428493, -0.944055870136255977963,
		-0.933288535043079545924, -0.921609298145333952667,
		-0.9090295709825296904671, -0.8955616449707269866985,
		-0.8812186793850184155733, -0.8660146884971646234107,
		-0.8499645278795912842934, -0.8330838798884008235429,
		-0.815389238339176254394, -0.79689789239031447639,
		-0.7776279096494954756276, -0.7575981185197071760357,
		-0.7368280898020207055124, -0.71533811757305644646,
		-0.6931491993558019659487, -0.670283015603141015803,
		-0.6467619085141292798326, -0.622608860203707771604,
		-0.5978474702471787212648, -0.5725019326213811913169,
		-0.546597012065094167468, -0.520158019881763056647,
		-0.493210789208190933569, -0.4657816497733580422492,
		-0.437897402172031513109, -0.4095852916783015425289,
		-0.3808729816246299567634, -0.3517885263724217209723,
		-0.3223603439005291517225, -0.2926171880384719647376,
		-0.2625881203715034791689, -0.2323024818449739696495,
		-0.2017898640957359972361, -0.1710800805386032748875,
		-0.140203137236113973208, -0.109189203580061115003,
		-0.078068582813436636695, -0.046871682421591631615,
		-0.0156289844215430828722, 0.0156289844215430828722,
		0.04687168242159163161492, 0.07806858281343663669482,
		0.1091892035800611150034, 0.1402031372361139732075,
		0.171080080538603274888, 0.2017898640957359972361,
		0.23230248184497396965, 0.2625881203715034791689,
		0.2926171880384719647376, 0.3223603439005291517225,
		0.351788526372421720972, 0.3808729816246299567634,
		0.409585291678301542529, 0.437897402172031513109,
		0.465781649773358042249, 0.4932107892081909335693,
		0.5201580198817630566468, 0.546597012065094167468,
		0.5725019326213811913169, 0.597847470247178721265,
		0.6226088602037077716042, 0.6467619085141292798326,
		0.6702830156031410158026, 0.6931491993558019659487,
		0.71533811757305644646, 0.736828089802020705512,
		0.7575981185197071760357, 0.777627909649495475628,
		0.7968978923903144763896, 0.815389238339176254394,
		0.833083879888400823543, 0.8499645278795912842934,
		0.8660146884971646234107, 0.8812186793850184155733,
		0.895561644970726986699, 0.9090295709825296904671,
		0.921609298145333952667, 0.9332885350430795459243,
		0.9440558701362559779628, 0.953900782925491742849,
		0.9628136542558155272937, 0.970785775763706331931,
		0.977809358486918288554, 0.9838775407060570154961,
		0.988984395242991748004, 0.993124937037443459652,
		0.9962951347331251491861, 0.9984919506395958184002,
		0.999713726773441233678 };

real integrate(const std::function<real(real)>& f, real a, real b) {
	const real dx = b - a;
	real sum = real(0);
	for (int i = 0; i < N; i++) {
		const real x = a + dx * 0.5 * (legendre_roots[i] + 1.0);
		sum += f(x) * quadrature_weights[i];
	}
	return sum * dx * 0.5;
}

std::function<real(real)> fd(real k, real eta, real beta) {
	return [k,eta,beta]( real x ) {
//		printf( "%e %e %e\n", x, eta, x - eta);
		return std::pow(x,k) * std::sqrt(1.0 + beta * x / 2.0) / (std::exp(x-eta)+1.0);
	};
}

std::function<real(real)> fd2(real k, real eta, real beta) {
	return [k,eta,beta]( real x ) {
//		printf( "%e %e %e\n", x, eta, x - eta);
		return std::pow(x,2*k+1) * std::sqrt(1.0 + beta * x*x / 2.0) / (std::exp(x*x-eta)+1.0);
	};
}

std::function<real(real)> fd_bigeta(real k, real eta, real beta) {
	return [k,eta,beta]( real x ) {
		if( x <= eta ) {
			return std::pow(x,k) * std::sqrt(1.0 + beta * x / 2.0);
		} else {
			return 0.0;
		}
	};
}

std::function<real(real)> dfd2_deta(real k, real eta, real beta) {
	return [k,eta,beta]( real x ) {
		const real exp = std::exp(x*x-eta);
		return std::pow(x,2*k+1) * std::sqrt(1.0 + beta * x*x / 2.0) / (exp+1.0/exp+2.0);
	};
}

std::function<real(real)> dfd_deta(real k, real eta, real beta) {
	return [k,eta,beta]( real x ) {
		const real exp = std::exp(x-eta);
		return std::pow(x,k) * std::sqrt(1.0 + beta * x / 2.0) / (exp+1.0/exp+2.0);
	};
}

real dfd_deta_bigeta(real k, real eta, real beta) {
	return std::pow(eta, k) * std::sqrt(1.0 + beta * eta / 2.0) / 4.0;
}

real integrate(const std::function<real(real)>& f, real a, real b, int n) {
	const real dx = (b - a) / n;
	const int nthreads = omp_get_max_threads();
	std::vector<real> sums(nthreads, 0.0);
	real x = a;
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		const real x = a + i * dx;
		const real this_sum = integrate(f, x, x + dx);
		sums[omp_get_thread_num()] += this_sum;
	}
	return std::accumulate(sums.begin(), sums.end(), 0.0);
}

constexpr real width = 30.0;
int M = 100;

real fd_ne(real eta, real beta) {
	real w = std::max(width, 1.0e-10 * eta);
	real a = 0.0;
	real b = std::max(w, eta + w);
	const real c0 = real(8) * M_PI * std::sqrt(2) * std::pow(me * c / h, three)
			* std::pow(beta, 1.5);
	const auto func = [eta,beta](real x ) {
		const real a = std::sqrt(x + 0.5 * x * x * beta);
		const real b = (1.0 + x * beta);
		const real c = std::sinh(eta+1.0/beta) * 0.5;
		const real d = std::cosh(eta+1.0/beta) + std::cosh(x+1.0/beta);
		return a * b * c / d;
	};
	return integrate(func, a, b, M);

}
/*
real FermiDirac(real k, real eta, real beta) {
	real w = std::max(width, 1.0e-10 * eta);
	real a = 0.0;
	real b = std::max(w, eta + w);
	if (eta > 1.0e+10) {
		return integrate(fd_bigeta(k, eta, beta), a, b, M);
	} else {
		return integrate(fd(k, eta, beta), a, b, M);
	}
}

real dFermiDirac_deta(real k, real eta, real beta) {
	if (eta > 1.0e+10) {
		return dfd_deta_bigeta(k, eta, beta);
	} else {
		real w = std::max(width, 1.0e-10 * eta);
		real a = std::max(0.0, eta - w);
		real b = std::max(w, eta + w);
		return integrate(dfd_deta(k, eta, beta), a, b, M);
	}
}
*/
void aparcio_cuts(real eta, real& x1, real& x2, real& x3) {
	constexpr real D = 3.3609;
	constexpr real sigma = 9.1186e-2;
	constexpr real a1 = 6.7774;
	constexpr real b1 = 1.1418;
	constexpr real c1 = 2.9826;
	constexpr real a2 = 3.7601;
	constexpr real b2 = 9.3719e-2;
	constexpr real c2 = 2.1063e-2;
	constexpr real d2 = 3.1084e+1;
	constexpr real e2 = 1.0056;
	constexpr real a3 = 7.5669;
	constexpr real b3 = 1.1695;
	constexpr real c3 = 7.5416e-1;
	constexpr real d3 = 6.6558;
	constexpr real e3 = -1.2819e-1;
	real xsi;
	if (eta < 7.739391e+03) {
		xsi = std::log(real(1) + std::exp(sigma * (eta - D))) / sigma;
	} else {
		xsi = eta;
	}
	const real xsi2 = xsi * xsi;
	real Xa = a1 + b1 * xsi + c1 * xsi2;
	real Xb = a2 + b2 * xsi + c2 * d2 * xsi2;
	real Xc = a3 + b3 * xsi + c3 * d3 * xsi2;
	Xa /= real(1) + c1 * xsi;
	Xb /= real(1) + e2 * xsi + c2 * xsi2;
	Xc /= real(1) + e3 * xsi + c3 * xsi2;
	x1 = Xa - Xb;
	x2 = Xa;
	x3 = Xa + Xc;
}

real FermiDirac(real k, real eta, real beta) {
	real x0, x1, x2, x3, x4;
	x0 = zero;
	aparcio_cuts(eta, x1, x2, x3);
	x4 = x3 + width;
	constexpr int N = 100;
	const real sum1 = integrate(fd2(k, eta, beta), std::sqrt(x0), std::sqrt(x1),
			M);
	const real sum2 = integrate(fd2(k, eta, beta), std::sqrt(x1), std::sqrt(x2),
			M);
	const real sum3 = integrate(fd(k, eta, beta), x2, x3, 1);
	const real sum4 = integrate(fd(k, eta, beta), x3, x4, 1);
	return two * (sum1 + sum2) + sum3 + sum4;
}

real dFermiDirac_deta(real k, real eta, real beta) {
	if (eta > 1.0e+10) {
		return dfd_deta_bigeta(k, eta, beta);
	} else {
		real x0, x1, x2, x3, x4;
		x0 = zero;
		aparcio_cuts(eta, x1, x2, x3);
		x4 = x3 + width;
		constexpr int N = 100;
		const real sum1 = integrate(dfd2_deta(k, eta, beta), std::sqrt(x0),
				std::sqrt(x1), 1);
		const real sum2 = integrate(dfd2_deta(k, eta, beta), std::sqrt(x1),
				std::sqrt(x2), 1);
		const real sum3 = integrate(dfd_deta(k, eta, beta), x2, x3, 1);
		const real sum4 = integrate(dfd_deta(k, eta, beta), x3, x4, 1);
		return two * (sum1 + sum2) + sum3 + sum4;
	}
}
/*

int main() {
	real x1, x2, x3;
	const real k = 0.5;
	real eta = 1.0e+7;
	real beta = 1.0;
	real last_f1, last_f2;
	real f1 = 1, f2 = 2;
	for (M = 1; M < 1000909000; M *= 2) {
		last_f1 = f1;
		last_f2 = f2;
		f1 = dFermiDirac_deta(k, eta, beta);
		f2 = dFermiDirac2_deta(k, eta, beta);
		real e1 = std::abs(last_f1 - f1) / f1;
		real e2 = std::abs(last_f2 - f2) / f2;
		printf("%i %.14e %.4e %.14e %.4e\n", M, f1, e1, f2, e2);
	}
}
*/

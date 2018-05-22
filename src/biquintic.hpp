/*
 * biquintic.hpp
 *
 *  Created on: May 9, 2018
 *      Author: dmarce1
 */

#ifndef BICUBIC_HPP_
#define BICUBIC_HPP_

#include "real.hpp"
#include <array>
#include <vector>
#include <functional>

class biquintic_table {
private:
	real xmin, xmax, ymin, ymax;
	real L1, L2, Linf;
	real dx, dy;
	int NX, NY;
	std::vector<std::array<real, 36>> C;
public:
	biquintic_table(const std::function<real(real, real)>& func, real _xmin,
			real _xmax, real _ymin, real _ymax, int _NX, int _NY);
	biquintic_table(const std::function<real(real, real)>& func, real _xmin,
			real _xmax, real _ymin, real _ymax, real toler,
			const char* filename = nullptr);

	real operator()(real x, real y, int Dx = 0, int Dy = 0) const;
	bool in_range(real x, real y) const;

	void save(FILE* fp) const {
		const auto& f = fwrite;
		f(&xmin, sizeof(real), 1, fp);
		f(&xmax, sizeof(real), 1, fp);
		f(&ymin, sizeof(real), 1, fp);
		f(&ymax, sizeof(real), 1, fp);
		f(&dx, sizeof(real), 1, fp);
		f(&dy, sizeof(real), 1, fp);
		f(&NX, sizeof(int), 1, fp);
		f(&NY, sizeof(int), 1, fp);
		f(C.data(), sizeof(std::array<real, 36>), NX * NY, fp);
	}
	biquintic_table(FILE* fp) {
		const auto& f = fread;
		f(&xmin, sizeof(real), 1, fp);
		f(&xmax, sizeof(real), 1, fp);
		f(&ymin, sizeof(real), 1, fp);
		f(&ymax, sizeof(real), 1, fp);
		f(&dx, sizeof(real), 1, fp);
		f(&dy, sizeof(real), 1, fp);
		f(&NX, sizeof(int), 1, fp);
		f(&NY, sizeof(int), 1, fp);
		C.resize(NX * NY);
		f(C.data(), sizeof(std::array<real, 36>), NX * NY, fp);
	}
};

#endif /* BICUBIC_HPP_ */

/*
 * bicubic.hpp
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

class bicubic_table {
private:
	std::vector<std::array<real, 16>> C;
	real xmin, xmax, ymin, ymax;
	int NX, NY;
	real L1, L2, Linf;
	real dx, dy;
public:
	bicubic_table(const std::function<real(real, real)>& func, real _xmin,
			real _xmax, real _ymin, real _ymax, int _NX, int _NY);
	bicubic_table(const std::function<real(real, real)>& func, real _xmin,
			real _xmax, real _ymin, real _ymax, real toler,
			const char* filename = nullptr);

	real operator()(real x, real y) const;
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
		f(C.data(), sizeof(std::array<real, 16>), NX * NY, fp);
	}
	bicubic_table(FILE* fp) {
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
		f(C.data(), sizeof(std::array<real, 16>), NX * NY, fp);
	}
};

#endif /* BICUBIC_HPP_ */

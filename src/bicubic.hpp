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
	real L1, L2;
	real dx, dy;
public:
	bicubic_table(const std::function<real(real, real)>& func, real _xmin,
			real _xmax, real _ymin, real _ymax, int _NX, int _NY);

	real operator()(real x, real y);
};


#endif /* BICUBIC_HPP_ */

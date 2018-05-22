/*
 * electron.hpp
 *
 *  Created on: May 21, 2018
 *      Author: dmarce1
 */

#ifndef ELECTRON_HPP_
#define ELECTRON_HPP_

#include "biquintic.hpp"

class electron_eos {
private:
	biquintic_table helmholtz;
	static real N_ele(real eta, real beta);
	static real dN_ele_deta(real eta, real beta);
	static real N_pos(real eta, real beta);
	static real dN_pos_deta(real eta, real beta);
	static real electron_chemical_potential(real ne, real beta);
	static real electron_free_energy(real ne, real T);
	static constexpr real n_min = 1.0e-4;
	static constexpr real n_max = 1.0e+35;
	static constexpr real T_min = 1.0e+1;
	static constexpr real T_max = 1.0e+12;
public:
	electron_eos();
	real P(real ne, real T) const;
	real eps(real ne, real T) const;
	real dP_dne(real ne, real T) const;
	real dP_T(real ne, real T) const;
	real deps_dne(real ne, real T) const;
	real deps_dT(real ne, real T) const;

};


#endif /* ELECTRON_HPP_ */

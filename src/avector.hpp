/*
 * avector.hpp
 *
 *  Created on: May 21, 2018
 *      Author: dmarce1
 */

#ifndef AVECTOR_HPP_
#define AVECTOR_HPP_

#include "real.hpp"
#include <array>

template<int N>
class avector {
	std::array<real, N> a;
public:
	avector() {
	}
	avector(real d) {
		for (int i = 0; i < N; i++) {
			a[i] = d;
		}
	}
	real operator[](int i) const {
		return a[i];
	}
	real& operator[](int i) {
		return a[i];
	}
	avector<N> operator*(real d) const {
		avector<N> b;
		for (int i = 0; i < N; i++) {
			b[i] = d * a[i];
		}
		return b;
	}
	avector<N> operator/(real d) const {
		avector<N> b;
		for (int i = 0; i < N; i++) {
			b[i] = a[i] / d;
		}
		return b;
	}
	avector<N> operator-(const avector<N>& d) const {
		avector<N> b;
		for (int i = 0; i < N; i++) {
			b[i] = a[i] - d[i];
		}
		return b;
	}
	avector<N> operator/(const avector<N>& d) const {
		avector<N> b;
		for (int i = 0; i < N; i++) {
			b[i] = a[i] / d[i];
		}
		return b;
	}
	avector<N> operator+(const avector<N>& d) const {
		avector<N> b;
		for (int i = 0; i < N; i++) {
			b[i] = a[i] + d[i];
		}
		return b;
	}
	avector<N> operator+=(const avector<N>& other) {
		for (int i = 0; i < N; i++) {
			a[i] += other[i];
		}
		return *this;
	}

};

template<int N>
real operator+=(real& a, const avector<N>& v) {
	for (int i = 0; i < N; i++) {
		a += v[i];
	}
	return a;
}

template<int N>
avector<N> operator*(real d, const avector<N>& v) {
	avector<N> b;
	for (int i = 0; i < N; i++) {
		b[i] = d * v[i];
	}
	return b;
}
template<int N>
avector<N> operator*(const avector<N>& v1, const avector<N>& v2) {
	avector<N> b;
	for (int i = 0; i < N; i++) {
		b[i] = v2[i] * v1[i];
	}
	return b;
}

namespace std {
template<int N>
avector<N> abs(const avector<N>& a) {
	avector<N> b;
	for (int i = 0; i < N; i++) {
		b[i] = abs(a[i]);
	}
	return b;
}

template<int N>
real max(const real& a, const avector<N>& b) {
	real c = a;
	for (int i = 0; i < N; i++) {
		c = max(c, b[i]);
	}
	return c;
}


template<int N>
real min(const real& a, const avector<N>& b) {
	real c = a;
	for (int i = 0; i < N; i++) {
		c = min(c, b[i]);
	}
	return c;
}


template<int N>
real min(const avector<N>& a, const avector<N>& b) {
	real c = a[0];
	for (int i = 0; i < N; i++) {
		c = min(c, b[i]);
		c = min(c, a[i]);
	}
	return c;
}
}

#endif /* AVECTOR_HPP_ */

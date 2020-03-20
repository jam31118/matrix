#include "../include/matrix.hh"

#include <complex>

template int tridiag_mul_forward<double, double>(
		double *, double *, double *, 
		double *, double *, long);

template int tridiag_mul_forward<double, std::complex<double>>(
		double *, double *, double *,
		std::complex<double> *, std::complex<double> *, long);

template int tridiag_mul_forward<std::complex<double>, std::complex<double>>(
		std::complex<double> *, std::complex<double> *, std::complex<double> *,
		std::complex<double> *, std::complex<double> *, long);


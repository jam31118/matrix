#include <iostream>
#include <complex>
#include "matrix.hh"

typedef std::complex<double> mt;

int main() {

  long N = 11;
  mt alpha[N], beta[N-1], gamma[N-1], v[N], out[N], v_reconstructed[N];
  long i;
  alpha[N-1] = 1;
  for (i=0; i<N-1; i++) {
    alpha[i] = 1;
    beta[i] = 0;
    gamma[i] = 1;
    v[i] = 0.1 * i;
  }
  v[N-1] = 0.1 * (N-1);

  // Test mat_vec_mul_tridiag()
  mat_vec_mul_tridiag(alpha, beta, gamma, v, out, N);

  for (i=0; i<N; i++) {
    std::cout << out[i] << " ";
  }
  std::cout << std::endl;


  // Test gaussian_elimination_tridiagonal()
  gaussian_elimination_tridiagonal(alpha, beta, gamma, v_reconstructed, out, N); 
 
  for (i=0; i<N; i++) {
    std::cout << v[i] << " ";
  }
  std::cout << "= v";
  std::cout << std::endl;

  for (i=0; i<N; i++) {
    std::cout << v_reconstructed[i] << " ";
  }
  std::cout << "= A^{-1} A v";
  std::cout << std::endl;

  return 0;
}


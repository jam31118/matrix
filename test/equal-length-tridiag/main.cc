#include <iostream>
#include <complex>
#include "matrix.hh"

typedef std::complex<double> mt;

int main() {

  //// Prepare arrays
  // [NOTE] all tridiagonals (`ld`, `d`, `ud`) have same length `N`
  long N = 11;
  mt d[N], ld[N], ud[N], v[N], b[N], v_reconstructed[N];
  long i;
//  d[N-1] = 1;
  for (i=0; i<N; i++) {
    d[i] = 1;
    ld[i] = 0;
    ud[i] = 1;
    v[i] = 0.1 * i;
  }
//  v[N-1] = 0.1 * (N-1);
  ld[0] = 0; ud[N-1] = 0;



  //// Forward-Backward operation
  // Test mat_vec_mul_tridiag()
  tridiag_mul_forward(ld, d, ud, v, b, N);

  // Test gaussian_elimination_tridiagonal()
  tridiag_mul_backward(ld, d, ud, v_reconstructed, b, N);



  //// Print results
  std::cout << "v[] = ";
  for (i=0; i<N; i++) {
    std::cout << v[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "b[] = ";
  for (i=0; i<N; i++) {
    std::cout << b[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "= A^{-1} A v = ";
  for (i=0; i<N; i++) {
    std::cout << v_reconstructed[i] << " ";
  }
  std::cout << std::endl;

  return 0;
}


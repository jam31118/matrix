#ifndef _MATRIX_HH_
#define _MATRIX_HH_

#include <cstdio>
#include <cstdlib>
#include <iostream>

template <typename T, typename type_vec>
int tridiag_mul_forward(T *ld, T *d, T *ud, type_vec *v, type_vec *b, long N) {
  
  long offset = 1, num_of_elements_in_loop = N - 2;

  T *p_d = d + offset, 
    *p_ld = ld + offset, 
    *p_ud = ud + offset;
  type_vec 
		*p_v = v + offset, 
    *p_b = b + offset;

  T *p_d_max = p_d + num_of_elements_in_loop;
 
  for ( ; p_d < p_d_max; ++p_d, ++p_ld, ++p_ud, ++p_v, ++p_b ) {
    *p_b = (*p_ld) * (*(p_v - 1)) + (*p_d) * (*p_v) + (*p_ud) * (*(p_v + 1)); 
  }

  b[0] = d[0] * v[0] + ud[0] * v[1];
  b[N-1] = ld[N-1] * v[N-2] + d[N-1] * v[N-1];

  return 0;
}

template <typename T, typename type_vec>
int tridiag_mul_backward(T *ld, T *d, T *ud, type_vec *v, type_vec *b, long N) {

  // Allocate arrays for intermediate results
  T *delta = (T *) malloc(sizeof(T) * N);
  
  long i;
  T *p_d, *p_ld, *p_ud, *p_delta; 
  T ld_temp;
	type_vec *p_v, *p_b;

  // copy the first element as a start
  *delta = *d;  
  *v = *b / (*delta);

  // iteration to bigger index
  for ( p_d=d+1, p_ld=ld+1, p_ud=ud, p_delta=delta, p_b=b+1, p_v=v, i=0; i < N-1;
      ++p_d, ++p_ld, ++p_ud, ++p_delta, ++p_v, ++i, ++p_b ) {
    
    ld_temp = *p_ld;
    *(p_delta+1) = *(p_d) - (*p_ud) * ld_temp / (*p_delta);
    *(p_v+1) = (*p_b - *p_v * ld_temp) / (*(p_delta+1)); 
    // [NOTE] Using `p_v+1` in place of `p_b` is not valid since `p_v+1` hasn't been set yet.

  }

  // iteration to smaller index
  // [NOTE] Check that the address pointed by `p_v`, `p_ud`, `p_delta` at previous loop can be used.
  for (--p_ud, --p_delta; i > 0; --p_ud, --p_delta, --p_v, --i) {
//    std::cout << "i: " << i << " *(p_v-1): " << *(p_v-1) << " *p_v * (*p_ud) / (*p_delta): " << (*p_v * (*p_ud) / (*p_delta)) << std::endl;
//    std::cout << "equal? : " << (*(p_v-1) == (*p_v * (*p_ud) / (*p_delta))) << "  diff: " << (*(p_v-1) - (*p_v * (*p_ud) / (*p_delta))) <<  std::endl;
    *(p_v-1) = *(p_v-1) - *p_v * (*p_ud) / (*p_delta);
    // [NOTE] Using `p_v-1` is valid since it has already been set.
  }
  
  //// Free allocated memory
  free(delta);

  return 0;
}




template <class T>
int mat_vec_mul_tridiag(T *alpha, T *beta, T *gamma, T *v, T *out, long N) {
  
  long offset = 1, num_of_elements_in_loop = N - 2;

  T *p_alpha = alpha + offset, 
    *p_beta = beta + offset, 
    *p_gamma = gamma + offset, 
    *p_v = v + offset, 
    *p_out = out + offset;

  T *p_alpha_max = p_alpha + num_of_elements_in_loop;
 
  for ( ; p_alpha < p_alpha_max; ++p_alpha, ++p_beta, ++p_gamma, ++p_v, ++p_out ) {
    *p_out = (*(p_beta - 1)) * (*(p_v - 1)) + (*p_alpha) * (*p_v) + (*p_gamma) * (*(p_v + 1)); 
  }

  out[0] = alpha[0] * v[0] + gamma[0] * v[1];
  out[N-1] = beta[N-2] * v[N-2] + alpha[N-1] * v[N-1];

  return 0;
}

template <class T>
int gaussian_elimination_tridiagonal(T *alpha, T *beta, T *gamma, T *v, T *b, long N) {

  // Allocate arrays for intermediate results
  T *delta = (T *) malloc(sizeof(T) * N);
  
  long i;
  T *p_alpha, *p_beta, *p_gamma, *p_b, *p_delta, *p_v;
  T beta_temp;

  // copy the first element as a start
  *delta = *alpha;  
  *v = *b / (*delta);

  // iteration to bigger index
  for ( p_alpha=alpha+1, p_beta=beta, p_gamma=gamma, p_delta=delta, p_b=b+1, p_v=v, i=0; i < N-1;
      ++p_alpha, ++p_beta, ++p_gamma, ++p_delta, ++p_v, ++i, ++p_b ) {
    
    beta_temp = *p_beta;
    *(p_delta+1) = *(p_alpha) - (*p_gamma) * beta_temp / (*p_delta);
    *(p_v+1) = (*p_b - *p_v * beta_temp) / (*(p_delta+1)); 
    // [NOTE] Using `p_v+1` in place of `p_b` is not valid since `p_v+1` hasn't been set yet.

  }

  // iteration to smaller index
  // [NOTE] Check that the address pointed by `p_v`, `p_gamma`, `p_delta` at previous loop can be used.
  for (--p_gamma, --p_delta; i > 0; --p_gamma, --p_delta, --p_v, --i) {
    *(p_v-1) = *(p_v-1) - *p_v * (*p_gamma) / (*p_delta);
    // [NOTE] Using `p_v-1` is valid since it has already been set.
  }
  
  //// Free allocated memory
  free(delta);

  return 0;
}

#endif // _MATRIX_HH_

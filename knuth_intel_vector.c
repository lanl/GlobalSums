#include <x86intrin.h>

static double sum[4] __attribute__ ((aligned (64)));

double do_knuth_sum_v(double* restrict var, long ncells)
{
#ifndef __PGI
   double const zero = 0.0;
   __m256d local_sum = _mm256_broadcast_sd((double const*) &zero);
   __m256d local_correction = _mm256_broadcast_sd((double const*) &zero);

#ifdef __INTEL_COMPILER
      #pragma ivdep
#else
      #pragma simd
#endif
   #pragma vector aligned
   for (long i = 0; i < ncells; i+=4) {
      __m256d u = local_sum;
      __m256d v = _mm256_load_pd(&var[i]) + local_correction;
      __m256d upt = u + v;
      __m256d up = upt - v;
      __m256d vpp = upt - up;
      local_sum = upt;
      local_correction = (u - up) + (v - vpp);
   }

   __m256d sum_v = local_sum + local_correction;
   _mm256_store_pd(sum, sum_v);

   // double to do final sum
   double ud, vd, uptd, upd, vppd;
   double final_sum = 0.0;
   double final_correction = 0.0;

   for (long i = 0; i < 4; i++) {
      ud = final_sum;
      vd = sum_v[i] + final_correction;
      uptd = ud + vd;
      upd = uptd - vd;
      vppd = uptd - upd;
      final_sum = uptd;
      final_correction = (ud - upd) + (vd - vppd);
   }
#else
   double final_sum = 0.0;
#endif

   return(final_sum);
}

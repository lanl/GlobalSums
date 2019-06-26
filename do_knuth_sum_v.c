#include <immintrin.h>
#include <x86intrin.h>

static double sum[4] __attribute__ ((aligned (64)));

double do_knuth_sum_v(double* restrict var, long ncells)
{
   double const zero = 0.0;
   double final_sum = 0.0;
   double final_correction = 0.0;

   __m256d u, v, upt, up, vpp;
   __m256d local_sum, local_correction, sum_v;
   
   local_sum = _mm256_broadcast_sd((double const*) &zero);
   local_correction = _mm256_broadcast_sd((double const*) &zero);
   sum_v = _mm256_broadcast_sd((double const*) &zero);   

   #pragma simd
   #pragma vector aligned
   for (long i = 0; i < ncells; i+=4) {
      u = local_sum;
      v = _mm256_load_pd(&var[i]) + local_correction;
      upt = u + v;
      up = upt - v;
      vpp = upt - up;
      local_sum = upt;
      local_correction = (u - up) + (v - vpp);
   }

   sum_v = local_sum + local_correction;
   _mm256_store_pd(sum, sum_v);

   // double to do final sum
   double ud, vd, uptd, upd, vppd;

   for (long i = 0; i < 4; i++) {
      ud = final_sum;
      vd = sum_v[i] + final_correction;
      uptd = ud + vd;
      upd = uptd - vd;
      vppd = uptd - upd;
      final_sum = uptd;
      final_correction = (ud - upd) + (vd - vppd);
   }

   return(final_sum);
}

#include <x86intrin.h>

static double sum[8] __attribute__ ((aligned (64)));

double do_knuth_sum_intel_v8(double* restrict var, long ncells)
{
   double const zero[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   __m512d local_sum = _mm512_load_pd((double const*) &zero);
   __m512d local_correction = _mm512_load_pd((double const*) &zero);

   #pragma simd
   #pragma vector aligned
   for (long i = 0; i < ncells; i+=8) {
      __m512d u = local_sum;
      __m512d v = _mm512_load_pd(&var[i]) + local_correction;
      __m512d upt = u + v;
      __m512d up = upt - v;
      __m512d vpp = upt - up;
      local_sum = upt;
      local_correction = (u - up) + (v - vpp);
   }

   __m512d sum_v = local_sum + local_correction;
   _mm512_store_pd(sum, sum_v);

   // double to do final sum
   double ud, vd, uptd, upd, vppd;
   double final_sum = 0.0;
   double final_correction = 0.0;

   for (long i = 0; i < 8; i++) {
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

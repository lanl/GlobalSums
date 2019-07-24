#ifdef HAVE_AVX512
#include <x86intrin.h>

static double sum[8] __attribute__ ((aligned (64)));

double do_kahan_sum_intel_v8(double* restrict var, long ncells)
{
   double const zero[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   __m512d local_sum = _mm512_load_pd(zero);
   __m512d local_correction = _mm512_load_pd(zero);

   #pragma simd
   #pragma vector aligned
   for (long i = 0; i < ncells; i+=8) {
       __m512d var_v = _mm512_load_pd(&var[i]);
       __m512d corrected_next_term = var_v + local_correction;
       __m512d new_sum = local_sum + local_correction;
       local_correction = corrected_next_term - (new_sum - local_sum);
       local_sum = new_sum;
   }
   __m512d sum_v;
   sum_v  = local_correction;
   sum_v += local_sum;
   _mm512_store_pd(sum, sum_v);

   struct esum_type{
      double sum;
      double correction;
   } local;
   local.sum = 0.0;
   local.correction = 0.0;

   for (long i = 0; i < 8; i++) {
      double corrected_next_term_s = sum[i] + local.correction;
      double new_sum_s             = local.sum + local.correction;
      local.correction   = corrected_next_term_s - (new_sum_s - local.sum);
      local.sum          = new_sum_s;
   }
   double final_sum = local.sum + local.correction;
   return(final_sum);
}
#endif

#include <x86intrin.h>

static double sum[4] __attribute__ ((aligned (64)));

double do_kahan_sum_v(double* restrict var, long ncells)
{
#ifndef __PGI
   double const zero = 0.0;
   __m256d local_sum = _mm256_broadcast_sd((double const*) &zero);
   __m256d local_correction = _mm256_broadcast_sd((double const*) &zero);
   __m256d var_v;

#ifdef __INTEL_COMPILER
   #pragma ivdep
#else
   #pragma simd
#endif
   #pragma vector aligned
   for (long i = 0; i < ncells; i+=4) {
       var_v = _mm256_load_pd(&var[i]);
       __m256d corrected_next_term = var_v + local_correction;
       __m256d new_sum = local_sum + local_correction;
       local_correction = corrected_next_term - (new_sum - local_sum);
       local_sum = new_sum;
   }
   __m256d sum_v;
   sum_v  = local_correction;
   sum_v += local_sum;
   _mm256_store_pd(sum, sum_v);

   struct esum_type{
      double sum;
      double correction;
   } local;
   local.sum = 0.0;
   local.correction = 0.0;

   for (long i = 0; i < 4; i++) {
      double corrected_next_term_s = sum[i] + local.correction;
      double new_sum_s             = local.sum + local.correction;
      local.correction   = corrected_next_term_s - (new_sum_s - local.sum);
      local.sum          = new_sum_s;
   }
   double final_sum = local.sum + local.correction;
#else
   double final_sum = 0.0;
#endif
   return(final_sum);
}

#include "vectorclass.h"

static double sum[4] __attribute__ ((aligned (64)));

extern "C" {
double do_kahan_sum_agner_v(double* var, long ncells);
}

double do_kahan_sum_agner_v(double* var, long ncells)
{
   Vec4d local_sum(0.0);
   Vec4d local_correction(0.0);
   Vec4d var_v;

   for (long i = 0; i < ncells; i+=4) {
       var_v.load(var+i);
       Vec4d corrected_next_term = var_v + local_correction;
       Vec4d new_sum = local_sum + local_correction;
       local_correction = corrected_next_term - (new_sum - local_sum);
       local_sum = new_sum;
   }
   Vec4d sum_v;
   sum_v  = local_correction;
   sum_v += local_sum;
   sum_v.store(sum);

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
   return(final_sum);
}

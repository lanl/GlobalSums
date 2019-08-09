#include "vectorclass.h"

static double sum[8] __attribute__ ((aligned (64)));

extern "C" {
double do_kahan_sum_agner_v8(double* var, long ncells);
}

double do_kahan_sum_agner_v8(double* var, long ncells)
{
   Vec8d local_sum(0.0);
   Vec8d local_correction(0.0);
   Vec8d var_v;

   int ncells_main=(ncells/8)*8;
   int ncells_remainder=ncells%8;
   for (long i = 0; i < ncells_main; i+=8) {
       var_v.load(var+i);
       Vec8d corrected_next_term = var_v + local_correction;
       Vec8d new_sum = local_sum + local_correction;
       local_correction = corrected_next_term - (new_sum - local_sum);
       local_sum = new_sum;
   }
   if (ncells_remainder > 0) {
       var_v.load_partial(ncells_remainder,var+ncells_main);
       Vec8d corrected_next_term = var_v + local_correction;
       Vec8d new_sum = local_sum + local_correction;
       local_correction = corrected_next_term - (new_sum - local_sum);
       local_sum = new_sum;
   }

   Vec8d sum_v;
   sum_v  = local_correction;
   sum_v += local_sum;
   sum_v.store(sum);

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

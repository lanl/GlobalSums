static double sum[8] __attribute__ ((aligned (64)));

double do_kahan_sum_gcc_v8(double* restrict var, long ncells)
{
   typedef double vec8d __attribute__ ((vector_size(8 * sizeof(double))));

   vec8d local_sum = {0.0};
   vec8d local_correction = {0.0};

   for (long i = 0; i < ncells; i+=8) {
       vec8d var_v = *(vec8d *)&var[i];
       vec8d corrected_next_term = var_v + local_correction;
       vec8d new_sum = local_sum + local_correction;
       local_correction = corrected_next_term - (new_sum - local_sum);
       local_sum = new_sum;
   }
   vec8d sum_v;
   sum_v  = local_correction;
   sum_v += local_sum;
   *(vec8d *)sum = sum_v;

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

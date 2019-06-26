double do_kahan_sum_omp(double* restrict var, long ncells)
{
   struct esum_type{
      double sum;
      double correction;
   };

   double sum = 0.0;

#pragma omp parallel reduction(+:sum)
   {
      double corrected_next_term, new_sum;
      struct esum_type local;

      local.sum = 0.0;
      local.correction = 0.0;
#pragma omp for
      for (long i = 0; i < ncells; i++) {
         corrected_next_term= var[i] + local.correction;
         new_sum      = local.sum + local.correction;
         local.correction   = corrected_next_term - (new_sum - local.sum);
         local.sum          = new_sum;
      }

      sum += local.correction;
#pragma omp barrier
      sum += local.sum;
   }
   return(sum);
}

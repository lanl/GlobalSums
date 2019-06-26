double do_kahan_sum(double* restrict var, long ncells)
{
   struct esum_type{
      double sum;
      double correction;
   } local;
   local.sum = 0.0;
   local.correction = 0.0;

   for (long i = 0; i < ncells; i++) {
      double corrected_next_term= var[i] + local.correction;
      double new_sum            = local.sum + local.correction;
      local.correction   = corrected_next_term - (new_sum - local.sum);
      local.sum          = new_sum;
   }

   double sum = local.sum + local.correction;
   return(sum);
}

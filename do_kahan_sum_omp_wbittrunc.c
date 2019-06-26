typedef unsigned int uint;
double bittruncate(double sum, uint nbits);

double do_kahan_sum_omp_wbittrunc(double* restrict var, long ncells, uint nbits)
{
   struct esum_type{
      double sum;
      double correction;
   };

   double sum = 0.0;
   double correction = 0.0;

#pragma omp parallel reduction(+:sum, correction)
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

//    sum += local.correction;
//    sum += local.sum;
         correction = local.correction;
         corrected_next_term = sum + correction;
         new_sum = sum + correction;
         correction = corrected_next_term - (new_sum - sum);
         sum = new_sum;
#ifdef _OPENMP
#pragma omp barrier
#endif
         correction = local.sum;
         corrected_next_term = sum + correction;
         new_sum = sum + correction;
         correction = corrected_next_term - (new_sum - sum);
         sum = new_sum;
   }

   sum = bittruncate(sum, nbits);

   return(sum);
}

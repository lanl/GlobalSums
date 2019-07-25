double do_sum_novec(double* restrict var, long ncells)
{
   // Serial sum
   double sum = 0.0;
   for (long i = 0; i < ncells; i++){
      sum += var[i];
   }

   return(sum);
}

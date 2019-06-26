double digitround(double var, int ndigits);

double do_sum_wdigittrunc(double *var, long ncells, int ndigits)
{
   // Serial sum
   double sum = 0.0;
   for (long i = 0; i < ncells; i++){
      sum += var[i];
   }

   sum = digitround(sum, ndigits);
   return(sum);
}

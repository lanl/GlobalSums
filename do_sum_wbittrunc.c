typedef unsigned int uint;
double bittruncate(double sum, uint nbits);

double do_sum_wbittrunc(double* restrict var, long ncells, uint nbits)
{
   // Serial sum
   double sum = 0.0;
   for (long i = 0; i < ncells; i++){
      sum += var[i];
   }

   sum = bittruncate(sum, nbits);
   return(sum);
}

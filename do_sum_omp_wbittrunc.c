typedef unsigned int uint;
double bittruncate(double sum, uint nbits);

double do_sum_omp_wbittrunc(double* restrict var, long ncells, uint nbits)
{
   // Serial sum
   double sum = 0.0;
#pragma omp parallel for reduction(+: sum)
   for (long i = 0; i < ncells; i++){
      sum += var[i];
   }

   sum = bittruncate(sum, nbits);
   return(sum);
}

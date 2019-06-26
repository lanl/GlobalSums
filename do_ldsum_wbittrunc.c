typedef unsigned int uint;
double bittruncate(double sum, uint nbits);

long double do_ldsum_wbittrunc(double* restrict var, long ncells, uint nbits)
{
   // Serial sum with long doubles
   long double ldsum = 0.0;
   for (long i = 0; i < ncells; i++){
      ldsum += var[i];
   }

   ldsum = bittruncate(ldsum, nbits);
   return(ldsum);
}

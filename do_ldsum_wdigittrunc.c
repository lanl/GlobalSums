double digitround(double var, int ndigits);

long double do_ldsum_wdigittrunc(double* restrict var, long ncells, int ndigits)
{
   // Serial sum with long doubles
   long double ldsum = 0.0;
   for (long i = 0; i < ncells; i++){
      ldsum += var[i];
   }

   ldsum = digitround(ldsum, ndigits);
   return(ldsum);
}

long double do_ldsum(double* restrict var, long ncells)
{
   // Serial sum with long doubles
   long double ldsum = 0.0;
   for (long i = 0; i < ncells; i++){
      ldsum += var[i];
   }
   return(ldsum);
}

#include <math.h>
#include <quadmath.h>

__float128 do_qdsum_wtrunc(double* restrict var, long ncells, int ndigits)
{
   // Serial sum with quad doubles
   __float128 qdsum = 0.0;
   for (long i = 0; i < ncells; i++){
      qdsum += (__float128)var[i];
   }

   int n = (int)log10((double)qdsum);
   __float128 mult = pow((double)10.0,(double)(ndigits-n));

   qdsum = round(qdsum*mult)/mult;
   return(qdsum);
}

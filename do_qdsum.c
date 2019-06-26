#include <quadmath.h>

__float128 do_qdsum(double *var, long ncells)
{
   // Serial sum with quad doubles
   __float128 qdsum = 0.0;
   for (long i = 0; i < ncells; i++){
      qdsum += (__float128)var[i];
   }

   return(qdsum);
}

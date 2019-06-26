#include <quadmath.h>

__float128 do_full_qdsum(__float128* restrict varq, long ncells)
{
   // Serial sum with quad doubles
   __float128 qdsum = 0.0;
   for (long i = 0; i < ncells; i++){
      qdsum += (__float128)varq[i];
   }
   return(qdsum);
}

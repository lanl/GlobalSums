#include <x86intrin.h>

int main(int argc, char *argv[])
{
   __m512d local_sum = {0.0};

   double var[8] = {0.0};
   __m512d var_v = _mm512_load_pd(var);
}

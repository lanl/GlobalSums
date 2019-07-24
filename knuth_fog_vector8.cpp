#include "vectorclass.h"

static double sum[8] __attribute__ ((aligned (64)));

extern "C" {
double do_knuth_sum_agner_v8(double* var, long ncells);
}

double do_knuth_sum_agner_v8(double* var, long ncells)
{
   Vec8d local_sum(0.0);
   Vec8d local_correction(0.0);
   Vec8d var_v;

   for (long i = 0; i < ncells; i+=8) {
      var_v.load(var+i);
      Vec8d u = local_sum;
      Vec8d v = var_v + local_correction;
      Vec8d upt = u + v;
      Vec8d up = upt - v;
      Vec8d vpp = upt - up;
      local_sum = upt;
      local_correction = (u - up) + (v - vpp);
   }

   Vec8d sum_v = local_sum + local_correction;
   sum_v.store(sum);

   // double to do final sum
   double ud, vd, uptd, upd, vppd;
   double final_sum = 0.0;
   double final_correction = 0.0;

   for (long i = 0; i < 8; i++) {
      ud = final_sum;
      vd = sum_v[i] + final_correction;
      uptd = ud + vd;
      upd = uptd - vd;
      vppd = uptd - upd;
      final_sum = uptd;
      final_correction = (ud - upd) + (vd - vppd);
   }

   return(final_sum);
}

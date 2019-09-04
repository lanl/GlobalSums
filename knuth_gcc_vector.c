static double sum[4] __attribute__ ((aligned (64)));

double do_knuth_sum_gcc_v(double* restrict var, long ncells)
{
   typedef double vec4d __attribute__ ((vector_size(4 * sizeof(double))));

   vec4d local_sum = {0.0};
   vec4d local_correction = {0.0};
   vec4d var_v;

   for (long i = 0; i < ncells; i+=4) {
      var_v = *(vec4d *)&var[i];
      vec4d u = local_sum;
      vec4d v = var_v + local_correction;
      vec4d upt = u + v;
      vec4d up = upt - v;
      vec4d vpp = upt - up;
      local_sum = upt;
      local_correction = (u - up) + (v - vpp);
   }

   vec4d sum_v = local_sum + local_correction;
   *(vec4d *)sum = sum_v;

   // double to do final sum
   double ud, vd, uptd, upd, vppd;
   double final_sum = 0.0;
   double final_correction = 0.0;

   for (long i = 0; i < 4; i++) {
      ud = final_sum;
      vd = sum[i] + final_correction;
      uptd = ud + vd;
      upd = uptd - vd;
      vppd = uptd - upd;
      final_sum = uptd;
      final_correction = (ud - upd) + (vd - vppd);
   }

   return(final_sum);
}

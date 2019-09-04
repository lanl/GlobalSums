static double sum[8] __attribute__ ((aligned (64)));

double do_knuth_sum_gcc_v8(double* restrict var, long ncells)
{
    typedef double vec8d __attribute__ ((vector_size(8 * sizeof(double))));

   vec8d local_sum = {0.0};
   vec8d local_correction = {0.0};
   vec8d var_v;

   for (long i = 0; i < ncells; i+=8) {
      var_v = *(vec8d *)&var[i];
      vec8d u = local_sum;
      vec8d v = var_v + local_correction;
      vec8d upt = u + v;
      vec8d up = upt - v;
      vec8d vpp = upt - up;
      local_sum = upt;
      local_correction = (u - up) + (v - vpp);
   }

   vec8d sum_v = local_sum + local_correction;
   *(vec8d *)sum = sum_v;

   // double to do final sum
   double ud, vd, uptd, upd, vppd;
   double final_sum = 0.0;
   double final_correction = 0.0;

   for (long i = 0; i < 8; i++) {
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

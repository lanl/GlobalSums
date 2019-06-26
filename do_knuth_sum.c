double do_knuth_sum(double* restrict var, long ncells)
{
   struct esum_type{
      double sum;
      double correction;
   };

   double u, v, upt, up, vpp;
   struct esum_type local;

   local.sum = 0.0;
   local.correction = 0.0;
   for (long i = 0; i < ncells; i++) {
      u = local.sum;
      v = var[i] + local.correction;
      upt = u + v;
      up = upt - v;
      vpp = upt - up;
      local.sum = upt;
      local.correction = (u - up) + (v - vpp);
   }

   double sum = local.sum + local.correction;
   return(sum);
}

/*
 *  Copyright (c) 2015, Los Alamos National Security, LLC.
 *  All rights Reserved.
 *
 *  Distributed under the OSI Certified Apache License 2.0
 *
 *  GlobalSums, Version 1.0.0 (C16001) -- LA-CC-15-087
 *
 *  Author -- Bob Robey, brobey@lanl.gov
 *
 *  ABSTRACT
 *  A demonstration code to support a paper Computational Reproducibility for
 *  Production Physics Applications submitted to the Numerical Reproducibility
 *  at Exascale (NRE 2015) workshop at the 2015 Supercomputing conference, Nov 20, 2015.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include <sys/time.h>
#include <immintrin.h>
#include <x86intrin.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define DEBUG 0
#define ORDERS_OF_MAGNITUDE 1.0e9;
#define QORDERS_OF_MAGNITUDE 1.0e9q;

#define MIN(a,b) (((a)<(b))?(a):(b))

typedef unsigned int uint;

double do_sum(double* restrict var, long ncells);
double do_sum_wdigittrunc(double* restrict var, long ncells, int ndigits);
double do_sum_wbittrunc(double* restrict var, long ncells, uint nbits);
long double do_ldsum(double* restrict var, long ncells);
long double do_ldsum_wdigittrunc(double* restrict var, long ncells, int ndigits);
long double do_ldsum_wbittrunc(double* restrict var, long ncells, uint nbits);
double do_kahan_sum(double* restrict var, long ncells);
double do_kahan_sum_v(double* restrict var, long ncells);
double do_kahan_sum_gcc_v(double* restrict var, long ncells);
double do_knuth_sum(double* restrict var, long ncells);
double do_knuth_sum_v(double* restrict var, long ncells);

void do_kahan_sum_omp(double *var, long ncells, double accurate_sum);
void do_kahan_sum_omp_wbittrunc(double *var, long ncells, double accurate_sum, uint nbits);
void do_pair_sum(double *var, long ncells, double accurate_sum);

void do_sum_omp(double *var, long ncells, double accurate_sum);
void do_sum_omp_wbittrunc(double *var, long ncells, double accurate_sum, uint nbits);

void do_qdsum(double *var, long ncells, __float128 accurate_qdsum);
void do_qdsum_wtrunc(double *var, long ncells, __float128 accurate_qdsum, int ndigits);
void do_full_qdsum(__float128 *var, long ncells, __float128 accurate_qdsum);
void do_full_qdsum_wtrunc(__float128 *var, long ncells, __float128 accurate_qdsum, int ndigits);

void cpu_timer_start(struct timeval *tstart_cpu);
double cpu_timer_stop(struct timeval tstart_cpu);

double digitround(double var, int ndigits);
double bittruncate(double sum, uint nbits);

int main(int argc, char *argv[])
{

#ifdef _OPENMP
   int nt = 0;
   int tid = 0;

   nt = omp_get_max_threads();
   tid = omp_get_thread_num();
   if (0 == tid) {
        printf("--- max num openmp threads: %d\n", nt);
   }
#pragma omp parallel
   {
      nt = omp_get_num_threads();
      tid = omp_get_thread_num();

#pragma omp master
      printf("--- num openmp threads in parallel region: %d\n", nt);
   }
#endif

  for (int pow_of_two = 4; pow_of_two < 31; pow_of_two++){
      long ncells = (long)pow((double)2,(double)pow_of_two);
      long ncellsdiv2 = ncells/2;

      uint nbits = 4+(uint)(30.0 * (log(ncells)/log(1073741824)));
      uint nbitsld = (uint)(18.0 * (log(ncells)/log(1073741824)));
      uint nbitsomp = 2+(uint)(28.0 * (log(ncells)/log(1073741824)));
      uint nbitskahan = 2;

      int ndigits = 3+(int)(6.0 * (log(ncells)/log(1073741824)));
      int ndigitsld = (int)(6.0 * (log(ncells)/log(1073741824)));

      printf("SETTINGS INFO -- ncells %ld log %d ndigits %d ndigitsld %d nbits %d nbitsld %d\n",ncells,(int)log2((double)ncells),ndigits,ndigitsld,nbits,nbitsld);
   
      double high_value = 1.0e-1;
      double low_value  = 1.0e-1/ORDERS_OF_MAGNITUDE;
      double accurate_sum = (double)ncellsdiv2 * high_value +
                            (double)ncellsdiv2 * low_value;

      long double accurate_ldsum = (long double)ncellsdiv2 * (long double)high_value +
                                   (long double)ncellsdiv2 * (long double)low_value;

      __float128 high_valueq = 1.0e-1q;
      __float128 low_valueq  = 1.0e-1q/QORDERS_OF_MAGNITUDE;
      __float128 accurate_qdsum = (__float128)ncellsdiv2 * high_valueq +
                                  (__float128)ncellsdiv2 * low_valueq;
   
      double *energy = (double *)malloc(ncells*sizeof(double));

      // Initialize with high values first
      printf("Initializing mesh with Leblanc problem, high values first\n");
      for (long i = 0; i < ncells; i++){
         energy[i] = (i < ncellsdiv2) ? high_value : low_value;
      }

      double test_sum, test_accurate_sum;
      long double test_ldsum, test_accurate_ldsum;
      struct timeval cpu_timer;
      double cpu_time;

//******************************************************

      cpu_timer_start(&cpu_timer);

      test_sum = do_sum(energy, ncells);

      cpu_time = cpu_timer_stop(cpu_timer);
      printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
             accurate_sum,test_sum,test_sum-accurate_sum,(test_sum-accurate_sum)/accurate_sum, cpu_time);
      printf("   Serial sum\n");

//******************************************************

      cpu_timer_start(&cpu_timer);

      test_sum = do_sum_wdigittrunc(energy, ncells, ndigits);

      cpu_time = cpu_timer_stop(cpu_timer);
      test_accurate_sum = digitround(accurate_sum, ndigits);
      printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
             test_accurate_sum,test_sum,test_sum-test_accurate_sum,(test_sum-test_accurate_sum)/test_accurate_sum, cpu_time);
      printf("   Serial sum with digit truncation\n");

//******************************************************

      cpu_timer_start(&cpu_timer);

      test_sum = do_sum_wbittrunc(energy, ncells, nbits);

      cpu_time = cpu_timer_stop(cpu_timer);

      test_accurate_sum = bittruncate(accurate_sum, nbits);
      printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
             test_accurate_sum,test_sum,test_sum-test_accurate_sum,(test_sum-test_accurate_sum)/test_accurate_sum, cpu_time);
      printf("   Serial sum with bit truncation\n");

//******************************************************

      cpu_timer_start(&cpu_timer);

      test_ldsum = do_ldsum(energy, ncells);

      cpu_time = cpu_timer_stop(cpu_timer);
      printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
             (double)accurate_ldsum,(double)test_ldsum,(double)(test_ldsum-accurate_ldsum),(double)((test_ldsum-accurate_ldsum)/accurate_ldsum), cpu_time);
      printf("   Serial sum with long double accumulator\n");

//******************************************************

      cpu_timer_start(&cpu_timer);

      test_ldsum = do_ldsum_wdigittrunc(energy, ncells, ndigitsld);

      test_accurate_ldsum = digitround(accurate_ldsum, ndigitsld);

      cpu_time = cpu_timer_stop(cpu_timer);
      printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
             (double)test_accurate_ldsum,(double)test_ldsum,(double)(test_ldsum-test_accurate_ldsum),(double)((test_ldsum-test_accurate_ldsum)/test_accurate_ldsum), cpu_time);
      printf("   Serial sum with long double accumulator with ndigit truncation\n");

//******************************************************

      cpu_timer_start(&cpu_timer);

      test_ldsum = do_ldsum_wbittrunc(energy, ncells, nbitsld);

      test_accurate_ldsum = bittruncate(accurate_ldsum, nbitsld);

      cpu_time = cpu_timer_stop(cpu_timer);
      printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
             (double)test_accurate_ldsum,(double)test_ldsum,(double)(test_ldsum-test_accurate_ldsum),(double)((test_ldsum-test_accurate_ldsum)/test_accurate_ldsum), cpu_time);
      printf("   Serial sum with long double accumulator with bit truncation\n");

//******************************************************

      cpu_timer_start(&cpu_timer);

      test_sum = do_kahan_sum(energy, ncells);

      cpu_time = cpu_timer_stop(cpu_timer);
      printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
             accurate_sum,test_sum,(test_sum-accurate_sum),((test_sum-accurate_sum)/accurate_sum), cpu_time);
      printf("   Serial sum with double double kahan sum accumulator\n");

//******************************************************

      cpu_timer_start(&cpu_timer);

      test_sum = do_kahan_sum_v(energy, ncells);

      cpu_time = cpu_timer_stop(cpu_timer);
      printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
          accurate_sum,test_sum,(test_sum-accurate_sum),((test_sum-accurate_sum)/accurate_sum), cpu_time);
      printf("   Vectorized sum with double double kahan sum accumulator\n");

//******************************************************

      cpu_timer_start(&cpu_timer);

      test_sum = do_kahan_sum_gcc_v(energy, ncells);

      cpu_time = cpu_timer_stop(cpu_timer);
      printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
             accurate_sum,test_sum,(test_sum-accurate_sum),((test_sum-accurate_sum)/accurate_sum), cpu_time);
      printf("   GCC Extensions Vectorized sum with double double kahan sum accumulator\n");

//******************************************************

      cpu_timer_start(&cpu_timer);

      test_sum = do_knuth_sum(energy, ncells);

      cpu_time = cpu_timer_stop(cpu_timer);
      printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
             accurate_sum,test_sum,(test_sum-accurate_sum),((test_sum-accurate_sum)/accurate_sum), cpu_time);
      printf("   Serial sum with double double knuth sum accumulator\n");

//******************************************************

      cpu_timer_start(&cpu_timer);

      test_sum = do_knuth_sum_v(energy, ncells);

      cpu_time = cpu_timer_stop(cpu_timer);
      printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
             accurate_sum,test_sum,(test_sum-accurate_sum),((test_sum-accurate_sum)/accurate_sum), cpu_time);
      printf("   Vectorized sum with double double knuth sum accumulator\n");

//******************************************************

      do_pair_sum(energy, ncells, accurate_sum);

//******************************************************

      printf("\n");

//******************************************************

      do_qdsum(energy, ncells, accurate_qdsum);

//******************************************************

//    do_qdsum_wtrunc(energy, ncells, accurate_qdsum, 17);

//******************************************************

      free(energy);

      __float128 *energyq = (__float128 *)malloc(ncells*sizeof(__float128));

      for (long i = 0; i < ncells; i++){
         energyq[i] = (i < ncellsdiv2) ? high_valueq : low_valueq;
      }

//******************************************************

      do_full_qdsum(energyq, ncells, accurate_qdsum);

//******************************************************

//    do_full_qdsum_wtrunc(energyq, ncells, accurate_qdsum, 28);

//******************************************************

      free(energyq);

      printf("\n");

      energy = (double *)malloc(ncells*sizeof(double));

      // Initialize with high values first
#pragma omp parallel for
      for (long i = 0; i < ncells; i++){
         energy[i] = (i < ncellsdiv2) ? high_value : low_value;
      }

//******************************************************

      do_sum_omp(energy, ncells, accurate_sum);

//******************************************************

      do_sum_omp_wbittrunc(energy, ncells, accurate_sum, nbitsomp);

//******************************************************

      do_kahan_sum_omp(energy, ncells, accurate_sum);

//******************************************************

      do_kahan_sum_omp_wbittrunc(energy, ncells, accurate_sum, nbitskahan);

//******************************************************

      free(energy);

      printf("\n");
   }
}

   
void do_sum_omp(double *var, long ncells, double accurate_sum)
{
   struct timeval cpu_timer;

   cpu_timer_start(&cpu_timer);

   // Serial sum
   double sum = 0.0;
#pragma omp parallel for reduction(+: sum)
   for (long i = 0; i < ncells; i++){
      sum += var[i];
   }

   double cpu_time = cpu_timer_stop(cpu_timer);

   printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
          accurate_sum,sum,sum-accurate_sum,(sum-accurate_sum)/accurate_sum, cpu_time);
   printf("   OpenMP sum\n");
}

void do_sum_omp_wbittrunc(double *var, long ncells, double accurate_sum, uint nbits)
{
   struct timeval cpu_timer;

   cpu_timer_start(&cpu_timer);

   // Serial sum
   double sum = 0.0;
#pragma omp parallel for reduction(+: sum)
   for (long i = 0; i < ncells; i++){
      sum += var[i];
   }

   sum = bittruncate(sum, nbits);
   accurate_sum = bittruncate(accurate_sum, nbits);

   double cpu_time = cpu_timer_stop(cpu_timer);

   printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
          accurate_sum,sum,sum-accurate_sum,(sum-accurate_sum)/accurate_sum, cpu_time);
   printf("   OpenMP sum with bit truncation\n");
}

void do_kahan_sum_omp(double *var, long ncells, double accurate_sum)
{
   struct timeval cpu_timer;

   cpu_timer_start(&cpu_timer);

   struct esum_type{
      double sum;
      double correction;
   };

   double sum = 0.0;

#pragma omp parallel reduction(+:sum)
   {
      double corrected_next_term, new_sum;
      struct esum_type local;

      local.sum = 0.0;
      local.correction = 0.0;
#pragma omp for
      for (long i = 0; i < ncells; i++) {
         corrected_next_term= var[i] + local.correction;
         new_sum      = local.sum + local.correction;
         local.correction   = corrected_next_term - (new_sum - local.sum);
         local.sum          = new_sum;
      }

      sum += local.correction;
#pragma omp barrier
      sum += local.sum;
   }

   double cpu_time = cpu_timer_stop(cpu_timer);
   
   printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
          accurate_sum,sum,(sum-accurate_sum),((sum-accurate_sum)/accurate_sum), cpu_time);
   printf("   OpenMP sum with double double kahan sum accumulator\n");
}

void do_kahan_sum_omp_wbittrunc(double *var, long ncells, double accurate_sum, uint nbits)
{
   struct timeval cpu_timer;

   cpu_timer_start(&cpu_timer);

   struct esum_type{
      double sum;
      double correction;
   };

   double sum = 0.0;
   double correction = 0.0;

#pragma omp parallel reduction(+:sum, correction)
   {
      double corrected_next_term, new_sum;
      struct esum_type local;

      local.sum = 0.0;
      local.correction = 0.0;
#pragma omp for
      for (long i = 0; i < ncells; i++) {
         corrected_next_term= var[i] + local.correction;
         new_sum      = local.sum + local.correction;
         local.correction   = corrected_next_term - (new_sum - local.sum);
         local.sum          = new_sum;
      }

//    sum += local.correction;
//    sum += local.sum;
         correction = local.correction;
         corrected_next_term = sum + correction;
         new_sum = sum + correction;
         correction = corrected_next_term - (new_sum - sum);
         sum = new_sum;
#ifdef _OPENMP
#pragma omp barrier
#endif
         correction = local.sum;
         corrected_next_term = sum + correction;
         new_sum = sum + correction;
         correction = corrected_next_term - (new_sum - sum);
         sum = new_sum;
   }

   sum = bittruncate(sum, nbits);
   accurate_sum = bittruncate(accurate_sum, nbits);

   double cpu_time = cpu_timer_stop(cpu_timer);
   
   printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
          accurate_sum,sum,(sum-accurate_sum),((sum-accurate_sum)/accurate_sum), cpu_time);
   printf("   OpenMP sum with double double kahan sum accumulator with bit truncation\n");
}

void do_pair_sum(double *var, long ncells, double accurate_sum)
{
   struct timeval cpu_timer;

   cpu_timer_start(&cpu_timer);

   // Pair-wise sum
   double *pwsum = (double *)malloc(ncells/2*sizeof(double));

   long nmax = ncells/2;
   for (long i = 0; i<nmax; i++){
      pwsum[i] = var[i*2]+var[i*2+1];
   }

   for (long j = 1; j<log2(ncells); j++){
      nmax /= 2;
      for (long i = 0; i<nmax; i++){
         pwsum[i] = pwsum[i*2]+pwsum[i*2+1];
      }
   }
   double sum = pwsum[0];

   free(pwsum);

   double cpu_time = cpu_timer_stop(cpu_timer);
   
   printf("  accurate sum %-17.16lg sum %-17.16lg diff %10.4lg relative diff %10.4lg runtime %lf",
          accurate_sum,sum,sum-accurate_sum,(sum-accurate_sum)/accurate_sum, cpu_time);
   printf("   Pair-wise sum\n");

}

void do_qdsum(double *var, long ncells, __float128 accurate_qdsum)
{
   struct timeval cpu_timer;

   cpu_timer_start(&cpu_timer);

   char quadstring1[40], quadstring2[40], quadstring3[40], quadstring4[40];

   // Serial sum with quad doubles
   __float128 qdsum = 0.0;
   for (long i = 0; i < ncells; i++){
      qdsum += (__float128)var[i];
   }

   double cpu_time = cpu_timer_stop(cpu_timer);
   
   quadmath_snprintf(quadstring1,24,"%-25.24Qg",accurate_qdsum);
   quadmath_snprintf(quadstring2,24,"%-25.24Qg",qdsum);
   quadmath_snprintf(quadstring3,24,"%-20.14Qg",qdsum-accurate_qdsum);
   quadmath_snprintf(quadstring4,24,"%-20.14Qg",(qdsum-accurate_qdsum)/accurate_qdsum);
   printf("  accurate sum %-24s sum %-24s diff %-20s relative diff %-20s runtime %lf",
          quadstring1,quadstring2,quadstring3,quadstring4,cpu_time);
   printf("   Serial sum with quad double accumulator\n");
}

void do_qdsum_wtrunc(double *var, long ncells, __float128 accurate_qdsum, int ndigits)
{
   struct timeval cpu_timer;

   cpu_timer_start(&cpu_timer);

   char quadstring1[40], quadstring2[40], quadstring3[40], quadstring4[40];

   // Serial sum with quad doubles
   __float128 qdsum = 0.0;
   for (long i = 0; i < ncells; i++){
      qdsum += (__float128)var[i];
   }

   int n = (int)log10((double)qdsum);
   __float128 mult = pow((double)10.0,(double)(ndigits-n));

   qdsum = round(qdsum*mult)/mult;
   accurate_qdsum = round(accurate_qdsum*mult)/mult;

   double cpu_time = cpu_timer_stop(cpu_timer);
   
   quadmath_snprintf(quadstring1,24,"%-25.24Qg",accurate_qdsum);
   quadmath_snprintf(quadstring2,24,"%-25.24Qg",qdsum);
   quadmath_snprintf(quadstring3,24,"%-20.14Qg",qdsum-accurate_qdsum);
   quadmath_snprintf(quadstring4,24,"%-20.14Qg",(qdsum-accurate_qdsum)/accurate_qdsum);
   printf("  accurate sum %-24s sum %-24s diff %-20s relative diff %-20s runtime %lf",
          quadstring1,quadstring2,quadstring3,quadstring4,cpu_time);
   printf("   Serial sum with quad double accumulator with truncation\n");
}

void do_full_qdsum(__float128 *varq, long ncells, __float128 accurate_qdsum)
{
   struct timeval cpu_timer;

   cpu_timer_start(&cpu_timer);

   char quadstring1[40], quadstring2[40], quadstring3[40], quadstring4[40];

   // Serial sum with quad doubles
   __float128 qdsum = 0.0;
   for (long i = 0; i < ncells; i++){
      qdsum += (__float128)varq[i];
   }

   double cpu_time = cpu_timer_stop(cpu_timer);
   
   quadmath_snprintf(quadstring1,24,"%-25.24Qg",accurate_qdsum);
   quadmath_snprintf(quadstring2,24,"%-25.24Qg",qdsum);
   quadmath_snprintf(quadstring3,24,"%-20.14Qg",qdsum-accurate_qdsum);
   quadmath_snprintf(quadstring4,24,"%-20.14Qg",(qdsum-accurate_qdsum)/accurate_qdsum);
   printf("  accurate sum %-24s sum %-24s diff %-20s relative diff %-20s runtime %lf",
          quadstring1,quadstring2,quadstring3,quadstring4,cpu_time);
   printf("   Serial sum with quad double accumulator and quad terms\n");
}

void do_full_qdsum_wtrunc(__float128 *varq, long ncells, __float128 accurate_qdsum, int ndigits)
{
   struct timeval cpu_timer;

   cpu_timer_start(&cpu_timer);

   char quadstring1[40], quadstring2[40], quadstring3[40], quadstring4[40];

   // Serial sum with quad doubles
   __float128 qdsum = 0.0;
   for (long i = 0; i < ncells; i++){
      qdsum += (__float128)varq[i];
   }

   int n = (int)log10((double)qdsum);
   __float128 mult = pow((double)10.0,ndigits-n);

   qdsum = round(qdsum*mult)/mult;
   accurate_qdsum = round(accurate_qdsum*mult)/mult;

   double cpu_time = cpu_timer_stop(cpu_timer);
   
   quadmath_snprintf(quadstring1,24,"%-25.24Qg",accurate_qdsum);
   quadmath_snprintf(quadstring2,24,"%-25.24Qg",qdsum);
   quadmath_snprintf(quadstring3,24,"%-20.14Qg",qdsum-accurate_qdsum);
   quadmath_snprintf(quadstring4,24,"%-20.14Qg",(qdsum-accurate_qdsum)/accurate_qdsum);
   printf("  accurate sum %-24s sum %-24s diff %-20s relative diff %-20s runtime %lf",
          quadstring1,quadstring2,quadstring3,quadstring4,cpu_time);
   printf("   Serial sum with quad double accumulator and quad terms with truncation\n");
}

void cpu_timer_start(struct timeval *tstart_cpu){
   gettimeofday(tstart_cpu, NULL);
}

double cpu_timer_stop(struct timeval tstart_cpu){
   double result;
   struct timeval tstop_cpu, tresult;

   gettimeofday(&tstop_cpu, NULL);
   tresult.tv_sec = tstop_cpu.tv_sec - tstart_cpu.tv_sec;
   tresult.tv_usec = tstop_cpu.tv_usec - tstart_cpu.tv_usec;
   result = (double)tresult.tv_sec + (double)tresult.tv_usec*1.0e-6;
   return(result);
}

double digitround(double var, int ndigits)
{
   int n = (int)log10(var);
   int nshift = 15 - ndigits - n;
   if (nshift >= 0) {
      double mult = pow((double)10.0,nshift);
      return(round(var*mult)/mult);
   } else {
      double div = pow((double)10.0,abs(nshift));
      return(round(var/div)*div);
   }
}

double bittruncate(double var, uint nbits)
{
   unsigned long long bitmask[41] = {
      0x00000000,
      0x00000001, 0x00000003, 0x00000007, 0x0000000F,
      0x0000001F, 0x0000003F, 0x0000007F, 0x000000FF,
      0x000001FF, 0x000003FF, 0x000007FF, 0x00000FFF,
      0x00001FFF, 0x00003FFF, 0x00007FFF, 0x0000FFFF,
      0x0001FFFF, 0x0003FFFF, 0x0007FFFF, 0x000FFFFF,
      0x001FFFFF, 0x003FFFFF, 0x007FFFFF, 0x00FFFFFF,
      0x01FFFFFF, 0x03FFFFFF, 0x07FFFFFF, 0x0FFFFFFF,
      0x1FFFFFFF, 0x3FFFFFFF, 0x7FFFFFFF, 0xFFFFFFFF,
      0x1FFFFFFFF, 0x3FFFFFFFF, 0x7FFFFFFFF, 0xFFFFFFFFF,
      0x1FFFFFFFFF, 0x3FFFFFFFFF, 0x7FFFFFFFFF, 0xFFFFFFFFFF
   };

   union twiddler {
      double dvalue;
      unsigned long long ivalue;
   } q;

   nbits = MIN(40,nbits);

   q.dvalue = var;
   q.ivalue &= ~bitmask[nbits];
   var = q.dvalue;

   return(var);
}

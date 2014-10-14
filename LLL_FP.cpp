
#include <NTL/LLL.h>
#include <iostream>
#include <NTL/fileio.h>
#include <NTL/vec_double.h>


#include <NTL/new.h>

NTL_START_IMPL

static inline 
void CheckFinite(double *p)
{
   if (!IsFinite(p)) Error("LLL_FP: numbers too big...use LLL_XD");
}

static double InnerProduct(double *a, double *b, long n)
{
   double s;
   long i;

   s = 0;
   for (i = 1; i <= n; i++) 
      s += a[i]*b[i];

   return s;
}

static void RowTransform(vec_ZZ& A, vec_ZZ& B, const ZZ& MU1)
// x = x - y*MU
{
   static ZZ T, MU;
   long k;

   long n = A.length();
   long i;

   MU = MU1;

   if (MU == 1) {
      for (i = 1; i <= n; i++)
         sub(A(i), A(i), B(i));

      return;
   }

   if (MU == -1) {
      for (i = 1; i <= n; i++)
         add(A(i), A(i), B(i));

      return;
   }

   if (MU == 0) return;

   if (NumTwos(MU) >= NTL_ZZ_NBITS) 
      k = MakeOdd(MU);
   else
      k = 0;


   if (MU.WideSinglePrecision()) {
      long mu1;
      conv(mu1, MU);

      if (k > 0) {

         for (i = 1; i <= n; i++) {
            mul(T, B(i), mu1);
            LeftShift(T, T, k);
            sub(A(i), A(i), T);
         }

      }
      else {

         for (i = 1; i <= n; i++) {
            MulSubFrom(A(i), B(i), mu1);
         }

      }
   }
   else {
      for (i = 1; i <= n; i++) {
         mul(T, B(i), MU);
         if (k > 0) LeftShift(T, T, k);
         sub(A(i), A(i), T);
      }
   }
}


#define TR_BND (NTL_FDOUBLE_PRECISION/2.0)
// Just to be safe!!

static double max_abs(double *v, long n)
{
   long i;
   double res, t;

   res = 0;

   for (i = 1; i <= n; i++) {
      t = fabs(v[i]);
      if (t > res) res = t;
   }

   return res;
}


static void RowTransformStart(double *a, long *in_a, long& in_float, long n)
{
   long i;
   long inf = 1;

   for (i = 1; i <= n; i++) {
      in_a[i] = (a[i] < TR_BND && a[i] > -TR_BND);
      inf = inf & in_a[i];
   }

   in_float = inf;
}


static void RowTransformFinish(vec_ZZ& A, double *a, long *in_a)
{
   long n = A.length();
   long i;

   for (i = 1; i <= n; i++) {
      if (in_a[i])  {
         conv(A(i), a[i]);
      }
      else {
         conv(a[i], A(i));
         CheckFinite(&a[i]);
      }
   }
}


static void RowTransform(vec_ZZ& A, vec_ZZ& B, const ZZ& MU1, 
                         double *a, double *b, long *in_a,
                         double& max_a, double max_b, long& in_float)
// x = x - y*MU
{
   static ZZ T, MU;
   long k;
   double mu;

   conv(mu, MU1);
   CheckFinite(&mu);

   long n = A.length();
   long i;

   if (in_float) {
      double mu_abs = fabs(mu);
      if (mu_abs > 0 && max_b > 0 && (mu_abs >= TR_BND || max_b >= TR_BND)) {
         in_float = 0;
      }
      else {
         max_a += mu_abs*max_b;
         if (max_a >= TR_BND) 
            in_float = 0;
      }
   }

   if (in_float) {
      if (mu == 1) {
         for (i = 1; i <= n; i++)
            a[i] -= b[i];

         return;
      }

      if (mu == -1) {
         for (i = 1; i <= n; i++)
            a[i] += b[i];

         return;
      }

      if (mu == 0) return;

      for (i = 1; i <= n; i++)
         a[i] -= mu*b[i];


      return;
   }


   MU = MU1;

   if (MU == 1) {
      for (i = 1; i <= n; i++) {
         if (in_a[i] && a[i] < TR_BND && a[i] > -TR_BND &&
             b[i] < TR_BND && b[i] > -TR_BND) {

            a[i] -= b[i];
         }
         else {
            if (in_a[i]) {
               conv(A(i), a[i]);
               in_a[i] = 0;
            }
         
            sub(A(i), A(i), B(i));
         }
      }
      return;
   }

   if (MU == -1) {
      for (i = 1; i <= n; i++) {
         if (in_a[i] && a[i] < TR_BND && a[i] > -TR_BND &&
             b[i] < TR_BND && b[i] > -TR_BND) {

            a[i] += b[i];
         }
         else {
            if (in_a[i]) {
               conv(A(i), a[i]);
               in_a[i] = 0;
            }
         
            add(A(i), A(i), B(i));
         }
      }
      return;
   }

   if (MU == 0) return;

   double b_bnd = fabs(TR_BND/mu) - 1;
   if (b_bnd < 0) b_bnd = 0; 

   if (NumTwos(MU) >= NTL_ZZ_NBITS) 
      k = MakeOdd(MU);
   else
      k = 0;


   if (MU.WideSinglePrecision()) {
      long mu1;
      conv(mu1, MU);

      if (k > 0) {
         for (i = 1; i <= n; i++) {
            if (in_a[i]) {
               conv(A(i), a[i]);
               in_a[i] = 0;
            }

            mul(T, B(i), mu1);
            LeftShift(T, T, k);
            sub(A(i), A(i), T);
         }
      }
      else {
         for (i = 1; i <= n; i++) {
            if (in_a[i] && a[i] < TR_BND && a[i] > -TR_BND &&
                b[i] < b_bnd && b[i] > -b_bnd) {
   
               a[i] -= b[i]*mu;
            }
            else {
               if (in_a[i]) {
                  conv(A(i), a[i]);
                  in_a[i] = 0;
               }
               MulSubFrom(A(i), B(i), mu1);
            }
         }
      }
   }
   else {
      for (i = 1; i <= n; i++) {
         if (in_a[i]) {
            conv(A(i), a[i]);
            in_a[i] = 0;
         }
         mul(T, B(i), MU);
         if (k > 0) LeftShift(T, T, k);
         sub(A(i), A(i), T);
      }
   }
}

static void RowTransform2(vec_ZZ& A, vec_ZZ& B, const ZZ& MU1)
// x = x + y*MU

{
   static ZZ T, MU;
   long k;

   long n = A.length();
   long i;

   MU = MU1;

   if (MU == 1) {
      for (i = 1; i <= n; i++)
         add(A(i), A(i), B(i));

      return;
   }

   if (MU == -1) {
      for (i = 1; i <= n; i++)
         sub(A(i), A(i), B(i));

      return;
   }

   if (MU == 0) return;

   if (NumTwos(MU) >= NTL_ZZ_NBITS) 
      k = MakeOdd(MU);
   else
      k = 0;

   if (MU.WideSinglePrecision()) {
      long mu1;
      conv(mu1, MU);

      for (i = 1; i <= n; i++) {
         mul(T, B(i), mu1);
         if (k > 0) LeftShift(T, T, k);
         add(A(i), A(i), T);
      }
   }
   else {
      for (i = 1; i <= n; i++) {
         mul(T, B(i), MU);
         if (k > 0) LeftShift(T, T, k);
         add(A(i), A(i), T);
      }
   }
}

static
void ComputeGS(mat_ZZ& B, double **B1, double **mu, double *b,
               double *c, long k, double bound, long st, double *buf)

{
   long n = B.NumCols();
   long i, j;
   double s, t1, y, t;

   ZZ T1;
   long test;

   double *mu_k = mu[k];

   if (st < k) {
      for (i = 1; i < st; i++)
         buf[i] = mu_k[i]*c[i];
   }

   for (j = st; j <= k-1; j++) {
      s = InnerProduct(B1[k], B1[j], n);

      // test = b[k]*b[j] >= NTL_FDOUBLE_PRECISION^2

      test = (b[k]/NTL_FDOUBLE_PRECISION >= NTL_FDOUBLE_PRECISION/b[j]);

      // test = test && s^2 <= b[k]*b[j]/bound,
      // but we compute it in a strange way to avoid overflow

      if (test && (y = fabs(s)) != 0) {
         t = y/b[j];
         t1 = b[k]/y;
         if (t <= 1)
            test = (t*bound <= t1);
         else if (t1 >= 1)
            test = (t <= t1/bound);
         else
            test = 0;
      }

      if (test) {
         InnerProduct(T1, B(k), B(j));
         conv(s, T1);
      }

      double *mu_j = mu[j];

      t1 = 0;
      for (i = 1; i <= j-1; i++) {
         t1 += mu_j[i]*buf[i];
      }
  
      mu_k[j] = (buf[j] = (s - t1))/c[j];
   }

#if (!NTL_EXT_DOUBLE)

   // Kahan summation 

   double c1;

   s = c1 = 0;
   for (j = 1; j <= k-1; j++) {
      y = mu_k[j]*buf[j] - c1;
      t = s+y;
      c1 = t-s;
      c1 = c1-y;
      s = t;
   }


#else

   s = 0;
   for (j = 1; j <= k-1; j++)
      s += mu_k[j]*buf[j];

#endif

   c[k] = b[k] - s;
}

static double red_fudge = 0;
static long log_red = 0;

static long verbose = 0;

double LLLStatusInterval = 900.0;
char *LLLDumpFile = 0;

static unsigned long NumSwaps = 0;
static double RR_GS_time = 0;
static double StartTime = 0;
static double LastTime = 0;



static void LLLStatus(long max_k, double t, long m, const mat_ZZ& B)
{
   cerr << "---- LLL_FP status ----\n";
   cerr << "elapsed time: ";
   PrintTime(cerr, t-StartTime);
   cerr << ", stage: " << max_k;
   cerr << ", rank: " << m;
   cerr << ", swaps: " << NumSwaps << "\n";

   ZZ t1;
   long i;
   double prodlen = 0;

   for (i = 1; i <= m; i++) {
      InnerProduct(t1, B(i), B(i));
      if (!IsZero(t1))
         prodlen += log(t1);
   }

   cerr << "log of prod of lengths: " << prodlen/(2.0*log(2.0)) << "\n";

   if (LLLDumpFile) {
      cerr << "dumping to " << LLLDumpFile << "...";

      ofstream f;
      OpenWrite(f, LLLDumpFile);
      
      f << "[";
      for (i = 1; i <= m; i++) {
         f << B(i) << "\n";
      }
      f << "]\n";

      f.close();

      cerr << "\n";
   }

   LastTime = t;
   
}

static void init_red_fudge()
{
   long i;

   log_red = long(0.50*NTL_DOUBLE_PRECISION);
   red_fudge = 1;

   for (i = log_red; i > 0; i--)
      red_fudge = red_fudge*0.5;
}

static void inc_red_fudge()
{

   red_fudge = red_fudge * 2;
   log_red--;

   
   cerr << "LLL_FP: warning--relaxing reduction (" << log_red << ")\n";

   if (log_red < 4)
      Error("LLL_FP: too much loss of precision...stop!");
}


#if 0

static void print_mus(double **mu, long k)
{
   long i;

   for (i = k-1; i >= 1; i--)
      cerr << mu[k][i] << " ";
   cerr << "\n";
}

#endif

void ComputeGS(const mat_ZZ& B, mat_RR& B1, 
               mat_RR& mu, vec_RR& b,
               vec_RR& c, long k, const RR& bound, long st,
               vec_RR& buf, const RR& bound2);



static void RR_GS(mat_ZZ& B, double **B1, double **mu, 
                  double *b, double *c, double *buf, long prec,
                  long rr_st, long k, long m_orig,
                  mat_RR& rr_B1, mat_RR& rr_mu, 
                  vec_RR& rr_b, vec_RR& rr_c)
{
   double tt;

   cerr << "LLL_FP: RR refresh " << rr_st << "..." << k << "...";
   tt = GetTime();

   if (rr_st > k) Error("LLL_FP: can not continue!!!");

   long old_p = RR::precision();
   RR::SetPrecision(prec);

   long n = B.NumCols();

   rr_B1.SetDims(k, n);
   rr_mu.SetDims(k, m_orig);
   rr_b.SetLength(k);
   rr_c.SetLength(k);

   vec_RR rr_buf;
   rr_buf.SetLength(k);

   long i, j;

   for (i = rr_st; i <= k; i++)
      for (j = 1; j <= n; j++)
         conv(rr_B1(i, j), B(i, j));

   for (i = rr_st; i <= k; i++)
      InnerProduct(rr_b(i), rr_B1(i), rr_B1(i));

   

   RR bound;
   power2(bound, 2*long(0.15*RR::precision()));

   RR bound2;
   power2(bound2, 2*RR::precision());

   for (i = rr_st; i <= k; i++)
      ComputeGS(B, rr_B1, rr_mu, rr_b, rr_c, i, bound, 1, rr_buf, bound2);

   for (i = rr_st; i <= k; i++)
      for (j = 1; j <= n; j++) {
         conv(B1[i][j], rr_B1(i,j));
         CheckFinite(&B1[i][j]);
      }

   for (i = rr_st; i <= k; i++)
      for (j = 1; j <= i-1; j++) {
         conv(mu[i][j], rr_mu(i,j));
      }

   for (i = rr_st; i <= k; i++) {
      conv(b[i], rr_b(i));
      CheckFinite(&b[i]);
   }
   

   for (i = rr_st; i <= k; i++) {
      conv(c[i], rr_c(i));
      CheckFinite(&c[i]);
   }

   for (i = 1; i <= k-1; i++) {
      conv(buf[i], rr_buf[i]);
   }


   RR::SetPrecision(old_p);

   tt = GetTime()-tt;
   RR_GS_time += tt;
   cerr << tt << " (" << RR_GS_time << ")\n";
}

void ComputeGS(const mat_ZZ& B, mat_RR& mu, vec_RR& c)
{
   long n = B.NumCols();
   long k = B.NumRows();

   mat_RR B1;
   vec_RR b;

   B1.SetDims(k, n);
   mu.SetDims(k, k);
   b.SetLength(k);
   c.SetLength(k);

   vec_RR buf;
   buf.SetLength(k);

   long i, j;

   for (i = 1; i <= k; i++)
      for (j = 1; j <= n; j++)
         conv(B1(i, j), B(i, j));

   for (i = 1; i <= k; i++)
      InnerProduct(b(i), B1(i), B1(i));

   

   RR bound;
   power2(bound, 2*long(0.15*RR::precision()));

   RR bound2;
   power2(bound2, 2*RR::precision());


   for (i = 1; i <= k; i++)
      ComputeGS(B, B1, mu, b, c, i, bound, 1, buf, bound2);

}






static
long ll_LLL_FP(mat_ZZ& B, mat_ZZ* U, double delta, long deep, 
           LLLCheckFct check, double **B1, double **mu, 
           double *b, double *c,
           long m, long init_k, long &quit)
{
   long n = B.NumCols();

   long i, j, k, Fc1;
   ZZ MU;
   double mu1;

   double t1;
   ZZ T1;
   double *tp;


   static double bound = 0;

   if (bound == 0) {
      // we tolerate a 15% loss of precision in computing
      // inner products in ComputeGS.

      bound = 1;
      for (i = 2*long(0.15*NTL_DOUBLE_PRECISION); i > 0; i--)
         bound = bound * 2;
   }

   double half_plus_fudge = 0.5 + red_fudge;

   quit = 0;
   k = init_k;


   vec_long st_mem;
   st_mem.SetLength(m+2);
   long *st = st_mem.elts();

   for (i = 1; i < k; i++)
      st[i] = i;

   for (i = k; i <= m+1; i++)
      st[i] = 1;

   double *buf;
   buf = NTL_NEW_OP double [m+1];
   if (!buf) Error("out of memory in lll_LLL_FP");

   vec_long in_vec_mem;
   in_vec_mem.SetLength(n+1);
   long *in_vec = in_vec_mem.elts();

   double *max_b;
   max_b = NTL_NEW_OP double [m+1];
   if (!max_b) Error("out of memory in lll_LLL_FP");

   for (i = 1; i <= m; i++)
      max_b[i] = max_abs(B1[i], n);

   long in_float;

   long rst;
   long counter;
   long start_over;

   long trigger_index;
   long small_trigger;
   long cnt;

   mat_RR rr_B1;
   mat_RR rr_mu;
   vec_RR rr_c;
   vec_RR rr_b;

   long m_orig = m;

   long rr_st = 1;

   long max_k = 0;

   long prec = RR::precision();

   double tt;

   long swap_cnt = 0;


   while (k <= m) {

      if (k > max_k) {
         max_k = k;
         swap_cnt = 0;
      }

      if (verbose) {
         tt = GetTime();

         if (tt > LastTime + LLLStatusInterval)
            LLLStatus(max_k, tt, m, B);
      }

      if (k < rr_st) rr_st = k;

      if (st[k] == k)
         rst = 1;
      else
         rst = k;

      if (st[k] < st[k+1]) st[k+1] = st[k];
      ComputeGS(B, B1, mu, b, c, k, bound, st[k], buf);
      CheckFinite(&c[k]);
      st[k] = k;

      if (swap_cnt > 200000) {
         cerr << "LLL_FP: swap loop?\n";
         RR_GS(B, B1, mu, b, c, buf, prec,
               rr_st, k, m_orig, rr_B1, rr_mu, rr_b, rr_c);
         if (rr_st < st[k+1]) st[k+1] = rr_st;
         rr_st = k+1;
         rst = k;
         swap_cnt = 0;
      }

      counter = 0;
      trigger_index = k;
      small_trigger = 0;
      cnt = 0;

      long thresh = 10;
      long sz=0, new_sz;

      long did_rr_gs = 0;


      do {
         // size reduction

         counter++;
         if ((counter & 127) == 0) {

            new_sz = 0;
            for (j = 1; j <= n; j++)
               new_sz += NumBits(B(k,j));

            if ((counter >> 7) == 1 || new_sz < sz) {
               sz = new_sz;
            }
            else {
               cerr << "LLL_FP: warning--infinite loop?\n";
            }
         }

         Fc1 = 0;
         start_over = 0;
   
         for (j = rst-1; j >= 1; j--) {
            t1 = fabs(mu[k][j]);
            if (t1 > half_plus_fudge) { 


               if (!Fc1) {
                  if (j > trigger_index || 
                      (j == trigger_index && small_trigger)) {

                     cnt++;

                     if (cnt > thresh) {
                        if (log_red <= 15) { 

                           while (log_red > 10)
                              inc_red_fudge();

                           half_plus_fudge = 0.5 + red_fudge;

                           if (!did_rr_gs) {
                              RR_GS(B, B1, mu, b, c, buf, prec,
                                    rr_st, k, m_orig, rr_B1, rr_mu, rr_b, rr_c);
                              if (rr_st < st[k+1]) st[k+1] = rr_st;
                              rr_st = k+1;
                              did_rr_gs = 1;
                              rst = k;
                              trigger_index = k;
                              small_trigger = 0;
                              start_over = 1;
                              break;
                           }
                        }
                        else {
                           inc_red_fudge();
                           half_plus_fudge = 0.5 + red_fudge;
                           cnt = 0;
                        }
                     }
                  }

                  trigger_index = j;
                  small_trigger = (t1 < 4);

                  Fc1 = 1;
                  if (k < rr_st) rr_st = k;
                  RowTransformStart(B1[k], in_vec, in_float, n);
               }
                  

               mu1 = mu[k][j];
               if (mu1 >= 0)
                  mu1 = ceil(mu1-0.5);
               else
                  mu1 = floor(mu1+0.5);
   
               double *mu_k = mu[k];
               double *mu_j = mu[j];
   
               if (mu1 == 1) {
                  for (i = 1; i <= j-1; i++)
                     mu_k[i] -= mu_j[i];
               }
               else if (mu1 == -1) {
                  for (i = 1; i <= j-1; i++)
                     mu_k[i] += mu_j[i];
               }
               else {
                  for (i = 1; i <= j-1; i++)
                     mu_k[i] -= mu1*mu_j[i];
               }
   
               mu_k[j] -= mu1;
   
               conv(MU, mu1);

               RowTransform(B(k), B(j), MU, B1[k], B1[j], in_vec,
                            max_b[k], max_b[j], in_float);
               if (U) RowTransform((*U)(k), (*U)(j), MU);
            }
         }


         if (Fc1) {
            RowTransformFinish(B(k), B1[k], in_vec);
            max_b[k] = max_abs(B1[k], n);

            if (!did_rr_gs) {
               b[k] = InnerProduct(B1[k], B1[k], n);
               CheckFinite(&b[k]);

               ComputeGS(B, B1, mu, b, c, k, bound, 1, buf);
               CheckFinite(&c[k]);
            }
            else {
               RR_GS(B, B1, mu, b, c, buf, prec,
                     rr_st, k, m_orig, rr_B1, rr_mu, rr_b, rr_c);
               rr_st = k+1;
            }

            rst = k;
         }
      } while (Fc1 || start_over);

      if (check && (*check)(B(k))) 
         quit = 1;

      if (b[k] == 0) {
         for (i = k; i < m; i++) {
            // swap i, i+1
            swap(B(i), B(i+1));
            tp = B1[i]; B1[i] = B1[i+1]; B1[i+1] = tp;
            t1 = b[i]; b[i] = b[i+1]; b[i+1] = t1;
            t1 = max_b[i]; max_b[i] = max_b[i+1]; max_b[i+1] = t1;
            if (U) swap((*U)(i), (*U)(i+1));
         }

         for (i = k; i <= m+1; i++) st[i] = 1;
         if (k < rr_st) rr_st = k;

         m--;
         if (quit) break;
         continue;
      }

      if (quit) break;

      if (deep > 0) {
         // deep insertions

         double cc = b[k];
         long l = 1;
         while (l <= k-1 && delta*c[l] <= cc) {
            cc = cc - mu[k][l]*mu[k][l]*c[l];
            l++;
         }
   
         if (l <= k-1 && (l <= deep || k-l <= deep)) {
            // deep insertion at position l
   
            for (i = k; i > l; i--) {
               // swap rows i, i-1
               swap(B(i), B(i-1));
               tp = B1[i]; B1[i] = B1[i-1]; B1[i-1] = tp;
               tp = mu[i]; mu[i] = mu[i-1]; mu[i-1] = tp;
               t1 = b[i]; b[i] = b[i-1]; b[i-1] = t1;
               t1 = max_b[i]; max_b[i] = max_b[i-1]; max_b[i-1] = t1;
               if (U) swap((*U)(i), (*U)(i-1));
            }
   
            k = l;
            NumSwaps++;
            swap_cnt++;
            continue;
         }
      } // end deep insertions

      // test LLL reduction condition

      if (k > 1 && delta*c[k-1] > c[k] + mu[k][k-1]*mu[k][k-1]*c[k-1]) {
         // swap rows k, k-1
         swap(B(k), B(k-1));
         tp = B1[k]; B1[k] = B1[k-1]; B1[k-1] = tp;
         tp = mu[k]; mu[k] = mu[k-1]; mu[k-1] = tp;
         t1 = b[k]; b[k] = b[k-1]; b[k-1] = t1;
         t1 = max_b[k]; max_b[k] = max_b[k-1]; max_b[k-1] = t1;
         if (U) swap((*U)(k), (*U)(k-1));

         k--;
         NumSwaps++;
         swap_cnt++;
         // cout << "-\n";
      }
      else {

         k++;
         // cout << "+\n";
      }

   }

   if (verbose) {
      LLLStatus(m+1, GetTime(), m, B);
   }


   delete [] buf;
   delete [] max_b;

   return m;
}




static
long LLL_FP(mat_ZZ& B, mat_ZZ* U, double delta, long deep, 
           LLLCheckFct check)
{
   long m = B.NumRows();
   long n = B.NumCols();

   long i, j;
   long new_m, dep, quit;
   ZZ MU;

   ZZ T1;

   init_red_fudge();

   if (U) ident(*U, m);

   double **B1;  // approximates B

   typedef double *doubleptr;

   B1 = NTL_NEW_OP doubleptr[m+1];
   if (!B1) Error("LLL_FP: out of memory");

   for (i = 1; i <= m; i++) {
      B1[i] = NTL_NEW_OP double[n+1];
      if (!B1[i]) Error("LLL_FP: out of memory");
   }

   double **mu;
   mu = NTL_NEW_OP doubleptr[m+1];
   if (!mu) Error("LLL_FP: out of memory");

   for (i = 1; i <= m; i++) {
      mu[i] = NTL_NEW_OP double[m+1];
      if (!mu[i]) Error("LLL_FP: out of memory");
   }

   double *c; // squared lengths of Gramm-Schmidt basis vectors

   c = NTL_NEW_OP double[m+1];
   if (!c) Error("LLL_FP: out of memory");

   double *b; // squared lengths of basis vectors

   b = NTL_NEW_OP double[m+1];
   if (!b) Error("LLL_FP: out of memory");


   for (i = 1; i <=m; i++)
      for (j = 1; j <= n; j++) {
         conv(B1[i][j], B(i, j));
         CheckFinite(&B1[i][j]);
      }

         
   for (i = 1; i <= m; i++) {
      b[i] = InnerProduct(B1[i], B1[i], n);
      CheckFinite(&b[i]);
   }

   new_m = ll_LLL_FP(B, U, delta, deep, check, B1, mu, b, c, m, 1, quit);
   dep = m - new_m;
   m = new_m;

   if (dep > 0) {
      // for consistency, we move all of the zero rows to the front

      for (i = 0; i < m; i++) {
         swap(B(m+dep-i), B(m-i));
         if (U) swap((*U)(m+dep-i), (*U)(m-i));
      }
   }


   // clean-up

   for (i = 1; i <= m+dep; i++) {
      delete [] B1[i];
   }

   delete [] B1;

   for (i = 1; i <= m+dep; i++) {
      delete [] mu[i];
   }

   delete [] mu;

   delete [] c;

   delete [] b;

   return m;
}

         

long LLL_FP(mat_ZZ& B, double delta, long deep, LLLCheckFct check, 
           long verb)
{
   verbose = verb;
   RR_GS_time = 0;
   NumSwaps = 0;
   if (verbose) {
      StartTime = GetTime();
      LastTime = StartTime;
   }

   if (delta < 0.50 || delta >= 1) Error("LLL_FP: bad delta");
   if (deep < 0) Error("LLL_FP: bad deep");
   return LLL_FP(B, 0, delta, deep, check);
}

long LLL_FP(mat_ZZ& B, mat_ZZ& U, double delta, long deep, 
           LLLCheckFct check, long verb)
{
   verbose = verb;
   RR_GS_time = 0;
   NumSwaps = 0;
   if (verbose) {
      StartTime = GetTime();
      LastTime = StartTime;
   }

   if (delta < 0.50 || delta >= 1) Error("LLL_FP: bad delta");
   if (deep < 0) Error("LLL_FP: bad deep");
   return LLL_FP(B, &U, delta, deep, check);
}



static vec_double BKZConstant;

static
void ComputeBKZConstant(long beta, long p)
{
   const double c_PI = 3.14159265358979323846264338328;
   const double LogPI = 1.14472988584940017414342735135;

   BKZConstant.SetLength(beta-1);

   vec_double Log;
   Log.SetLength(beta);


   long i, j, k;
   double x, y;

   for (j = 1; j <= beta; j++)
      Log(j) = log(double(j));

   for (i = 1; i <= beta-1; i++) {
      // First, we compute x = gamma(i/2)^{2/i}

      k = i/2;

      if ((i & 1) == 0) { // i even
         x = 0;
         for (j = 1; j <= k; j++)
            x = x + Log(j);
          
         x = x * (1/double(k));

         x = exp(x);
      }
      else { // i odd
         x = 0;
         for (j = k + 2; j <= 2*k + 2; j++)
            x = x + Log(j);

         x = 0.5*LogPI + x - 2*(k+1)*Log(2);

         x = x * (2.0/double(i));

         x = exp(x);
      }

      // Second, we compute y = 2^{2*p/i}

      y = -(2*p/double(i))*Log(2);
      y = exp(y);

      BKZConstant(i) = x*y/c_PI;
   }
}


static vec_double BKZThresh;

static 
void ComputeBKZThresh(double *c, long beta)
{
   BKZThresh.SetLength(beta-1);

   long i;
   double x;

   x = 0;

   for (i = 1; i <= beta-1; i++) {
      x += log(c[i-1]);
      BKZThresh(i) = exp(x/double(i))*BKZConstant(i);
      if (!IsFinite(&BKZThresh(i))) BKZThresh(i) = 0;
   }
}



static 
void BKZStatus(double tt, double enum_time, unsigned long NumIterations, 
               unsigned long NumTrivial, unsigned long NumNonTrivial, 
               unsigned long NumNoOps, long m, 
               const mat_ZZ& B)
{
   cerr << "---- BKZ_FP status ----\n";
   cerr << "elapsed time: ";
   PrintTime(cerr, tt-StartTime);
   cerr << ", enum time: ";
   PrintTime(cerr, enum_time);
   cerr << ", iter: " << NumIterations << "\n";
   cerr << "triv: " << NumTrivial;
   cerr << ", nontriv: " << NumNonTrivial;
   cerr << ", no ops: " << NumNoOps;
   cerr << ", rank: " << m;
   cerr << ", swaps: " << NumSwaps << "\n";



   ZZ t1;
   long i;
   double prodlen = 0;

   for (i = 1; i <= m; i++) {
      InnerProduct(t1, B(i), B(i));
      if (!IsZero(t1))
         prodlen += log(t1);
   }

   cerr << "log of prod of lengths: " << prodlen/(2.0*log(2.0)) << "\n";


   if (LLLDumpFile) {
      cerr << "dumping to " << LLLDumpFile << "...";

      ofstream f;
      OpenWrite(f, LLLDumpFile);
      
      f << "[";
      for (i = 1; i <= m; i++) {
         f << B(i) << "\n";
      }
      f << "]\n";

      f.close();

      cerr << "\n";
   }

   LastTime = tt;
   
}



static
long BKZ_FP(mat_ZZ& BB, mat_ZZ* UU, double delta, 
         long beta, long prune, LLLCheckFct check)
{

   long m = BB.NumRows();
   long n = BB.NumCols();
   long m_orig = m;
   
   long i, j;
   ZZ MU;

   double t1;
   ZZ T1;
   double *tp;

   init_red_fudge();

   mat_ZZ B;
   B = BB;

   B.SetDims(m+1, n);


   double **B1;  // approximates B

   typedef double *doubleptr;

   B1 = NTL_NEW_OP doubleptr[m+2];
   if (!B1) Error("BKZ_FP: out of memory");

   for (i = 1; i <= m+1; i++) {
      B1[i] = NTL_NEW_OP double[n+1];
      if (!B1[i]) Error("BKZ_FP: out of memory");
   }

   double **mu;
   mu = NTL_NEW_OP doubleptr[m+2];
   if (!mu) Error("LLL_FP: out of memory");

   for (i = 1; i <= m+1; i++) {
      mu[i] = NTL_NEW_OP double[m+1];
      if (!mu[i]) Error("BKZ_FP: out of memory");
   }


   double *c; // squared lengths of Gramm-Schmidt basis vectors

   c = NTL_NEW_OP double[m+2];
   if (!c) Error("BKZ_FP: out of memory");

   double *b; // squared lengths of basis vectors

   b = NTL_NEW_OP double[m+2];
   if (!b) Error("BKZ_FP: out of memory");

   double cbar;

   double *ctilda;
   ctilda = NTL_NEW_OP double[m+2];
   if (!ctilda) Error("BKZ_FP: out of memory");

   double *vvec;
   vvec = NTL_NEW_OP double[m+2];
   if (!vvec) Error("BKZ_FP: out of memory");

   double *yvec;
   yvec = NTL_NEW_OP double[m+2];
   if (!yvec) Error("BKZ_FP: out of memory");

   double *uvec;
   uvec = NTL_NEW_OP double[m+2];
   if (!uvec) Error("BKZ_FP: out of memory");

   double *utildavec;
   utildavec = NTL_NEW_OP double[m+2];
   if (!utildavec) Error("BKZ_FP: out of memory");


   long *Deltavec;
   Deltavec = NTL_NEW_OP long[m+2];
   if (!Deltavec) Error("BKZ_FP: out of memory");

   long *deltavec;
   deltavec = NTL_NEW_OP long[m+2];
   if (!deltavec) Error("BKZ_FP: out of memory");

   mat_ZZ Ulocal;
   mat_ZZ *U;

   if (UU) {
      Ulocal.SetDims(m+1, m);
      for (i = 1; i <= m; i++)
         conv(Ulocal(i, i), 1);
      U = &Ulocal;
   }
   else
      U = 0;

   long quit;
   long new_m;
   long z, jj, kk;
   long s, t;
   long h;
   double eta;


   for (i = 1; i <=m; i++)
      for (j = 1; j <= n; j++) {
         conv(B1[i][j], B(i, j));
         CheckFinite(&B1[i][j]);
      }

         
   for (i = 1; i <= m; i++) {
      b[i] = InnerProduct(B1[i], B1[i], n);
      CheckFinite(&b[i]);
   }



   m = ll_LLL_FP(B, U, delta, 0, check, B1, mu, b, c, m, 1, quit);

   double tt;

   double enum_time = 0;
   unsigned long NumIterations = 0;
   unsigned long NumTrivial = 0;
   unsigned long NumNonTrivial = 0;
   unsigned long NumNoOps = 0;

   long verb = verbose;

   verbose = 0;

   long clean = 1;

   if (m < m_orig) {
      for (i = m_orig+1; i >= m+2; i--) {
         // swap i, i-1

         swap(B(i), B(i-1));
         if (U) swap((*U)(i), (*U)(i-1));
      }
   }

   if (!quit && m > 1) {
      if (beta > m) beta = m;

      if (prune > 0) 
         ComputeBKZConstant(beta, prune);

      z = 0;
      jj = 0;
   
      while (z < m-1) {
         jj++;
         kk = min(jj+beta-1, m);
   
         if (jj == m) {
            jj = 1;
            kk = beta;
            clean = 1;
         }

         if (verb) {
            tt = GetTime();
            if (tt > LastTime + LLLStatusInterval)
               BKZStatus(tt, enum_time, NumIterations, NumTrivial,
                         NumNonTrivial, NumNoOps, m, B);
         }

   
         // ENUM

         double tt1;

         if (verb) {
            tt1 = GetTime();
         }


         if (prune > 0)
            ComputeBKZThresh(&c[jj], kk-jj+1);

   
         cbar = c[jj];
         utildavec[jj] = uvec[jj] = 1;
   
         yvec[jj] = vvec[jj] = 0;
         Deltavec[jj] = 0;
   
   
         s = t = jj;
         deltavec[jj] = 1;
   
         for (i = jj+1; i <= kk+1; i++) {
            ctilda[i] = uvec[i] = utildavec[i] = yvec[i] = 0;
            Deltavec[i] = 0;
            vvec[i] = 0;
            deltavec[i] = 1;
         }

         long enum_cnt = 0;
   
         while (t <= kk) {
            if (verb) {
               enum_cnt++;
               if (enum_cnt > 100000) {
                  enum_cnt = 0;
                  tt = GetTime();
                  if (tt > LastTime + LLLStatusInterval) {
                     enum_time += tt - tt1;
                     tt1 = tt;
                     BKZStatus(tt, enum_time, NumIterations, NumTrivial,
                               NumNonTrivial, NumNoOps, m, B);
                  }
               }
            }

            ctilda[t] = ctilda[t+1] + 
               (yvec[t]+utildavec[t])*(yvec[t]+utildavec[t])*c[t];

            ForceToMem(&ctilda[t]);  // prevents an infinite loop
   
            if (prune > 0 && t > jj) {
               eta = BKZThresh(t-jj);
            }
            else
               eta = 0;
   
            if (ctilda[t] < cbar - eta) {
               if (t > jj) {
                  t--;
                  t1 = 0;
                  for (i = t+1; i <= s; i++)
                     t1 += utildavec[i]*mu[i][t];
                  yvec[t] = t1;
                  t1 = -t1;
                  if (t1 >= 0)
                     t1 = ceil(t1-0.5);
                  else
                     t1 = floor(t1+0.5);
                  utildavec[t] = vvec[t] = t1;
                  Deltavec[t] = 0;
                  if (utildavec[t] > -yvec[t]) 
                     deltavec[t] = -1;
                  else
                     deltavec[t] = 1;
               }
               else {
                  cbar = ctilda[jj];
                  for (i = jj; i <= kk; i++) {
                     uvec[i] = utildavec[i];
                  }
               }
            }
            else {
               t++;
               s = max(s, t);
               if (t < s) Deltavec[t] = -Deltavec[t];
               if (Deltavec[t]*deltavec[t] >= 0) Deltavec[t] += deltavec[t];
               utildavec[t] = vvec[t] + Deltavec[t];
            }
         }

         if (verb) {
            tt1 = GetTime() - tt1;
            enum_time += tt1;
         }
         
         NumIterations++;
   
         h = min(kk+1, m);
   
         if ((delta - 8*red_fudge)*c[jj] > cbar) {

            clean = 0;

            // we treat the case that the new vector is b_s (jj < s <= kk)
            // as a special case that appears to occur most of the time.
   
            s = 0;
            for (i = jj+1; i <= kk; i++) {
               if (uvec[i] != 0) {
                  if (s == 0)
                     s = i;
                  else
                     s = -1;
               }
            }
   
            if (s == 0) Error("BKZ_FP: internal error");
   
            if (s > 0) {
               // special case

               NumTrivial++;
   
               for (i = s; i > jj; i--) {
                  // swap i, i-1
                  swap(B(i-1), B(i));
                  if (U) swap((*U)(i-1), (*U)(i));
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }
   
               // cerr << "special case\n";
               new_m = ll_LLL_FP(B, U, delta, 0, check, 
                                B1, mu, b, c, h, jj, quit);
               if (new_m != h) Error("BKZ_FP: internal error");
               if (quit) break;
            }
            else {
               // the general case

               NumNonTrivial++;
   
               for (i = 1; i <= n; i++) conv(B(m+1, i), 0);

               if (U) {
                  for (i = 1; i <= m_orig; i++)
                     conv((*U)(m+1, i), 0);
               }

               for (i = jj; i <= kk; i++) {
                  if (uvec[i] == 0) continue;
                  conv(MU, uvec[i]);
                  RowTransform2(B(m+1), B(i), MU);
                  if (U) RowTransform2((*U)(m+1), (*U)(i), MU);
               }
      
               for (i = m+1; i >= jj+1; i--) {
                  // swap i, i-1
                  swap(B(i-1), B(i));
                  if (U) swap((*U)(i-1), (*U)(i));
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }
      
               for (i = 1; i <= n; i++) {
                  conv(B1[jj][i], B(jj, i));
                  CheckFinite(&B1[jj][i]);
               }
      
               b[jj] = InnerProduct(B1[jj], B1[jj], n);
               CheckFinite(&b[jj]);
      
               if (b[jj] == 0) Error("BKZ_FP: internal error"); 
      
               // remove linear dependencies
   
               // cerr << "general case\n";
               new_m = ll_LLL_FP(B, U, delta, 0, 0, B1, mu, b, c, kk+1, jj, quit);
              
               if (new_m != kk) Error("BKZ_FP: internal error"); 

               // remove zero vector
      
               for (i = kk+2; i <= m+1; i++) {
                  // swap i, i-1
                  swap(B(i-1), B(i));
                  if (U) swap((*U)(i-1), (*U)(i));
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }
      
               quit = 0;
               if (check) {
                  for (i = 1; i <= kk; i++)
                     if ((*check)(B(i))) {
                        quit = 1;
                        break;
                     }
               }

               if (quit) break;
   
               if (h > kk) {
                  // extend reduced basis
   
                  new_m = ll_LLL_FP(B, U, delta, 0, check, 
                                   B1, mu, b, c, h, h, quit);
   
                  if (new_m != h) Error("BKZ_FP: internal error");
                  if (quit) break;
               }
            }
   
            z = 0;
         }
         else {
            // LLL_FP
            // cerr << "progress\n";

            NumNoOps++;

            if (!clean) {
               new_m = 
                  ll_LLL_FP(B, U, delta, 0, check, B1, mu, b, c, h, h, quit);
               if (new_m != h) Error("BKZ_FP: internal error");
               if (quit) break;
            }
   
            z++;
         }
      }
   }


   if (verb) {
      BKZStatus(GetTime(), enum_time, NumIterations, NumTrivial, NumNonTrivial, 
                NumNoOps, m, B);
   }

   // clean up


   if (m_orig > m) {
      // for consistency, we move zero vectors to the front

      for (i = m+1; i <= m_orig; i++) {
         swap(B(i), B(i+1));
         if (U) swap((*U)(i), (*U)(i+1));
      }

      for (i = 0; i < m; i++) {
         swap(B(m_orig-i), B(m-i));
         if (U) swap((*U)(m_orig-i), (*U)(m-i));
      }
   }

   B.SetDims(m_orig, n);
   BB = B;

   if (U) {
      U->SetDims(m_orig, m_orig);
      *UU = *U;
   }

   for (i = 1; i <= m_orig+1; i++) {
      delete [] B1[i];
   }

   delete [] B1;

   for (i = 1; i <= m_orig+1; i++) {
      delete [] mu[i];
   }

   delete [] mu;

   delete [] c;
   delete [] b;
   delete [] ctilda;
   delete [] vvec;
   delete [] yvec;
   delete [] uvec;
   delete [] utildavec;
   delete [] Deltavec;
   delete [] deltavec;

   return m;
}




const double c_PI = 3.14159265358979323846264338328;

static
double logGamma(double x)
{
	 double tmp = (x - 0.5) * log(x + 4.5) - (x + 4.5);  // -z = tmp
	 double ser = 1.0  + 76.18009173    / (x + 0)   - 86.50532033    / (x + 1)
                       + 24.01409822    / (x + 2)   -  1.231739516   / (x + 3)
                       +  0.00120858003 / (x + 4)   -  0.00000536382 / (x + 5); 
	 // see http://en.wikipedia.org/wiki/Gamma_function
	 return tmp + log(ser * sqrt(2 * c_PI)); //s = ser 
}

static
double gamma(double x) { return exp(logGamma(x)); }

static
double UnitVolume(int dim)
{
	double n = dim;
	//return sqrt(c_PI)/gamma(n/2+1);
	return std::pow(c_PI,n/2)/gamma(n/2+1);
	//std::pow(gamma(n/2 +1),1/n)/sqrt(c_PI);
	//return sqrt(c_PI)/std::pow(gamma(n/2 +1),1/n);
}

static
vec_double interpolateBoundingFunction(const vec_double& R, int newDimension)
{
	int i;
	vec_double retVal;
	int R_size = R.length();
	retVal.SetLength(newDimension);
	double ratio = (double)R_size/(double)newDimension;
	double targetIndex;
	int i_1,i_2;
		for(i=1;i<=newDimension;i++)
		{
			targetIndex = i*ratio;
			i_1 = floor(targetIndex+0.0000001);
			i_2 = ceil(targetIndex-0.0000001);
			if (i_1==i_2)
				retVal[i-1]=R[i_1-1];
			else
				retVal[i-1]= R[i_1-1]+((targetIndex-(double)i_1)*(R[i_2-1]-R[i_1-1]));
		}
		// put LAST element =1 ???
		//retVal[newDimension-1]=1;
	return retVal;
}


static double boundingFunctions[18][50]={
//{0.134539,0.150018,0.191214,0.207611,0.245203,0.259457,0.281809,0.292381,0.315642,0.325311,0.352234,0.360752,0.387897,0.400402,0.425706,0.437607,0.463207,0.473418,0.49261,0.502326,0.521107,0.534928,0.548977,0.566764,0.591693,0.597316,0.607629,0.618864,0.632766,0.642476,0.661942,0.666902,0.69268,0.696924,0.715517,0.726498,0.763853,0.771181,0.790907,0.796366,0.81573,0.82723,0.839154,0.848473,0.859563,0.863654,0.899907,0.901125,0.930598,0.93098},
{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
{0.158413,0.172804,0.211677,0.226714,0.260111,0.27964,0.319103,0.329964,0.355499,0.366772,0.387223,0.412774,0.430973,0.44228,0.474222,0.487395,0.513515,0.520666,0.533657,0.544764,0.570794,0.587126,0.598154,0.605092,0.62471,0.63875,0.663992,0.672189,0.690816,0.704527,0.717844,0.727764,0.735513,0.742691,0.768265,0.774202,0.796958,0.802141,0.820619,0.830508,0.852281,0.858906,0.871832,0.876767,0.891519,0.894277,0.916979,0.919513,0.940237,0.940532},
{0.168717,0.192896,0.212165,0.223837,0.269002,0.29436,0.324413,0.335089,0.360159,0.376202,0.398621,0.413805,0.435585,0.442884,0.46233,0.474343,0.512765,0.525039,0.544372,0.560863,0.57871,0.595154,0.601528,0.6142,0.641589,0.648855,0.680304,0.688691,0.704985,0.708833,0.723931,0.732585,0.755866,0.771072,0.784548,0.791462,0.829924,0.845811,0.854695,0.858965,0.871851,0.880376,0.885188,0.887802,0.911756,0.913988,0.918325,0.919486,0.952625,0.952944},
{0.162795,0.172546,0.21972,0.232156,0.263248,0.281081,0.321292,0.343326,0.385825,0.407185,0.422466,0.44503,0.462569,0.471939,0.495201,0.503983,0.51867,0.529763,0.557827,0.570885,0.590724,0.596891,0.625297,0.637895,0.656608,0.665881,0.692311,0.70807,0.727727,0.739999,0.755656,0.76816,0.783623,0.793689,0.809143,0.817057,0.840927,0.851919,0.872609,0.881441,0.891491,0.89606,0.903935,0.906936,0.91628,0.918562,0.934729,0.935649,0.957438,0.957574},
{0.160753,0.171057,0.238764,0.257386,0.275973,0.288633,0.31866,0.328857,0.354133,0.369052,0.405067,0.42271,0.439321,0.450674,0.497827,0.512384,0.534056,0.547787,0.563069,0.580596,0.596256,0.608994,0.622425,0.630646,0.653145,0.666755,0.683478,0.694023,0.713599,0.731693,0.759534,0.767484,0.78864,0.803547,0.824698,0.831788,0.841811,0.851879,0.876969,0.888981,0.895573,0.900428,0.915111,0.920683,0.934029,0.938484,0.943847,0.945464,0.963699,0.963877},
{0.181307,0.205524,0.270077,0.302604,0.327004,0.343813,0.374518,0.394141,0.406366,0.419862,0.447825,0.457671,0.472359,0.482899,0.508067,0.514505,0.545275,0.552956,0.570477,0.586794,0.604676,0.614673,0.624623,0.638259,0.663041,0.669309,0.696355,0.707035,0.732818,0.744551,0.760053,0.767742,0.806987,0.825652,0.838184,0.85067,0.864639,0.873307,0.893896,0.905854,0.916891,0.922505,0.933265,0.936277,0.9404,0.941999,0.949786,0.950804,0.969494,0.969578},
{0.214426,0.252644,0.301593,0.315764,0.335732,0.353957,0.390616,0.409404,0.427716,0.443599,0.475309,0.484126,0.514655,0.528104,0.554715,0.566422,0.577757,0.588242,0.604978,0.613869,0.633516,0.640547,0.678054,0.686453,0.711526,0.722405,0.735872,0.74986,0.772968,0.787552,0.796693,0.802808,0.832201,0.839982,0.854536,0.86729,0.87828,0.885318,0.902765,0.90861,0.91537,0.920798,0.935634,0.940547,0.946108,0.947174,0.948455,0.949012,0.971616,0.971696},
{0.204701,0.240387,0.276414,0.299504,0.334309,0.352551,0.37281,0.389629,0.413331,0.428174,0.440727,0.453433,0.488252,0.504115,0.520134,0.538391,0.55996,0.568876,0.593274,0.599307,0.628964,0.649904,0.663326,0.677531,0.704053,0.722448,0.735008,0.740596,0.75278,0.76171,0.770346,0.777448,0.8181,0.826353,0.834074,0.843326,0.85768,0.867836,0.884881,0.889949,0.912026,0.91899,0.935691,0.938342,0.948112,0.950986,0.958091,0.960158,0.978743,0.978896},
{0.182391,0.199508,0.244234,0.266255,0.293529,0.314042,0.330391,0.346381,0.384306,0.405835,0.442872,0.464816,0.489261,0.501425,0.524397,0.536828,0.557513,0.569291,0.588235,0.602353,0.629467,0.648533,0.661023,0.678428,0.692115,0.702899,0.717018,0.726852,0.746996,0.760181,0.784492,0.797725,0.806128,0.817072,0.842486,0.850411,0.867677,0.877229,0.896666,0.906926,0.919645,0.927463,0.935733,0.94055,0.950234,0.953106,0.96556,0.96751,0.980585,0.981388},
{0.221479,0.246223,0.287889,0.309424,0.3497,0.364573,0.377501,0.390184,0.422914,0.430694,0.454035,0.465466,0.491072,0.515517,0.531574,0.544786,0.570671,0.585461,0.605503,0.614945,0.630411,0.647181,0.668105,0.676661,0.710412,0.721894,0.737622,0.744476,0.764671,0.789081,0.796656,0.808393,0.817955,0.829949,0.839019,0.849478,0.868268,0.875362,0.89404,0.900554,0.920056,0.925856,0.938312,0.942894,0.957746,0.961082,0.9693,0.971201,0.985506,0.985877},
{0.177671,0.205005,0.246867,0.26754,0.29338,0.30976,0.328699,0.341923,0.393551,0.413278,0.432373,0.449221,0.482518,0.491072,0.505419,0.527836,0.544006,0.564458,0.595512,0.613555,0.627067,0.649972,0.664191,0.673047,0.680266,0.708926,0.725191,0.737323,0.754498,0.767752,0.780278,0.793694,0.808352,0.822188,0.831658,0.84146,0.857292,0.868746,0.891272,0.904124,0.914426,0.919861,0.932702,0.938986,0.950711,0.954552,0.965296,0.966747,0.989303,0.989806},
{0.154034,0.193018,0.246192,0.264836,0.300638,0.31625,0.33708,0.367313,0.380351,0.408273,0.424585,0.446907,0.464301,0.483562,0.498847,0.516204,0.535292,0.548498,0.564431,0.591635,0.598879,0.611605,0.643339,0.651438,0.663432,0.672468,0.703531,0.724195,0.735549,0.747452,0.760654,0.772919,0.789906,0.808783,0.837484,0.847901,0.877348,0.885651,0.895365,0.904115,0.91994,0.928246,0.943878,0.950028,0.960767,0.967968,0.973383,0.975453,0.991967,0.992591},
{0.17276,0.188923,0.216311,0.24064,0.269172,0.295842,0.321507,0.360126,0.385203,0.391821,0.415051,0.437722,0.451838,0.487603,0.509601,0.520109,0.555338,0.56291,0.58579,0.59997,0.609548,0.626575,0.645857,0.665574,0.692844,0.70691,0.725508,0.734851,0.747465,0.76057,0.781334,0.79581,0.806456,0.811926,0.836929,0.85605,0.873116,0.878852,0.892641,0.897289,0.915637,0.92292,0.939011,0.946932,0.959548,0.962614,0.982744,0.984086,0.994241,0.994805},
{0.188904,0.232397,0.260438,0.290067,0.314869,0.331418,0.361056,0.382012,0.406197,0.420188,0.44926,0.467564,0.482209,0.495676,0.507137,0.529239,0.55073,0.563029,0.583046,0.608377,0.630584,0.645432,0.66261,0.674199,0.690529,0.705593,0.724674,0.748691,0.757543,0.771951,0.794019,0.807969,0.832595,0.843828,0.858985,0.867436,0.877417,0.885254,0.90424,0.910588,0.923153,0.934634,0.944425,0.956036,0.967457,0.972795,0.980526,0.982416,0.993653,0.994274},
{0.152285,0.187803,0.222114,0.244365,0.303979,0.331484,0.36283,0.38691,0.406782,0.42373,0.438033,0.45676,0.496573,0.518113,0.539617,0.558762,0.577649,0.586643,0.605987,0.619717,0.641366,0.658549,0.665726,0.678271,0.691363,0.703322,0.724296,0.737341,0.751806,0.767468,0.78333,0.794921,0.80803,0.818741,0.849467,0.859481,0.874713,0.888351,0.900336,0.91566,0.929864,0.939213,0.9527,0.962877,0.970394,0.976771,0.984556,0.987076,0.995956,0.996204},
{0.176202,0.203937,0.265673,0.295211,0.322367,0.340329,0.381,0.41542,0.438665,0.460492,0.479147,0.496237,0.517581,0.543305,0.558308,0.572928,0.590301,0.610036,0.626231,0.643547,0.661526,0.672739,0.688525,0.706701,0.726085,0.734609,0.749387,0.759141,0.781332,0.795197,0.804233,0.814887,0.834892,0.853353,0.869237,0.883001,0.897798,0.909861,0.919175,0.926389,0.93808,0.946631,0.954494,0.961779,0.972773,0.977747,0.985679,0.988462,0.995963,0.996544},
{0.174593,0.207913,0.254526,0.296637,0.327998,0.366564,0.394117,0.417932,0.440638,0.462083,0.474244,0.503974,0.524182,0.534963,0.54887,0.556578,0.571783,0.586801,0.600913,0.618281,0.645643,0.657758,0.668543,0.691589,0.709503,0.718954,0.732838,0.754332,0.776896,0.788715,0.79965,0.816726,0.82759,0.849075,0.868893,0.885648,0.895641,0.901746,0.921308,0.930002,0.942744,0.95202,0.964283,0.968641,0.976894,0.982095,0.99184,0.994183,0.998253,0.998767},
{0.197048,0.23833,0.261482,0.266641,0.290854,0.318289,0.362154,0.390196,0.430542,0.467162,0.487569,0.504246,0.525159,0.54563,0.580928,0.58751,0.618524,0.638077,0.646086,0.67159,0.68387,0.694005,0.711823,0.729891,0.738479,0.750235,0.763228,0.778772,0.787494,0.806435,0.82693,0.835249,0.852136,0.858957,0.86698,0.876982,0.888717,0.903496,0.915532,0.926721,0.949551,0.958758,0.967987,0.971912,0.986628,0.989587,0.992259,0.994535,0.998671,0.999565},
};


ZZ InnerProduct(const vec_ZZ& x, const vec_ZZ& y)
{
	ZZ retVal;
	InnerProduct(retVal, x,y);
	return retVal;
}


ZZ Norm(const vec_ZZ& u)
{
	return InnerProduct(u,u);
}



static
	long BKZ_FP_v2(mat_ZZ& BB, mat_ZZ* UU, long alfa,  long beta, double gamma, double delta, 
				   double prunningProbability, long vNumRetry, BKZCheckFct check)
{

   long m = BB.NumRows();
   long n = BB.NumCols();
   long m_orig = m;
   
   long i, j;
   ZZ MU;

   double t1;
   ZZ T1;
   double *tp;

   init_red_fudge();

   mat_ZZ B;
   mat_ZZ gsoVectors; // keeps the gso VECTORS !
   B = BB;

   B.SetDims(m+1, n);
   gsoVectors.SetDims(m+1,n);


   double **B1;  // approximates B

   typedef double *doubleptr;

   B1 = NTL_NEW_OP doubleptr[m+2];
   if (!B1) Error("BKZ_FP: out of memory");

   double **mu;
   mu = NTL_NEW_OP doubleptr[m+2];
   if (!mu) Error("LLL_FP: out of memory");

   double *c; // squared lengths of Gramm-Schmidt basis vectors

   c = NTL_NEW_OP double[m+2];
   if (!c) Error("BKZ_FP: out of memory");

   double *b; // squared lengths of basis vectors

   b = NTL_NEW_OP double[m+2];
   if (!b) Error("BKZ_FP: out of memory");

   double cbar;

   double *ctilda;
   ctilda = NTL_NEW_OP double[m+2];
   if (!ctilda) Error("BKZ_FP: out of memory");

   double *vvec;
   vvec = NTL_NEW_OP double[m+2];
   if (!vvec) Error("BKZ_FP: out of memory");

   double *yvec;
   yvec = NTL_NEW_OP double[m+2];
   if (!yvec) Error("BKZ_FP: out of memory");

   double *uvec;
   uvec = NTL_NEW_OP double[m+2];
   if (!uvec) Error("BKZ_FP: out of memory");

   double *utildavec;
   utildavec = NTL_NEW_OP double[m+2];
   if (!utildavec) Error("BKZ_FP: out of memory");


   long *Deltavec;
   Deltavec = NTL_NEW_OP long[m+2];
   if (!Deltavec) Error("BKZ_FP: out of memory");

   long *deltavec;
   deltavec = NTL_NEW_OP long[m+2];
   if (!deltavec) Error("BKZ_FP: out of memory");

   double **sigma; // optimizer !
   sigma = NTL_NEW_OP doubleptr[m+2];
   if (!sigma) Error("LLL_FP: out of memory");

   long *r; // optimizer !!
   r = NTL_NEW_OP long[m+2];
   if (!r) Error("BKZ_FP: out of memory");
	
   
   for (i = 1; i <= m+1; i++) {
   
	  mu[i] = NTL_NEW_OP double[m+1];
      if (!mu[i]) Error("BKZ_FP: out of memory");
  
      B1[i] = NTL_NEW_OP double[n+1];
      if (!B1[i]) Error("BKZ_FP: out of memory");

	  sigma[i] = NTL_NEW_OP double[m+1];
	  if (!sigma[i]) Error("BKZ_FP: out of memory");
   }

   // get the already computed bounding function for the given probability !!!
   vec_double R_BoundingFunction,initialBoundingFunction;
   int idxBoundingFunction = ((int)(prunningProbability*100))/5;
   if (idxBoundingFunction>17) // -> we only have up to 0.9 probability 
	   idxBoundingFunction = 17;
   initialBoundingFunction.SetLength(50);
   for(i=0;i<50;i++)
	   initialBoundingFunction[i]=boundingFunctions[idxBoundingFunction][i];
   //cout << initialBoundingFunction << endl;
   
   mat_ZZ Ulocal;
   mat_ZZ *U;

   if (UU) {
      Ulocal.SetDims(m+1, m);
      for (i = 1; i <= m; i++)
         conv(Ulocal(i, i), 1);
      U = &Ulocal;
   }
   else
      U = 0;

   long quit;
   long new_m;
   long Z,z, JJ, jj,KK, kk;
   long s, t;
   long h;
   double eta;


   for (i = 1; i <=m; i++)
      for (j = 1; j <= n; j++) {
         conv(B1[i][j], B(i, j));
         CheckFinite(&B1[i][j]);
      }

         
   for (i = 1; i <= m; i++) {
      b[i] = InnerProduct(B1[i], B1[i], n);
      CheckFinite(&b[i]);
   }
   c[0]=0; c[m]=0;

   m = ll_LLL_FP(B, U, delta, 0, 0, B1, mu, b, c, m, 1, quit);



   unsigned long NumIterations = 0;
   unsigned long NumTrivial = 0;
   unsigned long NumNonTrivial = 0;
   unsigned long NumNoOps = 0;

   double g_enum_radius = gamma;

   long clean = 1;

   if (m < m_orig) {
      for (i = m_orig+1; i >= m+2; i--) {
         // swap i, i-1

         swap(B(i), B(i-1));
         if (U) swap((*U)(i), (*U)(i-1));
      }
   }

   if (!quit && m > 1) {
      if (beta > m) beta = m;
	  if (alfa > m) alfa = m;

      Z = 0;
      JJ = 0;
	  bool main_enum = true;

	  int shuffleCount=0;
	  bool needShuffle=false;
	  bool foundOneSol=true;

      while (Z < m-1) 
	  {
		  if (main_enum)
		  {
			  if (JJ ==1)
			  {
				  cout << "Full round passed ! " << endl;
				  /*for(i=0;i<=m;i++)
					   cout << Norm(B[i]) << endl;*/

			  }
			  JJ = (JJ%(m-1))+1;
			  KK = min (JJ+beta-1,m);
			  // find a new optimum alfa !!!!????
			  
			  jj = JJ;
			  kk = min(JJ+alfa-1,KK);
			  z = 0;
			  if (kk<KK) main_enum = false;
		  }
		  else
		  {
			  if (z<beta-1)
			  {
				  jj = (jj%(KK-1))+1; 
				  kk = min(jj+alfa-1,KK); 
			  }
			  else
			  {
				  jj = JJ;kk = KK; main_enum = true;
			  }
		  }
         // ENUM
		  R_BoundingFunction = interpolateBoundingFunction(initialBoundingFunction,kk-jj+1);
		  //cout << R_BoundingFunction << endl;

		  if (kk-jj+1>30)
		  {
			double det;
			det = 1;
			for(j=jj;j<=kk;j++)
				det *=sqrt(c[j]);
			
			cbar = std::pow(det/UnitVolume(kk-jj+1),(2/double(kk-jj+1)));
			//cout << "cbar=" << cbar*g_enum_radius << "   c[jj]=" << c[jj] << endl;
			cbar = std::min(cbar*g_enum_radius,c[jj]);
		  }
		  else
			cbar = c[jj];
		  
		  //for(int enumRetry=0;enumRetry<vNumRetry;enumRetry++)
		  {
			
		    foundOneSol=false;
			double tt1;
          
			utildavec[jj] = uvec[jj] = 1;
   
			 yvec[jj] = vvec[jj] = 0;
			 Deltavec[jj] = 0;
   
			 s = t = jj;
			 deltavec[jj] = 1;
   
			 r[jj]=jj;
			 for (i = jj+1; i <= kk+1; i++) {
				ctilda[i] = uvec[i] = utildavec[i] = yvec[i] = 0;
				Deltavec[i] = 0;
				vvec[i] = 0;
				deltavec[i] = 1;
				r[i]=i;
				for(j=0;j<=m;j++) sigma[i][j]=0;
			 }
		
			long enumIterations = 0;
			 while (t <= kk)
			 {
				ctilda[t] = ctilda[t+1] + 
				   (yvec[t]+utildavec[t])*(yvec[t]+utildavec[t])*c[t];

				ForceToMem(&ctilda[t]);  // prevents an infinite loop

				eta = R_BoundingFunction[kk-t]*cbar; 

				if (eta==0) Error("BKZ_FP: Bounding function error !!");
			
				if (ctilda[t] < eta) {
				   if (t > jj) {
					  t--;
					  r[t-1]=max(r[t-1],r[t]);
					  for(i=r[t];i>=t+1;i--) 
						  sigma[i][t]=sigma[i+1][t]+utildavec[i]*mu[i][t];
					  t1 = sigma[t+1][t];
					  yvec[t] = t1;
					  t1 = -t1;
					  if (t1 >= 0)
						 t1 = ceil(t1-0.5);
					  else
						 t1 = floor(t1+0.5);
					  utildavec[t] = vvec[t] = t1;
					  Deltavec[t] = 0;
					  if (utildavec[t] > -yvec[t]) 
						 deltavec[t] = -1;
					  else
						 deltavec[t] = 1;
				   }
				   else {
					   enumIterations=0;
					   foundOneSol = true;
					  cbar = ctilda[jj];

					  //cout << endl << "New sv found:" << cbar << endl;
					  for (i = jj; i <= kk; i++) {
						 uvec[i] = utildavec[i];
					  }
				   }
				}
				else {
				   t++;
				   r[t-1]=t;
				   s = max(s, t);
				   if (t < s) Deltavec[t] = -Deltavec[t];
				   if (Deltavec[t]*deltavec[t] >= 0) Deltavec[t] += deltavec[t];
				   utildavec[t] = vvec[t] + Deltavec[t];
				}
				enumIterations++;
				/*
				early abort enum -> WHEN ???
				if ((!main_enum)&&(enumIterations>100 000))
							break;
				*/
			 }
			 
			 /*if (foundOneSol)
				break;*/
			 //else
			 {
				 /*
				 TO DO: EXPAND local projected block (kk,jj) HERE !
				 Then apply ll_LLL_FP (lll for local projected block [jj,kk])
				 */
			 }
		 }
         NumIterations++;
   
         h = min(kk+1, m);
   
		 if (foundOneSol/*&&(delta - 8*red_fudge)*c[jj] > cbar*/) 
		 {

			 //cout << "Adding vector " << endl;

            clean = 0;

            // we treat the case that the new vector is b_s (jj < s <= kk)
            // as a special case that appears to occur most of the time.
   
            s = 0;
            for (i = jj+1; i <= kk; i++) {
               if (uvec[i] != 0) {
                  if (s == 0)
                     s = i;
                  else
                     s = -1;
               }
            }
   
            if (s == 0) Error("BKZ_FP: internal error");
   
            if (s > 0) {
               // special case
				/*cout << "Trivial: " << endl;
				 for (i = jj; i <= kk; i++)  cout << uvec[i] << "  ";
				 cout << endl;*/

               NumTrivial++;
   
               for (i = s; i > jj; i--) {
                  // swap i, i-1
                  swap(B(i-1), B(i));
                  if (U) swap((*U)(i-1), (*U)(i));
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }
   
               // cerr << "special case\n";
               new_m = ll_LLL_FP(B, U, delta, 0, 0, 
                                B1, mu, b, c, h, jj, quit);
               if (new_m != h) Error("BKZ_FP: internal error");
               if (quit) break;
            }
            else {
               // the general case

               NumNonTrivial++;
   
               for (i = 1; i <= n; i++) conv(B(m+1, i), 0);

               if (U) {
                  for (i = 1; i <= m_orig; i++)
                     conv((*U)(m+1, i), 0);
               }

               for (i = jj; i <= kk; i++) {
                  if (uvec[i] == 0) continue;
                  conv(MU, uvec[i]);
                  RowTransform2(B(m+1), B(i), MU);
                  if (U) RowTransform2((*U)(m+1), (*U)(i), MU);
               }
      
               for (i = m+1; i >= jj+1; i--) {
                  // swap i, i-1
                  swap(B(i-1), B(i));
                  if (U) swap((*U)(i-1), (*U)(i));
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }
      
               for (i = 1; i <= n; i++) {
                  conv(B1[jj][i], B(jj, i));
                  CheckFinite(&B1[jj][i]);
               }
      
               b[jj] = InnerProduct(B1[jj], B1[jj], n);
               CheckFinite(&b[jj]);
      
               if (b[jj] == 0) Error("BKZ_FP: internal error"); 
      
               // remove linear dependencies
   
               // cerr << "general case\n";
               new_m = ll_LLL_FP(B, U, delta, 0, 0, B1, mu, b, c, kk+1, jj, quit);
              
               if (new_m != kk) Error("BKZ_FP: internal error"); 

               // remove zero vector
      
               for (i = kk+2; i <= m+1; i++) {
                  // swap i, i-1
                  swap(B(i-1), B(i));
                  if (U) swap((*U)(i-1), (*U)(i));
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }
      
               quit = 0;
               if (check) { // callback the checker !!
                  for (i = 1; i <= kk; i++)
                     if ((*check)(B(i))) {
                        quit = 1;
                        break;
                     }
               }

               if (quit) break;
   
               if (h > kk) {
                  // extend reduced basis
   
                  new_m = ll_LLL_FP(B, U, delta, 0, 0, 
                                   B1, mu, b, c, h, h, quit);
   
                  if (new_m != h) Error("BKZ_FP: internal error");
                  if (quit) break;
               }
            }
   
			if (main_enum)
				Z = 0;
			else
				z = 0;
         }
         else {
            // LLL_FP
            // cerr << "progress\n";
			//cout << "Progress " << z << endl << endl;

            NumNoOps++;

            if (!clean) {
               new_m = 
                  ll_LLL_FP(B, U, delta, 0, 0, B1, mu, b, c, h, h, quit);
               if (new_m != h) Error("BKZ_FP: internal error");
               if (quit) break;
            }
			if (main_enum) 
				Z++; 
			else
				z++;
         }
      }
	}
   // clean up


   if (m_orig > m) {
      // for consistency, we move zero vectors to the front

      for (i = m+1; i <= m_orig; i++) {
         swap(B(i), B(i+1));
         if (U) swap((*U)(i), (*U)(i+1));
      }

      for (i = 0; i < m; i++) {
         swap(B(m_orig-i), B(m-i));
         if (U) swap((*U)(m_orig-i), (*U)(m-i));
      }
   }

   B.SetDims(m_orig, n);
   BB = B;

   if (U) {
      U->SetDims(m_orig, m_orig);
      *UU = *U;
   }

   for (i = 1; i <= m_orig+1; i++) {
      delete [] B1[i];
   }

   delete [] B1;

   for (i = 1; i <= m_orig+1; i++) {
      delete [] mu[i];
   }

   delete [] mu;

   delete [] c;
   delete [] b;
   delete [] ctilda;
   delete [] vvec;
   delete [] yvec;
   delete [] uvec;
   delete [] utildavec;
   delete [] Deltavec;
   delete [] deltavec;

   return m;
}

long BKZ_FP(mat_ZZ& BB, mat_ZZ& UU, double delta, 
         long beta, long prune, LLLCheckFct check, long verb)
{
   verbose = verb;
   RR_GS_time = 0;
   NumSwaps = 0;
   if (verbose) {
      StartTime = GetTime();
      LastTime = StartTime;
   }

   if (delta < 0.50 || delta >= 1) Error("BKZ_FP: bad delta");
   if (beta < 2) Error("BKZ_FP: bad block size");

   return BKZ_FP(BB, &UU, delta, beta, prune, check);
}

long BKZ_FP(mat_ZZ& BB, double delta, 
         long beta, long prune, LLLCheckFct check, long verb)
{
   verbose = verb;
   RR_GS_time = 0;
   NumSwaps = 0;
   if (verbose) {
      StartTime = GetTime();
      LastTime = StartTime;
   }

   if (delta < 0.50 || delta >= 1) Error("BKZ_FP: bad delta");
   if (beta < 2) Error("BKZ_FP: bad block size");

   return BKZ_FP(BB, 0, delta, beta, prune, check);
}


long BKZ_FP_v2(mat_ZZ& BB, long alfa,  long beta, double gamma, double delta, 
				   double prunningProbability, long vNumRetry, BKZCheckFct check)
{
   RR_GS_time = 0;
   NumSwaps = 0;
   if (verbose) {
      StartTime = GetTime();
      LastTime = StartTime;
   }

   if (delta < 0.50 || delta >= 1) Error("BKZ_FP: bad delta");
   if (beta < 2) Error("BKZ_FP: bad block size");

   return BKZ_FP_v2(BB, 0,alfa,beta,gamma, delta,prunningProbability,vNumRetry,check);
}




NTL_END_IMPL

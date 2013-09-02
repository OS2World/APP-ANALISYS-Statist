/* This file is part of statist
**
** It is distributed under the GNU General Public License.
** See the file COPYING for details.
**
** (c) 1997 Dirk Melcher
**  Doerper Damm 4
**  49134 Wallenhorst
**  GERMANY
**  Tel. 05407/7636
**  email: Dirk.Melcher@usf.Uni-Osnabrueck.DE
**
***************************************************************/

/* funcs.c for STATIST */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>

#include "statist.h"
#include "funcs.h"
#include "gettext.h"

/* ==================================================================== */


BOOLEAN equal(REAL x, REAL y) {
/*  REAL d=REAL_EPSILON;    */
  REAL d = x/1.0e09;
  if ( (x>=(y-d)) && (x<=(y+d)) ) {
    return TRUE;
  }
  else {
    return FALSE;
  }
}

/* ==================================================================== */

int  real_compar_down(const void *x, const void *y) {
   if (*((REAL*) x) <  *((REAL*) y)) {
      return  1;
   }
   else if (*((REAL*) x) >  *((REAL*) y)) {
      return -1;
   }
   else {
      return  0;
   }
}

/* ==================================================================== */

int real_compar_up(const void *x, const void *y) {
   if (*((REAL*) x) >  *((REAL*) y)) {
      return  1;
   }
   else if (*((REAL*) x) <  *((REAL*) y)) {
      return -1;
   }
   else {
      return  0;
   }
}

/* ==================================================================== */

int rank_compar(const void *x, const void *y) {
  if (((SORTREC*) x)->val <  ((SORTREC*) y)->val) {
     return  1;
  }
  else if (((SORTREC*) x)->val >  ((SORTREC*) y)->val) {
     return -1;
  }
  else {
     return  0;
   }
}

/* ==================================================================== */


int u_rank_compar(const void *x, const void *y) {
  if (((SORTREC*) x)->val >  ((SORTREC*) y)->val) {
     return  1;
  }
  else if (((SORTREC*) x)->val <  ((SORTREC*) y)->val) {
     return -1;
  }
  else {
     return  0;
   }
}

/* ==================================================================== */

int wilcoxon_rank_compar(const void *x, const void *y) {
  if ( fabs(((SORTREC*) x)->val) >  fabs(((SORTREC*) y)->val) ) {
     return  1;
  }
  else if ( fabs(((SORTREC*) x)->val) <  fabs(((SORTREC*) y)->val) ) {
     return -1;
  }
  else {
     return  0;
   }
}

/* ==================================================================== */


REAL get_rank_correlation(REAL x[], REAL y[], int n) {
  SORTREC *xrec, *yrec;
  int i, k, found, j;
  REAL  d=0., n_n, rho, m, tx, ty;

  xrec = (SORTREC*)m_calloc(n, sizeof(SORTREC));
  yrec = (SORTREC*)m_calloc(n, sizeof(SORTREC));

  for (i=0; i<n; i++) {
    xrec[i].val = x[i];
    yrec[i].val = y[i];
    xrec[i].ind = i;
    yrec[i].ind = i;
  }

  qsort(xrec, n, sizeof(SORTREC), rank_compar);
  qsort(yrec, n, sizeof(SORTREC), rank_compar);

  i = 0;
  k = 0;
  m = 0.;
  tx = 0.;
  while (i<n) {
    yrec[i].rank = (REAL)i;
    if ( (i<(n-1)) && (equal(yrec[i].val, yrec[i+1].val)) ) {
      k++;
      m += (REAL)i;
    }
    else {
      if (k!=0) {
	m += (REAL)i;
	k ++;
	/* tx += pow((REAL)k, 3.) - (REAL)k; */
	tx += (REAL)k * (SQR((REAL)k)-1.);
	m = m/ (REAL)(k);
	for (j=i; j>(i-k); j--) {
	  yrec[j].rank = m;
	}
      }
      k=0;
      m=0.;
    }
    i++;
  }

  i = 0;
  k = 0;
  m = 0.;
  ty = 0.;
  while (i<n) {
    xrec[i].rank = (REAL)i;
    if ( (i<(n-1)) && (equal(xrec[i].val, xrec[i+1].val)) ) {
      k++;
      m += (REAL)i;
    }
    else {
      if (k!=0) {
	m += (REAL)i;
	k ++;
	ty += (REAL)k * (SQR((REAL)k)-1.);
	m = m/ (REAL)(k);
	for (j=i; j>(i-k); j--) {
	  xrec[j].rank = m;
	}
      }
      k=0;
      m=0.;
    }
    i++;
  }

  for (i=0; i<n; i++) {
    found = FALSE;
    for (k=0; k<n; k++) {
      if (xrec[i].ind == yrec[k].ind) {
	d += SQR((xrec[i].rank-yrec[k].rank));
	found = TRUE;
	break;
      }
    }
    if (!found) {
      out_err(MAT, ERR_FILE, ERR_LINE, _("One value could not be found!") );
    }
  }
  tx *= 0.5;
  ty *= 0.5;
  d *= 6.;
  /*  FALSCH: n_n = (REAL)(n) * SQR((REAL)n) -1.;  */
  n_n = (REAL)n * (SQR((REAL)n) - 1.);
  RDIV(rho, d, (n_n-(tx+ty)));
  rho = 1. - rho;
  /* printf("rho=%f\n", rho);  */
  /* rho = 1. - (d/(n_n-(tx+ty)));  */
  return rho;
}

/* ==================================================================== */


REAL get_sum(REAL x[], int n) {
  int i;
  REAL sum=0.;
  for (i=0; i<n; i++) {
    sum += x[i];
  }
  return sum;
}

REAL get_qsum(REAL x[], int n) {
  int i;
  REAL q_sum=0.;
  for (i=0; i<n; i++) {
    q_sum += SQR(x[i]);
  }
  return q_sum;
}


REAL get_xysum(REAL x[], REAL y[], int n) {
  int i;
  REAL xy_sum=0.;
  for (i=0; i<n; i++) {
    xy_sum += (x[i] * y[i]);
  }
  return xy_sum;
}

REAL get_sdv(REAL x[], int n) {
  REAL q_sum, sum, sdv;
  q_sum = get_qsum(x, n);
  sum = get_sum(x, n);
  sdv =  sqrt((q_sum - (sum*sum)/(REAL) n)/
		 (REAL) (n-1));
  return sdv;
}

REAL get_mean(REAL x[], int n) {
  REAL mean;
  mean = get_sum(x,n)/(REAL) n;
  return mean;
}

REAL get_min(REAL x[], int n) {
  REAL min=REAL_MAX;
  int i;
  for (i=0; i<n; i++) {
     if (x[i]<min) {min = x[i];}
  }
  return min;
}

REAL get_max(REAL x[], int n) {
  REAL max=REAL_MIN;

  int i;
  for (i=0; i<n; i++) {
     if (x[i]>max) {max = x[i];}
  }
  return max;
}


int get_maxint(int x[], int n) {
  int max=INT_MIN;
  int i;
  for (i=0; i<n; i++) {
     if (x[i]>max) {
       max = x[i];
     }
  }
  return max;
}


REAL get_norm_int(REAL x) {
  REAL nint, t, a[5], p, arg, erf, sqrttwo;
  BOOLEAN positiv=TRUE;

  if (x < 0.) {
    x *= -1.;
    positiv = FALSE;
  }
  a[0] = 0.254829592;
  a[1] = -0.284496736;
  a[2] = 1.421413741;
  a[3] = -1.453152027;
  a[4] = 1.061405429;
  p = 0.3275911;
  sqrttwo = 1.414213562373095;
  arg = x/sqrttwo;
  t = 1./(1.+p*arg);
  erf = 1. - (a[0]*t + a[1]*t*t + a[2]*t*t*t + a[3]*t*t*t*t +
	 a[4]*t*t*t*t*t)*exp(-arg*arg);
  nint = 0.5 * (1. + erf);
  if (!positiv) {
    nint = 1. - nint;
  }
  return nint;
}



REAL get_norm_ord(REAL x) {
  const REAL pi=3.141592654;
  REAL y;
  y = (1./sqrt(2.*pi)) * exp(-0.5*SQR(x));
  return y;
}


REAL get_median(REAL x[], int n) {
  REAL *xh, median;
  int  i;

  xh = (REAL*)m_calloc(n, sizeof(REAL));

  for (i=0; i<n; i++) {
    xh[i] = x[i];
  }

  qsort(xh, n, sizeof(REAL), real_compar_up);

  if (n % 2 == 1) {
    i = (n-1)/2;
    median = xh[i];
  }
  else {
    i = (n/2) -1;
    median = xh[i];
    i++;
    median = (median + xh[i])/2.;
  }

  return median;
}


int get_round(REAL x) {
  int z;
  REAL rem;
  z = (int) x;
  rem = x - (REAL)z;
  if (rem >= 0.5) {
    z ++;
  }
  return z;
}

REAL get_t_int(REAL t, int f)  {
  REAL s, s3, s4, s5, s6, s7, s8, t1, df;

  df = (REAL) f;
  s = 1.0;
  /* t2 = get_sgn(t); */
  t = SQR(t);
  if (t < 1.) goto label1;
  s3 = 1.;
  s4 = df;
  t1 = t;
  goto label2;
  label1: s3 = df;
  s4 = 1.;
  t1 = 1./t;
  label2: s5 = 2./9./s3;
  s6 = 2./9./s4;
  s7 = fabs((1.-s6)*pow(t1,(1./3.))-1.+s5) / sqrt(s6*pow(t1,(2./3.))+s5);
  if (s4 > 4.) goto label3;
  s7 *= 1. + .08 * pow(s7,4.) / pow(s4,3.);
  label3: s8 = 0.5/pow((1.+s7*(.196854+s7*(.115194+s7*(.000344+s7*.019527)))),4.);
  if (t < 1.) goto label4;
  s8 = 1. - s8;
  label4: ;
  s =  floor(1.e6*s8)/1.e6;
/*  s1 = .5 * (1.+s*t2);
    s2 = .5 * (1.-s*t2);
  alpha = 1.-s;
*/
  return s;
}


REAL get_chi_int(REAL chi, int f)  {
  REAL k, s4, s5, s6, s7, w, df, alpha;
  int i;

  if (chi > 100.) {
    alpha = 0.0;
    return alpha;
  }

  df = (REAL) f;
  k=1.;
  for (i=f; i>1; i-=2) {
    k *= (REAL)i;
  }
  s4 = pow(chi, floor((df+1.)/2.)) * exp(-chi/2.)/k;
  if (f%2 == 0) goto label1;
  s5 = sqrt(2./chi/3.1415927);
  goto label2;
  label1: s5 = 1.;
  label2: s6 = 1.;
  s7 = 1.;
  label4: df += 2.;
  s7 *= chi/df;
  if (s7 < .00001) goto label3;
  s6 += s7;
  goto label4;
  label3: w = s4*s5*s6;
  alpha = 1.-w;
  return alpha;
}



REAL get_f_int(REAL f, int f1, int f2) {
  REAL  df1, df2, s1,s3,s4,s5,s6,s7,s8,t1;

  df2 = (REAL) f2;
  df1 = (REAL) f1;
  /* s = 1.; */
  if (f < 1.) goto label1;
  s3 = df1;
  s4 = df2;
  t1 = f;
  goto label2;
  label1: s3 = df2;
  s4 = df1;
  t1 = 1./f;
  label2: s5 = 2./9./s3;
  s6 = 2./9./s4;
  s7 = fabs((1.-s6)*pow(t1,(1./3.))-1.+s5) / sqrt(s6*pow(t1,(2./3.))+s5);
  if (s4 >= 4.) goto label3;
  s7 *= 1.+.08*pow(s7,4.)/pow(s4,3.);
  label3: ;
  s8 = .5/pow((1.+s7*(.196854+s7*(.115194+s7*(.000344+s7*.019527)))),4.);
  if (f < 1.) goto label4;
  s8 = 1.-s8;
  label4: s1 = floor(1.e6*s8)/1.e6;
/*    s2 = 1.-s1;
    s = fabs(2.*s1-1.);
    alpha = 1.-s1;
*/
  return s1;
}



REAL rise(REAL x, int exp) {
  REAL y=1.;
  int i;

  for (i=0; i<exp; i++) {
    y *= x;
  }

  return y;
}



REAL get_z(REAL alpha) {
  REAL a[3], b[4], z, t;

  a[0] = 2.515517;
  a[1] = 0.802853;
  a[2] = 0.010328;
  b[1] = 1.432788;
  b[2] = 0.189269;
  b[3] = 0.001308;

  t = sqrt(-2. * log(1.-alpha));
  z = t - (a[0] + a[1]*t + a[2]*rise(t,2)) /
    (1. + b[1]*t + b[2]*rise(t,2) + b[3]*rise(t,3));

  return z;
}



REAL get_t(REAL alpha, int df) {
  REAL t, z, c9, c7, c5, c3, c1, n;

  z = get_z(alpha);
  n = (REAL)df;
  c9 = 79.;
  c7 = 720.*n;
  c5 = 4800.*rise(n,2) + 4560.*n + 1482.;
  c3 = 23040.*rise(n,3) + 15360.*rise(n,2);
  c1 = 92160.*rise(n,4) + 23040.*rise(n,3) + 2880.*rise(n,2) -3600.*n - 945.;

  t = ( c9*rise(z,9) + c7*rise(z,7) + c5*rise(z,5) + c3*rise(z,3) + c1*z ) /
    (92160.*rise(n,4));

  return t;
}



REAL get_sgn(REAL x) {
  if (x > 0.0) {
    return 1. ;
  }
  else if (x < 0.0) {
    return -1.;
  }
  else  {
    return 0.;
  }
}


REAL get_ln_0(REAL x) {
  if (x > 0.0) {
    return log(x);
  }
  else if (x < 0.0) {
    out_err(FAT, ERR_FILE, ERR_LINE,
	_("Log argument < 0!"));
    return 0.0;
  }
  else {
    return 0.0;
  }
}


int pks(REAL d, int n) {
  const REAL sig[30][4]= {
    { 0.0000, 0.0000, 0.0000, 0.0000 },
    { 0.0000, 0.0000, 0.0000, 0.0000 },
    { 0.3666, 0.3758, 0.3812, 0.3830 },
    { 0.3453, 0.3753, 0.4007, 0.4134 },
    { 0.3189, 0.3431, 0.3755, 0.3970 },
    { 0.2972, 0.3234, 0.3523, 0.3708 },
    { 0.2802, 0.3043, 0.3321, 0.3509 },
    { 0.2652, 0.2880, 0.3150, 0.3332 },
    { 0.2523, 0.2741, 0.2999, 0.3174 },
    { 0.2411, 0.2619, 0.2869, 0.3037 },
    { 0.2312, 0.2514, 0.2754, 0.2916 },
    { 0.2225, 0.2420, 0.2651, 0.2810 },
    { 0.2148, 0.2336, 0.2559, 0.2714 },
    { 0.2077, 0.2261, 0.2476, 0.2627 },
    { 0.2013, 0.2192, 0.2401, 0.2549 },
    { 0.1954, 0.2129, 0.2332, 0.2476 },
    { 0.1901, 0.2071, 0.2270, 0.2410 },
    { 0.1852, 0.2017, 0.2212, 0.2349 },
    { 0.1807, 0.1968, 0.2158, 0.2292 },
    { 0.1765, 0.1921, 0.2107, 0.2238 },
    { 0.1725, 0.1878, 0.2060, 0.2188 },
    { 0.1688, 0.1838, 0.2015, 0.2141 },
    { 0.1653, 0.1800, 0.1947, 0.2079 },
    { 0.1620, 0.1764, 0.1936, 0.2056 },
    { 0.1589, 0.1730, 0.1889, 0.2018 },
    { 0.1560, 0.1699, 0.1865, 0.1981 },
    { 0.1533, 0.1670, 0.1833, 0.1947 },
    { 0.1507, 0.1642, 0.1802, 0.1915 },
    { 0.1483, 0.1615, 0.1773, 0.1884 },
    { 0.1460, 0.1589, 0.1746, 0.1855 }
  };

  const REAL approx[4] = {0.8255, 0.8993, 0.9885, 1.0500};

  int alpha, i, k;
  REAL critical[4];

  for (i=0; i<4; i++) {
    if (n<=30) {
      critical[i] = sig[n-1][i];
    }
    else {
      critical[i] = approx[i]/sqrt((REAL)n);
    }
  }

  k = 4;
  for (i=0; i<4; i++) {
    if (d < critical[i]) {
      k = i;
      break;
    }
  }

  out_d(_("Critical values for d (two-sided):\n") );
  out_d("  10%%     5%%      2%%      1%%\n");
  out_d("%6.4f  %6.4f  %6.4f  %6.4f\n", critical[0], critical[1],
	critical[2], critical[3]);

  switch (k) {
  case 0: alpha = 10;
    break;
  case 1: alpha = 5;
    break;
  case 2: alpha = 2;
    break;
  case 3: alpha = 1;
    break;
  default: alpha = 0;
    break;
  }
/*  printf("alpha = %i\n", alpha);  */
  return alpha;
}

/* =================================================================== */

REAL get_multiple_reg(REAL y[], PREAL x[], int nrow, int ncol,
		      REAL b[], REAL *sdv, REAL *f_calc) {
   PREAL q[MCOL], e, a;
   REAL  br, cr, sr, tr, jor, fr, kr, help;
   int   k, i, j, jo, t, s, f1, f2;

   e = (REAL*)m_calloc((ncol+2), sizeof(REAL));
   a = (REAL*)m_calloc((ncol+2), sizeof(REAL));
   for (i=0; i<(ncol+1); i++) {
     q[i] = (REAL*)m_calloc((ncol+2), sizeof(REAL));
   }

   e[ncol+1] = 0.0;

   for (i=0; i<(ncol+1); i++) {
      for (j=0; j<(ncol+2); j++) {
         q[i][j] = 0.0;
      }
   }

   for (k=0; k<nrow; k++) {
     e[ncol+1] +=  SQR(y[k]);
     q[0][ncol+1] += y[k];
     e[0] = q[0][ncol+1];
     for (i=0; i<ncol; i++) {
       q[0][i+1] += x[i][k];
       q[i+1][0] = q[0][i+1];
       q[i+1][ncol+1] += x[i][k] * y[k];
       e[i+1] = q[i+1][ncol+1];
       for (j=i; j<ncol; j++) {
         q[j+1][i+1] = q[i+1][j+1] + x[i][k]*x[j][k];
         q[i+1][j+1] = q[j+1][i+1];
       }
     }
   }

   q[0][0] = (REAL) nrow;

   for (i=1; i<(ncol+1); i++) {
      a[i] = q[0][i];
   }

/* Begin of Gauss-Elimination */
   for (s=0; s<(ncol+1); s++) {
     t = s;

     label1: ;
     if (q[t][s] != 0.0) {
       goto label2;
     }
     t ++;
     if (t < ncol) {
       goto label1;
     }
     out_err(MAT, ERR_FILE, ERR_LINE,
     		_("Gauss-Elimination: No solution exists!") );
     return REAL_MIN;

     label2: ;
     for (jo=0; jo<(ncol+2); jo++) {
       br = q[s][jo];
       q[s][jo] = q[t][jo];
       q[t][jo] = br;
     }
     cr = 1./q[s][s];
     for (jo=0; jo<(ncol+2); jo++) {
       q[s][jo] *= cr ;
     }
     for (t=0; t<(ncol+1); t++) {
       if (t == s) {
         goto label3;
       }
       cr = -q[t][s];
       for (jo=0; jo<(ncol+2); jo++) {
         q[t][jo] += (cr * q[s][jo]);
       }
     label3: ;
     }
   }
/* End of Gauss-Elimination */

   sr = 0.0;
   for (i=1; i<(ncol+1); i++) {
     sr +=  q[i][ncol+1]*(e[i] - a[i]*e[0]/(REAL)nrow);
   }
   tr = e[ncol+1] - SQR(e[0])/(REAL)nrow;
   cr = tr - sr;
   f2 = nrow-ncol-1;
   jor = sr/(REAL)ncol;
   kr = cr/(REAL)f2;
   f1 = ncol;
   fr = jor/kr;
   br = sr/tr;
   RSQRT(kr, br);
   /*kr = sqrt(br);  */
   help = cr/(REAL)f2;
   RSQRT(sr, help);

   /* sr = sqrt(cr/(REAL)f2); */

   for (i=0; i<(ncol+1); i++) {
     b[i] = q[i][ncol+1];
   }
   *sdv = sr;
   *f_calc = fr;

   return br;
}

/* =================================================================== */

REAL get_cross_validate(REAL y[], PREAL x[], int nrow, int ncol,
			REAL ypred[]) {
  REAL sdv, f_calc;
  int i, j, k, n;
  PREAL xout[MCOL];
  REAL *yout, *b;
  REAL rquad, qquad, press=0.0, ssy=0.0;
  REAL y_mean;


  yout = (REAL*) m_calloc(nrow, sizeof(REAL));
  b = (REAL*) m_calloc(ncol+1, sizeof(REAL));
  for (i=0; i<ncol; i++) {
    xout[i] = m_calloc(nrow, sizeof(REAL));
  }

  y_mean = get_mean(y, nrow);
  for (i=0; i<nrow; i++) {
    for (j=0; j<nrow; j++) {
      if (i != j) {
        if (j<i) {
          n = j;
        }
        else {
          n = j-1;
        }
        yout[n] = y[j];
        for (k=0; k<ncol; k++) {
          xout[k][n] = x[k][j];
          /* printf(" xout[%i][%i]=%6.2f", k, n, xout[k][n]); */
        }
        /* printf(" yout[%i]=%6.2f\n", j, yout[n]); */
      }
    }
    if ((rquad=get_multiple_reg(yout,xout,nrow-1,ncol,b,&sdv,&f_calc))==REAL_MIN) {
      return REAL_MIN;
    }

    ypred[i] = b[0];
    for (j=1; j<(ncol+1); j++) {
      ypred[i] += b[j] * x[j-1][i];
    }

    /* printf("i=%i y_ori=%f  y_pred=%f\n", i, y[i], ypred[i]);  */
    /* printf("yorig=%f  ypred=%f\n", y[i],ypred[i]);  */
    press += SQR(y[i]-ypred[i]);
    ssy += SQR(y[i] - y_mean);
  }
  qquad = 1.0 - (press/ssy);

  return qquad;
}


/* =================================================================== */



void copy_tupel(TUPEL* dest, const TUPEL* src) {
  /* uses m_calloc so the copy is only temporary */
  int i;
  dest->a = (unsigned short int*) m_calloc(src->n, sizeof(unsigned short int));
  dest->n = src->n;
  for (i=0; i<dest->n; i++) {
    dest->a[i] = src->a[i];
  }
  return;
}


void print_tupel(TUPEL atupel) {
  int i;
  printf("%i->", atupel.n);
  for (i=0; i<atupel.n; i++) {
    printf("%3i", atupel.a[i]);
  }
}


TUPEL* create_tupel(TUPEL *atupel, int ndata) {
  BOOLEAN again;
  int i, a, j;
  /* static int seed = 1; */
  /* TUPEL atupel;  */

  /* srand(seed++); */
  (*atupel).n = (short unsigned int) ndata;
  for (i=0; i<ndata; i++) {
    do {
      if (ndata<1000) {
	a = (rand()/13) % ndata;
      }
      else {
	a = rand() % ndata;
      }
      again = FALSE;
      for (j=0; j<i; j++) {
	if ((*atupel).a[j]==a) {
	  again = TRUE;
	  break;
	}
      }
    } while (again);
    (*atupel).a[i] = a;
  }
  return atupel;
}



BOOLEAN equal_tupel(TUPEL t1, TUPEL t2) {
  int i;

  if (t1.n != t2.n) {
    return FALSE;
  }
  for (i=0; i<t1.n; i++) {
    if (t1.a[i] != t2.a[i]) {
      return FALSE;
    }
  }
  return TRUE;
}


/* =================================================================== */

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
**  small Changes by Bernhard Reiter  http://www.usf.Uni-Osnabrueck.DE/~breiter
** $Id: procs.c,v 1.35 2006/11/15 17:17:34 jakson Exp $
**
**
** further changes: Andreas Beyer, 1999, abeyer@usf.uni-osnabrueck.de
***************************************************************/

/* procs.c for STATIST */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

#include "statist.h"
#include "funcs.h"
#include "plot.h"
#include "procs.h"
#include "data.h"
#include "gettext.h"


/* =================================================================== */


void part_corr(PREAL xx[], int nrow, int ncol) {
  PREAL r[6], x, y;
  REAL  x_sdv, y_sdv, s_xy;
  int   i, k, df;
  REAL r12s3, r12s4, r53s4, r25s4, r15s4, r12s5, r13s5,
        r23s5, r12s34, r15s34, r25s34, r14s5, r24s5, r12s35, r12s45,
        r12s345, r13s4, r23s4, t, alpha, a, r14s3, r13s2;
/*
   calculate partial correlation;
   at most 3 constant correlations
*/

  for (i=0; i<ncol; i++) {
    r[i+1] = (REAL*)m_calloc((ncol+1), sizeof(REAL));
    for (k=0; k<ncol; k++) {
      x = xx[i];
      y = xx[k];

      s_xy = (get_xysum(x,y,nrow) - ((get_sum(x,nrow)*get_sum(y,nrow))/(REAL) nrow))/
             (REAL) (nrow-1);
      x_sdv = get_sdv(x,nrow);
      y_sdv = get_sdv(y,nrow);
      DIV(r[i+1][k+1], s_xy, (x_sdv*y_sdv));
    }
  }

  DIV(r12s3, (r[1][2]-r[1][3]*r[2][3]), sqrt((1.-SQR(r[1][3]))*(1.-SQR(r[2][3]))));
  DIV(r13s2, (r[1][3]-r[1][2]*r[3][2]),sqrt((1.-SQR(r[1][2]))*(1.-SQR(r[3][2]))));

  out_start();
  colorize(ClHeader);
  out_r(_("\nResult partial correlation:\n") );
  colorize(ClDefault);
  out_r(" r12   = %f\n\n", r[1][2]);

  out_r(_("One constant correlation:") );
  out_r("\n r12_3 = %f  diff= %f\n", r12s3, (r[1][2]-r12s3));
  out_r(" r13_2 = %f  diff= %f\n", r13s2, (r[1][3]-r13s2));

  if (ncol >= 4) {
    DIV(r12s4, (r[1][2]-r[1][4]*r[2][4]), sqrt((1.-SQR(r[1][4]))*(1.-SQR(r[2][4]))));
    DIV(r13s4, (r[1][3]-r[1][4]*r[3][4]), sqrt((1.-SQR(r[1][4]))*(1.-SQR(r[3][4]))));
    DIV(r23s4, (r[2][3]-r[2][4]*r[3][4]), sqrt((1.-SQR(r[2][4]))*(1.-SQR(r[3][4]))));

    DIV(r14s3, (r[1][4]-r[1][3]*r[4][3]), sqrt((1.-SQR(r[1][3]))*(1.-SQR(r[4][3]))));

    DIV(r12s34, (r12s4-r13s4*r23s4), sqrt((1.-SQR(r13s4))*(1.-SQR(r23s4))));

    out_r(" r12_4 = %f  diff= %f\n", r12s4, r[1][2]-r12s4);
    out_r(" r13_4 = %f  diff= %f\n", r13s4, r[1][3]-r13s4);

    out_r(" r14_3 = %f  diff= %f\n\n", r14s3, r[1][4]-r14s3);
    out_r(_("Two constant correlations:"));
    out_r("\n r12_34 = %f  diff= %f\n", r12s34, r[1][2]-r12s34);

    if (ncol == 5) {
      DIV(r53s4, (r[5][3]-r[5][4]*r[3][4]), sqrt((1.-SQR(r[5][4]))*(1.-SQR(r[3][4]))));
      DIV(r25s4, (r[2][5]-r[2][4]*r[4][5]), sqrt((1.-SQR(r[2][4]))*(1.-SQR(r[4][5]))));
      DIV(r15s4, (r[1][5]-r[1][4]*r[4][5]), sqrt((1.-SQR(r[1][4]))*(1.-SQR(r[4][5]))));
      DIV(r12s5, (r[1][2]-r[1][5]*r[2][5]), sqrt((1.-SQR(r[1][5]))*(1.-SQR(r[2][5]))));
      DIV(r13s5, (r[1][3]-r[1][5]*r[3][5]), sqrt((1.-SQR(r[1][5]))*(1.-SQR(r[3][5]))));
      DIV(r23s5, (r[2][3]-r[2][5]*r[3][5]), sqrt((1.-SQR(r[2][5]))*(1.-SQR(r[3][5]))));

      DIV(r15s34, (r15s4-r13s4*r53s4), sqrt((1.-SQR(r13s4))*(1.-SQR(r53s4))));
      DIV(r25s34, (r25s4-r23s4*r53s4), sqrt((1.-SQR(r23s4))*(1.-SQR(r53s4))));
      DIV(r14s5, (r[1][4]-r[1][5]*r[4][5]), sqrt((1.-SQR(r[1][5]))*(1.-SQR(r[4][5]))));
      DIV(r24s5, (r[2][4]-r[2][5]*r[4][5]), sqrt((1.-SQR(r[2][5]))*(1.-SQR(r[4][5]))));
      DIV(r12s35, (r12s5-r13s5*r23s5), sqrt((1.-SQR(r13s5))*(1.-SQR(r23s5))));
      DIV(r12s45, (r12s5-r14s5*r24s5), sqrt((1.-SQR(r14s5))*(1.-SQR(r24s5))));

      DIV(r12s345, (r12s34-r15s34*r25s34), sqrt((1.-SQR(r15s34))*(1.-SQR(r25s34))));

      out_r(" r12_5 = %f  diff= %f\n", r12s5, r[1][2]-r12s5);
      out_r(" r12_35 = %f  diff= %f\n", r12s35, r[1][2]-r12s35);
      out_r(" r12_45 = %f  diff= %f\n\n", r12s45, r[1][2]-r12s45);
      out_r(_("Three constant correlations:") );
      out_r("\n r12_345 = %f  diff= %f\n", r12s345, r[1][2]-r12s345);
    }
  }

   df = nrow - 3;
   t = r12s3/(sqrt(1.-SQR(r12s3))) * sqrt((REAL)df);
   /* f = SQR(r12s3)/(1.-SQR(r12s3)) * (REAL)(nrow-2); */
   alpha = get_t_int(fabs(t), df);

   out_r( _("\nt-Test:\n"));
   out_r(_("Hypothesis H0: r12_3 = 0  against hypothesis H1: r12_3 # 0  "
	 "(->two sided)\n") );
   a = 1.0-alpha;
   out_r(_("Probability of H0 = %6.4f\n\n"), a);
/*
   out_r("F-Test:\n");
   alpha = get_f_int(f, 1, df);
   a = alpha;
   out_r("Hypothese H0: r=0  gegen  Hypothese H1: r#0  (->zweiseitig)\n");
   out_r("Irrtumswahrscheinlichkeit alpha fuer H0 = %6.4f\n\n", a);
*/
  out_end();

}

/* =================================================================== */


void multiple_reg(REAL y[], PREAL x[], int nrow, int ncol) {
   REAL *b, *resi, sdv, rquad, alpha, f_calc, reg, theo;
   int i, j, rows;
   FILE *rfile;
   BOOLEAN hasmiss = FALSE;

   b = (REAL*)m_calloc(ncol+1, sizeof(REAL));
   if ((rquad=get_multiple_reg(y, x, nrow, ncol, b, &sdv, &f_calc))==REAL_MIN) {
     return;
   }

   out_start();
   colorize(ClHeader);
   out_r(_("\nResults multiple linear regression:\n") );
   colorize(ClDefault);
   out_r(_("Regressed equation: y = B[0] + B[1]*x[1] + B[2]*x[2] +...+") );
   out_r(" B[n]*x[n] \n\n");

   for (i=0; i<(ncol+1); i++) {
     out_r("B[%i] = %f\n", i, b[i]);
   }
   out_r("\n");


   if (!noplot) {
     if (ncol==2) {
       plot_tripl(x[0], x[1], y, nrow, b[0], b[1], b[2],
	   get_label(x[0]), get_label(x[1]), get_label(y));
     }

     if (ncol==1) {
       plot_pair(x[0], y, nrow, b[0], b[1], get_label(x[0]), get_label(y));
     }
   }

   /* If there is any missing value, the residues will be in the wrong rows.
    * To avoid this, we need to realloc cols before calculating "resi" */
   rows= nrow;
   for(j = 0; j < ncol; j++)
     if(nn[acol[j]] != vn[acol[j]]){
       alloc_cols((ncol + 1), 3);
       rows = nn[acol[0]];
       break;
     }

   resi = (REAL*)m_calloc(rows, sizeof(REAL));

   for(j = 0; j < rows; j++){
     theo = b[0];
     for(i = 1; i <= ncol; i++){
       if(xx[acol[i]][j] == SYSMIS){
	 hasmiss = TRUE;
	 break;
       }
       theo += b[i] * xx[acol[i]][j];
     }
     if(hasmiss){
       resi[j] = SYSMIS;
       hasmiss = FALSE;
     } else
       resi[j] = xx[acol[0]][j] - theo;
   }

   if ((j=col_exist("resi", FALSE)) == -1) {
     if (!(make_new_col("resi", rows))) {
       out_end();
       return;
     }
   } else{
#ifndef STATIST_X
     out_r(_("Column 'resi' updated!\n") );
#endif
   }
   if ((j=col_exist("resi", FALSE)) != -1) {
     rfile = tmpptr[j];
   } else{
     return;
   }
   rewind(rfile);
   FWRITE(resi, sizeof(REAL), rows, rfile);
   rewind(rfile);

   SQRT(reg, rquad);
   out_r(_("Coefficient of determination r^2 = %f\n"), rquad);
   out_r(_("Correlation coefficient r = %f\n"), reg);
   out_r(_("Standard deviation s = %f\n"), sdv);
   out_r(_("Number of data points n = %i\n"), nrow);

   if (!equal(1.0, reg)) {
     out_r(_("\nF-value = %f\n"), f_calc);
     out_r(_("Degree of freedom f1 = %i\n"), ncol);
     out_r(_("Degree of freedom f2 = %i\n"), nrow-ncol-1);

     out_r(_("\nF-Test:\n"));
     out_r(_("Hypothesis H0: r=0  against hypothesis H1: r#0 \n") );
     alpha = get_f_int(f_calc, ncol, nrow-ncol-1);
     if (reg < 0.0) {
       alpha = 1.-alpha;
     }
     out_r(_("Probability of H0: %6.4f\n\n"), 1.0-alpha);
   }
   else {
     out_r(_("F-Test not possible since r = 1\n\n") );
   }
   out_end();
}

/* =================================================================== */



void cross_validate(REAL y[], PREAL x[], int nrow, int nc) {
  char *pred_label;
  REAL rquad_ori_pred, qquad, rquad, sdv, f_calc;
  REAL *b, *ypred, *ypredm;
  char analias[80];
  int i, j, k;
  BOOLEAN hasmiss = FALSE;


  b = (REAL*) m_calloc(nc+1, sizeof(REAL));
  ypred = (REAL*) m_calloc(nrow, sizeof(REAL));
  if ((qquad=get_cross_validate(y, x, nrow, nc, ypred))==REAL_MIN) {
    return;
  }
  if ((rquad=get_multiple_reg(y,x,nrow,nc,b,&sdv,&f_calc))==REAL_MIN) {
    return;
  }
  if ((rquad_ori_pred=get_multiple_reg(ypred,&y,nrow,1,b,&sdv,&f_calc))==REAL_MIN) {
    return;
  }

  if (!noplot) {
    pred_label = m_calloc(1,strlen(get_label(y))+strlen(" vorhergesagt")+1);
    strcpy(pred_label, get_label(y));
    strcat(pred_label, _(" predicted") );
    plot_pair(y, ypred, nrow, b[0], b[1], get_label(y), pred_label);
  }

  out_start();
  out_r(_("Coefficient of determination r^2 = %f\n"), rquad);
  out_r(_("r^2 of regression y versus y_pred = %f\n"), rquad_ori_pred);
  out_r(_("Variance of prediction q^2 = %f\n"), qquad);

  /* If we have missing values, the ypred values might be in the wrong rows.
   * In this case, we need to insert the missing values in "ypred". */
  ypredm = ypred;
  for(i = 0; i <= nc; i++)
    if(nn[acol[i]] != vn[acol[i]]){
      hasmiss = TRUE;
      nrow = nn[acol[0]];
      nc++;
      alloc_cols(nc, 3);
      break;
    }
  if(hasmiss){
    hasmiss = FALSE;
    k = 0;
    ypredm = (REAL*) m_calloc(nrow, sizeof(REAL));
    for(i = 0; i < nrow; i++){
      for(j = 0; j < nc; j++)
	if(xx[acol[j]][i] == SYSMIS){
	  hasmiss = TRUE;
	  break;
	}
      if(hasmiss){
	ypredm[i] = SYSMIS;
	hasmiss = FALSE;
      }
      else{
	ypredm[i] = ypred[k];
	k++;
      }
    }
  }
  
  /* create new data record: */
  strcpy(analias, "pred_");
  strncat(analias, get_name(xx[acol[0]]), 79-strlen(analias));
  if (!(make_new_col(analias, nrow))) {
    out_end();
    return;
  }
  
  FWRITE(ypredm, sizeof(REAL), nrow, tmpptr[ncol - 1]);
  out_end();
}

/* ================================================================== */


void random_tupel(REAL y[], PREAL x[], int nrow, int nc, int ntup) {
  int  i, k, j;
  BOOLEAN found;
  TUPEL *tupel_list, atupel;
  REAL rquad, sdv, f_calc;
  REAL *yperm, *b, *ypred, qquad;
  FILE *fp1, *fp2;

  yperm = (REAL*) m_calloc(nrow, sizeof(REAL));
  ypred = (REAL*) m_calloc(nrow, sizeof(REAL));
  atupel.a = (unsigned short int*) m_calloc(nrow, sizeof(unsigned short int));
  atupel.n = nrow;
  tupel_list = (TUPEL*) m_calloc(ntup, sizeof(TUPEL));
  b = (REAL*)m_calloc(nc+1, sizeof(REAL));

  if ((rquad=get_multiple_reg(y, x, nrow, nc, b, &sdv, &f_calc))==REAL_MIN) {
    return;
  };
  if ((qquad=get_cross_validate(y, x, nrow, nc, ypred))==REAL_MIN) {
    return;
  }
  out_start();
  out_r(_("\nOriginal y-Vector: r^%6.4f  q^2=%6.4f\n\n"), rquad, qquad);

  if(!(make_new_col("rquad", ntup))){
    out_end();
    return;
  }
  fp1 = tmpptr[ncol - 1];
  if(!(make_new_col("qquad", ntup))){
    out_end();
    return;
  }
  fp2 = tmpptr[ncol - 1];   

  out_d(_("Starting with randomization of y-vector. Please be patient ...\n"));
  i=0;
  do {
    create_tupel(&atupel, nrow);
    /* print_tupel(atupel); */
    found = FALSE;
    for (k=0; k<i; k++) {
      if (equal_tupel(atupel, tupel_list[k])) {
	found = TRUE;
	break;
      }
    }

    if (!found) {
      copy_tupel(&tupel_list[i], &atupel);
      i++;
    }
    if (i%100==0) {
      out_d(".");
      fflush(stdout);
    }
  } while (i<ntup);
  out_d("\n");

  out_d(_("Calculating q^2 and r^2 of randomized vectors ...") );

#define _FINISH_RANDOM_TUPLE 	\
  delete_last_columns(2);	\
  out_end();
  	
	
  for (i=0; i<ntup; i++) {
    if (i%100==0) {
      out_d(".");
      fflush(stdout);
    }
    for (k=0; k<nrow; k++) {
      j = tupel_list[i].a[k];
      yperm[k] = y[j];
    }
    if ((rquad=get_multiple_reg(yperm, x, nrow, nc, b, &sdv, &f_calc))==REAL_MIN) {
	_FINISH_RANDOM_TUPLE
	return;
    }
    if ((qquad=get_cross_validate(yperm, x, nrow, nc, ypred))==REAL_MIN) {
	_FINISH_RANDOM_TUPLE
	return;
    }
    FWRITE(&rquad, sizeof(REAL), 1, fp1);
    FWRITE(&qquad, sizeof(REAL), 1, fp2);
  }
  out_d(_("\ndone!\n\n") );

#undef _FINISH_RANDOM_TUPLE
}

/* ================================================================= */

void poly_reg(REAL x[], REAL y[], int n, int m) {

   REAL a[2*MPOLY+1], e[MPOLY+2], q[MPOLY+1][MPOLY+2];

   REAL  br, cr, sr, tr, jor, fr, kr, pr, alpha;
   int   i, jo, t, s, f1, f2;

   for (i=1; i<(2*m+1); i++) {
     a[i] = 0.0;
   }

   for (i=0; i<(m+2); i++) {
     e[i] = 0.0;
   }

   a[0] = (REAL) n;

   for (i=0; i<n; i++) {
     pr = 1.;
     for (jo=1; jo<(2*m+1); jo++) {
       pr *= x[i];
       a[jo] += pr;
     }
     pr = y[i];
     for (jo=0; jo<(m+1); jo++) {
       e[jo] += pr;
       pr *= x[i];
       q[jo][m+1] = e[jo];
     }
     e[m+1] += SQR(y[i]);
   }

   for (i=0; i<(m+1); i++) {
     for (jo=0; jo<(m+1); jo++) {
       q[i][jo] = a[i+jo];
     }
   }

/* Begin of Gauss-Elimination */
   for (s=0; s<(m+1); s++) {
     t = s;
     label1: ;
     if (q[t][s] != 0.0) {
       goto label2;
     }
     t ++;
     if (t < m) {
       goto label1;
     }
     out_err(MAT, ERR_FILE, ERR_LINE,
	 _("Gauss-Elimination: No possible solution!") );
     return;
     label2: ;
     for (jo=0; jo<(m+2); jo++) {
       br = q[s][jo];
       q[s][jo] = q[t][jo];
       q[t][jo] = br;
     }
     cr = 1./q[s][s];
     for (jo=0; jo<(m+2); jo++) {
       q[s][jo] *= cr ;
     }
     for (t=0; t<(m+1); t++) {
       if (t == s) {
         goto label3;
       }
       cr = -q[t][s];
       for (jo=0; jo<(m+2); jo++) {
         q[t][jo] += (cr * q[s][jo]);
       }
     label3: ;
     }
   }
/* End of Gauss-Elimination */

   sr = 0.0;
   for (i=1; i<(m+1); i++) {
     sr += q[i][m+1] * (e[i]-a[i]*e[0]/(REAL)n);
   }

   tr = e[m+1] - e[0]*e[0]/(REAL)n;
   cr = tr-sr;
   f2 = n-m-1;
   DIV(jor, sr, (REAL)m);
   DIV(kr, cr, (REAL)f2);
   f1 = m;
   DIV(fr, jor, kr);
   DIV(br, sr, tr);
   SQRT(kr, br);
   DIV(sr, cr, (REAL)f2);
   SQRT(sr, sr);

   out_start();
   for (i=0; i<(m+1); i++) {
     a[i] = q[i][m+1];
     out_r(_("Coefficient a%i = %15e\n"), i, a[i]);

   }

   colorize(ClHeader);
   out_r(_("\nResult polynomial regression:\n") );
   colorize(ClDefault);
   out_r(_("Regressed equation: y = a0 + a1*x + a2*x^2 +...+ an*x^n\n\n") );

   if (!noplot)
     plot_poly(x, y, n, a, m, get_label(x), get_label(y));

   out_r(_("Coefficient of determination r^2 = %f\n"), br);
   out_r(_("Correlation coefficient r = %f\n"), kr);
   out_r(_("Standard deviation s = %f\n"), sr);

   if (!equal(1.0, kr)) {
     out_r(_("F-value = %f\n"), fr);
     out_r(_("Degree of freedom f1 = %i\n"), f1);
     out_r(_("Degree of freedom f2 = %i\n"), f2);

     out_r(_("\nF-Test:\n"));
     out_r(_("Hypothesis H0:  r=0  against hypothesis H1: r>0 or r<0\n") );
     alpha = get_f_int(fr, f1, f2);
     if (kr < 0.0) {
       alpha = 1.-alpha;
     }
     out_r(_("Probability of H0: %6.4f\n\n"), 1.0-alpha);
   }
   else {
     out_r(_("F-Test not possible since r = 1\n\n") );
   }
   out_end();
}

/* =================================================================== */


void correl_matrix(PREAL xx[], int nrow, int ncol) {
   int i, k;
   REAL *x, *y, x_sdv, y_sdv, s_xy;
   /*   PREAL *r;    */
   PREAL r[MCOL];
   char label[10];

   /* r = (PREAL*)m_calloc(ncol, sizeof(PREAL)); */
   for (i=0; i<ncol; i++) {
     r[i] = (REAL*) m_calloc(ncol, sizeof(REAL));
   }

   for (i=0; i<ncol; i++) {
     r[i][i] = 1.0;
     for (k=0; k<i; k++) {
       x = xx[i];
       y = xx[k];
       s_xy = (get_xysum(x,y,nrow) - ((get_sum(x,nrow)*get_sum(y,nrow))/
              (REAL) nrow))/ (REAL) (nrow-1);
       x_sdv = get_sdv(x,nrow);
       y_sdv = get_sdv(y,nrow);
       DIV(r[i][k], s_xy, (x_sdv*y_sdv));
       r[k][i] = r[i][k];
     }
   }

   out_start();
   for(i = 0; i < 10; i++)
     label[i] = 0;
   colorize(ClHeader);
   out_r(_("Matrix of linear correlation of variables:\n\n") );
   out_r(_("Variable     ") );
   for (i=0; i<ncol; i++) {
     /*  out_r("%8i", (i+1)); */
     strncpy(label, get_name(xx[i]), 9);
     if (strlen(label)>6) {
       label[6] = '.';
       label[7] = '\0';
     }
     out_r("%-8s", label);
   }
   colorize(ClDefault);
   out_r("\n");
   for (i=0; i<ncol; i++) {
     out_r("--------");
   }
   out_r("------------\n");

   for (i=0; i<ncol; i++) {
     /*  out_r("      %2i |  ", (i+1)); */
     strncpy(label, get_name(xx[i]), 9);
     if (strlen(label)>6) {
       label[6] = '.';
       label[7] = '\0';
     }
     colorize(ClHeader);
     out_r(" %-7s", label);
     colorize(ClDefault);
     out_r(" | ");
     for (k=0; k<ncol; k++) {
       out_r("%8.4f", r[i][k]);
     }
     out_r("\n");
   }
   out_r("\n");
   out_end();
}


/* =================================================================== */

void rank_matrix(PREAL xx[], int nrow, int ncol) {
   int i, k;
   REAL *x, *y;
   /*   PREAL *r; */
   PREAL r[MCOL];

   char label[64];

   /* r = (PREAL*)m_calloc(ncol, sizeof(PREAL)); */
   for (i=0; i<ncol; i++) {
     r[i] = (REAL*) m_calloc(ncol, sizeof(REAL));
   }

   for (i=0; i<ncol; i++) {
     r[i][i] = 1.0;
     for (k=0; k<i; k++) {
       x = xx[i];
       y = xx[k];
       if ((r[i][k]=get_rank_correlation(x, y, nrow))==REAL_MIN) {
	 return;
       }
       r[k][i] = r[i][k];
     }
   }

   out_start();
   colorize(ClHeader);
   out_r(_("Matrix of SPEARMAN'S Rank-Correlation-Coefficients\n") );
   out_r(_("of the variables:\n\n") );
   out_r(_("Variable     ") );
   for (i=0; i<ncol; i++) {
     /*  out_r("%8i", (i+1)); */
     strncpy(label, get_name(xx[i]), 9);
     if (strlen(label)>6) {
       label[6] = '.';
       label[7] = '\0';
     }
     out_r("%-8s", label);
   }
   colorize(ClDefault);
   out_r("\n");
   for (i=0; i<ncol; i++) {
     out_r("--------");
   }
   out_r("------------\n");

   for (i=0; i<ncol; i++) {
     /*  out_r("      %2i |  ", (i+1)); */
     strncpy(label, get_name(xx[i]), 9);
     if (strlen(label)>6) {
       label[6] = '.';
       label[7] = '\0';
     }
     colorize(ClHeader);
     out_r(" %-7s", label);
     colorize(ClDefault);
     out_r(" | ");
     for (k=0; k<ncol; k++) {
       out_r("%8.4f", r[i][k]);
     }
     out_r("\n");
   }
   out_r("\n");
  out_end();
}


/*    out_r("Variable "); */
/*    for (i=0; i<ncol; i++) { */
/*      out_r("%8i", (i+1)); */
/*    } */
/*    out_r("\n"); */
/*    for (i=0; i<ncol; i++) { */
/*      out_r("--------"); */
/*    } */
/*    out_r("------------\n"); */

/*    for (i=0; i<ncol; i++) { */
/*      out_r("      %2i |  ", (i+1)); */
/*      for (k=0; k<ncol; k++) { */
/*        out_r("%8.4f", r[i][k]); */
/*      } */
/*      out_r("\n"); */
/*    } */
/*    out_r("\n"); */
/*  } */


/* =================================================================== */

void probit(REAL dose[], REAL num[], REAL effect[], int n) {
  REAL propcalc, chisq, ycalc;
  REAL *aprobit, *logx, neq, yh, w, oc;
  int i, np, offset, ind, old, nerr = 0;
  BOOLEAN positiv, equal_num, equal_dose, *infinity, show_errors = TRUE;
  REAL x_mean, x_sdv, y_mean, y_sdv;
  REAL s_xy, r_xy, a1, a0, ed50, ed90, alpha, tau, arg;
  REAL  nw_sum, nwx_sum, nwxsqr_sum, s_m, conf_l, conf_u, mq_m;

  REAL dose0, dose1;

  aprobit = (REAL*)m_calloc(n, sizeof(REAL));
  logx = (REAL*)m_calloc(n, sizeof(REAL));
  infinity = (BOOLEAN*)m_calloc(n, sizeof(BOOLEAN));

  np = 0;

  out_start();
  for (i=0; i<n; i++) {
    yh = effect[i]/num[i];
    if (yh < 0.50) {
      positiv = TRUE;
      yh = 1. - yh;
    }
    else {
      positiv = FALSE;
    }

    if (yh >= 1.0) {
      infinity[i] = TRUE;
      nerr++;
      if(nerr == 30){
	out_d("\n\n");
	out_i(_("Continue showing these math warnings? (%s) "), _("y/N"));
	GETNLINE;
	if(!(empty) && (line[0] != _("y")[0] || line[0] != _("Y")[0])){
	  show_errors = FALSE;
	}
      }
      if(show_errors){
	out_err(MWA, ERR_FILE, ERR_LINE,
	    _("Can not calculate probit: dose %g count %g effect %g"),
	    dose[i], num[i], effect[i]);
      }
    }
    else {
      infinity[i] = FALSE;

      neq = 0.;
      if (!equal(yh, 0.5)) {
          neq = get_z(yh);
      }

      if (positiv) {
        neq = -1. * neq;
      }
      aprobit[np] = neq + 5.;

      if (dose[i] <= 0.0) {
        out_err(MAT, ERR_FILE, ERR_LINE, _("Dose %i <= 0!\n"), i);
        out_end();
        return;
      }

      logx[np] = log10(dose[i]);
      out_r(_("dose=%6g  num=%g effect=%2f prop=%4.2f probit=%5.2f\n"),
             dose[i], num[i], effect[i], (effect[i]/num[i]), aprobit[np]);
      np++;
    }
    if(show_errors){
      if((np + nerr + 1) % 12 == 0)
	mywait();
    } else
      if((np + 1) % 12 == 0)
	mywait();
  }


  if (np < 2) {
    out_err(MAT, ERR_FILE, ERR_LINE,
	_("Less than 2 valid dose-effect pairs!") );
    out_end();
    return;
  }


/* Test for equal probabilities of effects and equal doses, respectively */
  equal_num = TRUE;
  equal_dose = TRUE;
  i=0;
  while(infinity[i]) {
    i++;
  }
  old=i;
  for (ind=0; ind<np; ind++) {
    while(infinity[i]) {
      i++;
    }
    if ( (ind>0) && ((effect[i]/num[i])!=(effect[old]/num[old])) ) {
      equal_num = FALSE;
    }
    if ( (ind>0) && (dose[i] != dose[old]) ) {
      equal_dose = FALSE;
    }
    old=i;
    i++;
  }

  if (equal_num) {
    out_err(MAT, ERR_FILE, ERR_LINE, _("All effect probabilities are equal!"));
    out_end();
    return;
  }

  if (equal_dose) {
    out_err(MAT, ERR_FILE, ERR_LINE, _("All doses are equal!") );
    out_end();
    return;
  }


  mywait();

  x_mean = get_mean(logx,np);
  y_mean = get_mean(aprobit,np);
  x_sdv = get_sdv(logx,np);
  y_sdv = get_sdv(aprobit,np);
  s_xy = (get_xysum(logx,aprobit,np) - ((get_sum(logx,np)*get_sum(aprobit,np))/
          (REAL) np))/(REAL) (np-1);
  r_xy = s_xy/(x_sdv*y_sdv);
/*  rr = r_xy * r_xy; */
  a1 = r_xy * (y_sdv/x_sdv);
  a0 = y_mean - a1*x_mean;
  DIV(arg,  SQR(r_xy), (1.-SQR(r_xy)));
  SQRT(tau, ((np-2.)*arg));

  chisq = 0.;
  nw_sum = 0.0;
  nwx_sum = 0.0;
  nwxsqr_sum = 0.0;
  offset = 0;

  for (i=0; i<np; i++) {
    if (infinity[i]) {
      offset++;
    }
    ind = i+offset;
    ycalc = a0 + a1*logx[i];
    neq = ycalc-5;
    if (neq > 0.) {
      positiv = TRUE;
      neq = neq * (-1.);
    }
    else {
      positiv = FALSE;
    }

    propcalc = get_norm_int(neq);
    if (positiv) {
      propcalc = 1. - propcalc;
    }
    w = SQR(get_norm_ord(ycalc-5.))/(propcalc*(1.-propcalc));
/*    printf("%g  %g  %g\n", logx[i], ycalc-5, w);    */
    nw_sum += num[ind]*w;
    nwx_sum += num[ind]*w*logx[i];
    nwxsqr_sum += num[ind]*w*SQR(logx[i]);
    oc = SQR(effect[ind] - (num[ind]*propcalc))/
      (num[ind]*propcalc*(1.-propcalc));
    chisq += oc;
    /*
    out_r("i=%i ind=%i num=%f effect=%f y=%f P=%f w=%f oc=%f logx=%g\n",
          i, ind, num[ind], effect[ind], ycalc, propcalc, w, oc, logx[i]);
    */
  }


  /*
  out_r("np = %i\n", np);
  out_r("nwxsqr_sum = %g\n", nwxsqr_sum);
  out_r("nwx_sum    = %g\n", nwx_sum);
  out_r("nw_sum     = %g\n", nw_sum);
  */

  nwxsqr_sum = nwxsqr_sum - ( SQR(nwx_sum)/nw_sum );
  x_mean = nwx_sum/nw_sum;
  a0 = y_mean - a1*x_mean;
  ed50 = (5.0-a0)/a1;
  mq_m = 1./SQR(a1) * (1./nw_sum + SQR(ed50-x_mean)/nwxsqr_sum);
  s_m = sqrt(mq_m);
  conf_l = ed50-(1.96*s_m);
  conf_u = ed50+(1.96*s_m);
  ed90 = (6.28-a0)/a1;

  /*
  out_r("nwxsqr = %g\n", nwxsqr_sum);
  out_r("x_mean = %g\n", x_mean);
  out_r("mq_m   = %g\n", mq_m);
  out_r("s_m    = %g\n", s_m);
  */

  colorize(ClHeader);
  out_r(_("Result probit analysis:\n") );
  colorize(ClDefault);
  if (a1<0.0) {
    out_err(MWA, ERR_FILE, ERR_LINE,
	_("Inverse probit curve due to negative slope a1!!!") );
  }
  out_r("ED50 = %g\n", pow(10., ed50));
  out_r(_("Confidence range (95%%) for ED50: [%g, %g]\n"),
         pow(10., conf_l), pow(10., conf_u));
/*
  out_r("nw_sum=%f nwx_sum=%f x_mean=%f s_m=%f\n",
         nw_sum, nwx_sum, x_mean, s_m);
*/
  out_r("ED90 = %g\n", pow(10., ed90));
  out_r("a0   = %g\n", a0);
  out_r("a1   = %g\n", a1);
  out_r( _("Chi^2= %g\n"), chisq);
  out_r(_("Degrees of freedom = %i\n"),(np-2));
  out_r(_("Correlation coefficient r = %f\n"), r_xy);
  out_r(_("Check value Tau = %f\n"), tau);
  if (tau>0.0) {
    out_r(_("\nt-Test with check value Tau:\n") );
    out_r(_("Hypothesis H0: Correlation according to probit-model exists\n") );
    alpha = get_t_int(tau, (np-2));
    out_r(_("Probability of H0: %f\n"), alpha);
  }
  else {
    out_r(_("t-Test with Tau not possible since Tau = 0\n") );
  }
  out_r(_("\nChi^2-test:\n") );
  out_r(_("Hypothesis H0: r=0  vs. H1: r#0\n") );
  alpha = get_chi_int(chisq, (np-2));
  out_r(_("Probability of H0: %6.4f\n\n"), 1.0-alpha);

  if (!noplot) {
    dose0 = pow(10., ((2.1-a0)/a1)); /* ~ 0% Probability */
    dose1 = pow(10., ((8.5-a0)/a1)); /* ~ 100 % Probability */
#ifndef STATIST_X
    out_r( _("doselab=|%s|, effectlab=|%s|\n"),
	get_name(dose), get_name(effect));
#endif
    plot_probit(dose, num, effect, n, a0, a1, dose0, dose1,
		get_label(dose), get_label(effect));
  }
  out_end();

}

/* =================================================================== */


void kolmo_test(REAL x[], int n) {
  REAL  mean, sdv;
  int   i, k, alpha;
  REAL  d, norm0, *cdf, *z, diff, fn1 = 0, fn2 = 0, z0 = 0;

  mean = get_mean(x,n);
  sdv = get_sdv(x,n);


  z =   (REAL*)m_calloc(n, sizeof(REAL));
  cdf = (REAL*)m_calloc(n, sizeof(REAL));
  for (i=0; i<n; i++) {
    DIV(z[i], (x[i]-mean), sdv);
    cdf[i] = (REAL)(i+1)/(REAL)n;
  }

  qsort(z, n, sizeof(REAL), real_compar_up);

  d = 0.0;
  for (i=0; i<n; i++) {
    norm0 = get_norm_int(z[i]);

    k = i+1;
    do {
      k--;
      diff = fabs(cdf[k]-norm0);
/*      printf("z=%g cdf=%g  norm=%g  d=%g\n", z[k], cdf[k], norm0, diff);  */
      if (diff > d) {
	d = diff;
	fn1 = cdf[k];
	fn2 = norm0;
	z0 = z[i];
      }
    }
    while (equal(z[k], z[i]));
  }

  if (!noplot) {
    plot_cdf_ks(z, n, z0, fn1, fn2, x, get_label(x));
  }

  alpha = pks(d, n);

  /* prob=pks(sqrt( SQR((REAL)n) / ((REAL)(2*n)))*(d));    */

  out_start();
  out_r(_("Hypothesis H0: Data is normally distributed versus\n") );
  out_r(_("Hypothesis H1: Data is not normally distributed\n\n") );

  colorize(ClHeader);
  out_r(_("Result KS-Lilliefors-Test on normal distribution:\n") );
  colorize(ClDefault);
  out_r(_("D (calculated)= %f\n"), d);
  out_r(_("Number of data = %i\n"), n);

  out_r(_("Mean = %g\n"), mean);
  out_r(_("Standard deviation = %g\n"), sdv);

  if (alpha > 0) {
    out_r(_("H0 accepted with a significance level of %i%%\n"), alpha);
  }
  else {
    out_r(_("H0 rejected\n") );
  }
  out_end();
}

void percentiles(REAL x[], int n) {
  REAL  mean, sdv;
  int   i, index;
  REAL  p, *z, alpha = 0;

  mean = get_mean(x,n);
  sdv = get_sdv(x,n);


  z =   (REAL*)m_calloc(n, sizeof(REAL));
  for (i=0; i<n; i++) {
    z[i] = x[i];
  }

  qsort(z, n, sizeof(REAL), real_compar_up);

  if (!noplot) {
    plot_cdf(z, n, get_label(x));
  }

  out_start();
  out_r(_("Percentiles for column \"%s\"\n"), get_label(x));
  alpha = 0;
  for(i = 0; i < 9; i++) {
	alpha += 0.1;
 	index = (int)(ceil((REAL)n * alpha));
 	p = z[index-1];
	out_r("%i%%\t%g\n", (int)(alpha*100.5), p);
  }
  alpha = 0.95;
  index = (int)(ceil((REAL)n * alpha));
  p = z[index-1];
  out_r("%d%%\t%g\n", (int)(alpha*100.0), p);
  alpha = 1.0;
  index = (int)(ceil((REAL)n * alpha));
  p = z[index-1];
  out_r("%d%%\t%g\n\n", (int)(alpha*100.0), p);

  out_end();
}


/* ===================================================================

void compare_dist(REAL x[], int nx, REAL y[], int ny) {
  REAL d, prob, alpha;

  ks(x, nx, y, ny, &d, &prob);

  out_start();
  out_r("Hypothese H0: x und y sind gleich verteilt gegen  \n");
  out_r("Hypothese H1: x und y sind nicht gleich verteilt  \n\n");

  out_r("Ergebnis Kolmogoroff-Smirnoff-Test auf Normalverteilung:\n");
  out_r("D (berechnet)= %f\n", d);
  out_r("Anzahl der x-Werte = %i\n", nx);
  out_r("Anzahl der y-Werte = %i\n", ny);

  if (prob == 0.0) {
    alpha = 0.0;
  }
  else {
    alpha = 1.0 - prob;
  }
  out_r("Irrtumswahrscheinlichkeit alpha fuer H0 = %f\n\n", alpha);
  out_end();
}

  =================================================================== */


void equal_freq(REAL x[], int n) {
  int *y, class[MCLASS], freq[MCLASS];
  int nclass, i, k, df;
  REAL chisq, phi, alpha, chi;
  BOOLEAN found;

  nclass = 0;
  y = (int*) m_calloc(n, sizeof(int));
  for (i=0; i<n; i++) {
    y[i] = get_round(x[i]);
  }

  for (i=0; i<n; i++) {
    found = FALSE;
    for (k=0; k<nclass; k++) {
      if (y[i]==class[k]) {
	found = TRUE;
	freq[k] ++;
	break;
      }
    }
    if (!found) {
      class[nclass] = y[i];
      freq[nclass] = 1;
      nclass ++;
    }
    if(nclass > MCLASS){
      out_err(MAT, ERR_FILE, ERR_LINE, 
	  _("Test aborted: There are more than %i classes of frequency."),
	  MCLASS);
      return;
    }
  }

  out_start();
  for (k=0; k<nclass; k++) {
    if (freq[k]<=5) {
      out_r(_("Warning: This test shouldn't be applied,\n"
	    "\tsince there are frequencies <= 5!\n\n") );
      break;
    }
  }


  chisq = 0;
  phi = (REAL)n/(REAL)nclass;

  if ( (nclass==2) && (n<200) ) {
    out_r(_("Correction according to YATES applied, since just 2 classes and n<200\n\n") );
    if (n<25) {
      out_r(_("Warning: FISCHER-Test shouldn't be applied,\n\tsince number of values <25\n\n") );
    }
    for (k=0; k<nclass; k++) {
      /* chisq += SQR(fabs(freq[k]-phi)-0.5)/phi; */
      DIV(chi, SQR(fabs(freq[k]-phi)-0.5), phi);
      chisq += chi;
    }
  }
  else {
    for (k=0; k<nclass; k++) {
      /* chisq += SQR( (REAL)freq[k] - phi)/phi; */
      DIV(chi, SQR( (REAL)freq[k] - phi), phi);
      chisq += chi;
    }
  }

  df = nclass -1;

  colorize(ClHeader);
  out_r(_("Result Chi^2-Test of equal frequency:\n") );
  colorize(ClDefault);
  out_r(_("Hypothesis H0: Values have equal frequency\n") );
  out_r(_("Hypothesis H1: Values don't have equal frequencies\n\n") );

  if (df >= 1) {
    out_r( _("Chi^2 = %f\n"), chisq);
    out_r(_("Degrees of freedom = %i\n"), df);
    alpha = 1. - get_chi_int(chisq, df);
    out_r(_("Probability of H0: %6.4f\n\n"), 1.0-alpha);
  }
  else {
    out_r(_("Chi^2-Test of normal distribution not possible since degrees of freedom < 1!\n\n") );
  }
  out_end();
}

/* ========================================================================== */


void compare_freq(REAL x[], int nx, REAL y[], int ny) {
  typedef struct {
    int class, x, y;
  } FREQUENCIES;

  int *ix, *iy;
  int nclass, i, k, df;
  REAL chisq, phi, alpha, chi;
  BOOLEAN found;
  FREQUENCIES freq[MCLASS];

  nclass = 0;
  for (i=0; i<MCLASS; i++) {
    freq[i].x = 0;
    freq[i].y = 0;
  }

  ix = (int*) m_calloc(nx, sizeof(int));
  for (i=0; i<nx; i++) {
    ix[i] = get_round(x[i]);
  }

  for (i=0; i<nx; i++) {
    found = FALSE;
    for (k=0; k<nclass; k++) {
      if (ix[i]==freq[k].class) {
	found = TRUE;
	freq[k].x ++;
	break;
      }
    }
    if (!found) {
      freq[nclass].class = ix[i];
      freq[nclass].x = 1;
      nclass ++;
    }
  }


  iy = (int*) m_calloc(ny, sizeof(int));
  for (i=0; i<ny; i++) {
    iy[i] = get_round(y[i]);
  }

  for (i=0; i<ny; i++) {
    found = FALSE;
    for (k=0; k<nclass; k++) {
      if (iy[i]==freq[k].class) {
	found = TRUE;
	freq[k].y ++;
	break;
      }
    }
    if (!found) {
      freq[nclass].class = iy[i];
      freq[nclass].y = 1;
      nclass ++;
    }
    if(nclass > MCLASS){
      out_err(MAT, ERR_FILE, ERR_LINE, 
	  _("Test aborted: There are more than %i classes of frequency."),
	  MCLASS);
      return;
    }
  }

  out_start();
  for (k=0; k<nclass; k++) {
    if (freq[k].x<=5) {
      out_r(_("Warning: This test shouldn't be applied,\n"
	    "\tsince there are frequencies <= 5!\n\n") );
      break;
    }
  }

  chisq = 0;

  if ( (nclass==2) && (nx<200) ) {
    out_r(_("Correction according to YATES applied, since just 2 classes and n<200\n\n") );
    if (nx<25) {
      out_r(_("Warning: FISCHER-Test shouldn't be applied,\n\tsince number of values <25\n\n") );
    }
    for (k=0; k<nclass; k++) {
      phi = ((REAL)freq[k].y/(REAL)ny)*(REAL)nx;
      /* chisq += SQR(fabs(freq[k].x-phi)-0.5)/phi; */
      DIV(chi, SQR(fabs(freq[k].x-phi)-0.5), phi);
      chisq += chi;
    }
  }
  else {
    for (k=0; k<nclass; k++) {
      phi = ((REAL)freq[k].y/(REAL)ny)*(REAL)nx;
      /* chisq += SQR( (REAL)freq[k].x - phi)/phi; */
      DIV(chi, SQR( (REAL)freq[k].x - phi), phi);
      chisq += chi;
    }
  }

  df = nclass -1;

  colorize(ClHeader);
  out_r(_("Result Chi^2-Test of equal frequency:\n") );
  colorize(ClDefault);
  out_r(_("Hypothesis H0: x and y are equally distributed\n") );
  out_r(_("Hypothesis H1: x and y are not equally distributed\n") );

  if (df >= 1) {
    out_r( _("Chi^2 = %f\n"), chisq);
    out_r(_("Degrees of freedom: %i\n"), df);
    alpha = 1. - get_chi_int(chisq, df);
    out_r(_("Probability of H0: %6.4f\n\n"), 1.0-alpha);
  }
  else {
    out_r(_("Chi^2-Test of normal distribution not possible since degrees of freedom < 1!\n\n") );
  }
  out_end();
}

/* ========================================================================== */



void lin_reg(REAL x[], REAL y[], int n) {
  REAL x_mean, x_sdv, y_mean, y_sdv;
  REAL s_xy, r_xy, a1, a0, rr, t, alpha;
  /*  REAL sse, qysum; */
  int  df;

  x_mean = get_mean(x,n);
  y_mean = get_mean(y,n);
  x_sdv = get_sdv(x,n);
  y_sdv = get_sdv(y,n);
  s_xy = (get_xysum(x,y,n) - ((get_sum(x,n)*get_sum(y,n))/(REAL) n))/
         (REAL) (n-1);
  DIV(r_xy, s_xy, x_sdv*y_sdv);
  rr = r_xy * r_xy;
  a1 = r_xy * (y_sdv/x_sdv);
  a0 = y_mean - a1*x_mean;

  df = n - 2;
  /* z = r_xy * sqrt((REAL)n-1.); */
/*
  sse = 0.0;
  qysum = 0.0;
  for (i=0; i<n; i++) {
    sse += SQR(y[i] - (a0 + a1*x[i]));
    qysum += SQR(y[i] - y_mean);
  }
*/
  out_start();
  colorize(ClHeader);
  out_r(_("\nResults of linear regression:\n") );
  colorize(ClDefault);
  out_r(_("number of data points n   : %i \n"),n);
  out_r(_("Intersection a0           : %g \n"),a0);
  out_r(_("Slope a1                  : %g \n"),a1);
  out_r(_("Correlation coefficient r : %g \n"),r_xy);
  out_r(_("Coefficient of determination r^2      : %g \n"),rr);
  out_r(_("Degr. of freedom df = n-2 : %i \n"), df);
/*  out_r("qysum=%f6.4, sse=%6.4f, r^2=%6.4f\n", qysum, sse, 1-sse/qysum); */
  if (fabs(r_xy) < 0.999999999) {
    t = r_xy * sqrt(((REAL)n-2.)/(1.-SQR(r_xy)));
    out_r(_("Estimated t-value         : %f\n"), t);
    alpha = get_t_int(fabs(t), df);
    out_r(_("\nt-Test:\n") );
    out_r(_("Hypothesis H0: r = 0  against hypothesis H1: r1 # 0  (->two-sided)\n") );
    /* a = alpha;  */
    out_r(_("Probability of H0 = %6.4f\n\n"), 1.0-alpha);
  } else{
    out_r(_("t-Test not possible since |r| = 1!\n") );
  }

  if (!noplot) {
    plot_pair(x, y, n, a0, a1, get_label(x), get_label(y));
  }

  out_r("\n");
  out_end();
}



/* =================================================================== */

void point_biserial_reg(REAL x[], REAL y[], int n) {
  REAL m_y0, m_y1, ysdv, *y0, *y1;
  REAL r_pb, rr, t, alpha;
  int *category;
  int i, n0=0, n1=0, df;

  y0 = (REAL*)m_calloc(n, sizeof(REAL));
  y1 = (REAL*)m_calloc(n, sizeof(REAL));
  category = (int*)m_calloc(n, sizeof(int));

  for (i=0; i<n; i++) {
    category[i] = get_round(x[i]);
    if ((category[i] != TRUE) && (category[i] != FALSE)) {
      out_err(ERR, ERR_FILE, ERR_LINE,
	  _("First column must contain only 0's and 1's!"));
      return;
    }
    if (category[i]) {
      y1[n1] = y[i];
      n1++;
    }
    else {
      y0[n0] = y[i];
      n0++;
    }
  }

  m_y0 = get_mean(y0, n0);
  m_y1 = get_mean(y1, n1);

  ysdv = get_sdv(y, n);
  DIV(r_pb, (m_y1-m_y0), ysdv);
  r_pb *= sqrt((REAL)(n1*n0)/(n*(n-1)));
  /* r_pb = (m_y1-m_y0)/ysdv*sqrt((REAL)(n1*n0)/(n*(n-1)));  */

  rr = r_pb*r_pb;
  df = n - 2;

  out_start();
  colorize(ClHeader);
  out_r(_("\nResult point biserial correlation:\n") );
  colorize(ClDefault);
  out_r(_("Number of data points n  : %i \n"), n);
  out_r(_("Correlation coefficient r: %20.12e \n"), r_pb);
  out_r(_("Coefficient of determination r^2     : %20.12e \n"), rr);
  out_r(_("Degrees of freedom df = n-2 : %i \n"), df);

  if (fabs(r_pb) < 1.) {
     t = r_pb * sqrt(((REAL)n-2.)/(1.-SQR(r_pb)));
     /* f = SQR(r_pb)/(1.-SQR(r_pb)) * (REAL)(n-2); */
    out_r(_("Calculated t-value      : %f \n"), t);
    alpha = get_t_int(fabs(t), df);
    out_r(_("\nt-Test:\n") );
    out_r(_("Hypothesis H0: r = 0  versus hypothesis H1: r1 # 0  (->two-sided)\n") );
    /* a = alpha; */
    out_r(_("Probability of H0 = %6.4f\n\n"), 1.0-alpha);
  } else{
    out_r(_("t-Test not possible since |r| = 1!\n") );
  }

  out_r("\n");
  out_end();
}


/* =================================================================== */


void histogram(REAL x[], int n, int mclass, REAL min, REAL max) {
			/* calculates the classes for a histrogram
                        calls: plot_hist() or plot_histo_discrete() with it*/

  int   aclass[MCLASS];	  /* counter for the values in class [i]	*/
  REAL  classval[MCLASS], fac, fac_short,arg, step;
  BOOLEAN zeroclass = FALSE;  /* true, if a class remains empty		*/
  BOOLEAN autoclass;	  /* whether autoclass is selected or not	*/
  BOOLEAN hit_max_mark;   /* TRUE, when a point hits the max mark	*/
  BOOLEAN discrete=FALSE; /* if we have only points, hitting the marks	*/
  int	lower=0,higher=0; /*counts, how many values are not counted in aclass*/
  int interval;		  /* width of one class				*/
  int i,ii;		  /* loop variables				*/

  autoclass=(mclass==0);	/* check for autoclass			*/

  if (autoclass) {  				/* autoclass was choosen*/
				/* start with a good approximation
				of the number of classes. Adjust that,
				as long, as we get empty classes	*/
    mclass = get_round(1. + 3.32*log10((REAL)n));
  }


  do {
    for (i=0; i<mclass; i++) {/* clear class counters			*/
      aclass[i] = 0;
    }
    lower=higher=0;		/* reset counters			*/
    hit_max_mark=FALSE;		/* resetting				*/

    fac = (max-min)/(REAL)(mclass);             /* size of each intervall*/
    for (i=0; i<n; i++) {
      if(x[i]<min) { lower++; continue; }     /* count lower values   */

      arg = (x[i]-min)/fac;
      interval = (int) arg;
      /* see if the value lower or equal to max
	 inc. aclass[intervall] or higher accordingly*/
      if(interval<mclass)
      {
	aclass[interval]++;
      }else
      {                         /* We have to include the right border for
                                cases with pure integers                */
	if(x[i]<=max)
	{
	  aclass[mclass-1]++;
				/* if a value hits the last point, there is
                                a chance that we are discrete           */
	  hit_max_mark=TRUE;
	}
	else
	  higher++;
      }
    }
				/* if we are in automatic class selection
                                check for for classes, which are empty
                                with the current number of classes	*/
    if(autoclass)
    {
      zeroclass = FALSE;
      for (i=0; i<mclass; i++) {
	if (aclass[i] == 0) {
	  zeroclass = TRUE;
	  mclass --;
	  break;

	}
      }
    }
  } while (autoclass && zeroclass);

				/* test if we are discrete, when a point has
				hit the max mark with the last mclass number*/
  if(hit_max_mark)
  {
    REAL help;
    discrete=TRUE;
    fac_short=(max-min)/(REAL)(mclass-1);

    for(ii=0; ii<n;ii++)
    {   help=(x[ii]-min)/fac_short;
      if( help!=(int)(help))
      { discrete=FALSE; break; }
    }
  }

					/* we have aclass[] properly filled*/
  out_start();
  if(discrete)
  {
    out_d(_("Discrete data with exactly %d classes!\n"), mclass);
    step = (max-min)/(REAL)(mclass-1);
    for (i=0; i<mclass; i++) {
      classval[i] = min + (step*(REAL)i);
    }
  } else
  {
    step = (max-min)/(REAL)(mclass);
    for (i=0; i<mclass; i++) {
      classval[i] = min + (step/2.) + (step*(REAL)i);
    }
  }

  if ((!thist) && (!noplot)){
    if(discrete)
      plot_histo_discrete(classval, aclass, mclass, step, x, get_label(x));
    else
      plot_histo(classval, aclass, mclass, step, x, get_label(x));
  }
  else {
    if( (!bernhard)|| ( bernhard && thist ))
    {
      out_r(_("\nHistogram (class center - number of values)\n") );
      print_histo(classval, aclass, mclass);
    }
  }
  out_end();
} /* end of histogram() */



void standard(REAL x[], int n, int mclass, REAL min, REAL max) {
  REAL mean, sdv, sdv_mean, var_coef, median;
  REAL *xh, t, conf;
  REAL q_l, q_u, index, no_u, no_l;
  int i, i_l, i_u;

  int   xcrit;

/* critical values for 95% confidence intervals of the Medians (from NEAVE):*/
  int critical[50] = {
    -1, -1, -1, -1, -1, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3,
    3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 8, 8, 9, 9, 9, 10, 10, 11,
    11, 12, 12, 12, 13, 13, 14, 14, 15, 15, 15, 16, 16, 17, 17
  };
  REAL med_conf_lo, med_conf_up;

/*  if (mclass==0) {
    min = get_min(x,n);
    max = get_max(x,n);
  }*/


  out_start();
  if (max==min) {
    if(!bernhard) {
      out_err(MAT, ERR_FILE, ERR_LINE, _("All values are equal!") );
    } else
    {
      colorize(ClHeader);
      out_r(_("#Result general statistical information in a table\n") );
      colorize(ClDefault);
      out_r(_("#All values are equal!\n") );
      out_r("#n\tmean\tm-conf\tm+conf\tmedian\tme_c_lo\tme_c_up"
	  "\tquar_lo\tquar_up\tsdv\tvarc(%%)\tsdv_err\tmin\tmax\n");
      out_r("%i\t",  n);
      out_r("%g\t%g\t%g\t", min, min, min);
      out_r("%g\t%g\t%g\t", min, min, min);
      out_r("%g\t", min);
      out_r("%g\t", min);
      out_r("%g\t", 0);
      out_r("%f\t", 0);
      out_r("%g\t", 0);
      out_r("%g\t", get_min(x,n));
      out_r("%g\t\n", get_max(x,n));

    }

    out_end();
    return;
  }
  if (min > max) {
    out_err(MAT, ERR_FILE, ERR_LINE,
	_("Minimum is larger than maximum!") );
    out_end();
    return;
  }

  mean = get_mean(x,n);
  sdv = get_sdv(x,n);
  sdv_mean = sdv/sqrt((REAL)n);
  var_coef = (sdv/mean) * 100.;
  t = get_t(0.975, (n-1));
  conf = (t*sdv)/sqrt((REAL)n);


  histogram(x,n,mclass,min,max);	/* calc and plot the histogram	*/


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

  index = (REAL)n*0.25;
  if (index==floor(index)) {
    i_l = (int)index-1;
    i_u = (int)index;
  }
  else {
    i_l = (int) floor(index-1.);
    i_u = (int) ceil(index-1.);
  }

  q_l = 0.5 * (xh[i_l]+xh[i_u]);
  q_u = 0.5 * (xh[n-i_l-1]+xh[n-i_u-1]);
  no_u = median + 1.58*(q_u-q_l)/sqrt((REAL)n);
  no_l = median - 1.58*(q_u-q_l)/sqrt((REAL)n);

  if (n<6) {
    med_conf_lo = med_conf_up = 0.0;
  }
  else if (n<=50) {
    med_conf_lo = xh[critical[n-1]];
    med_conf_up = xh[n-critical[n-1]-1];
  }
  else {
    xcrit = (int)floor(-(0.5*sqrt((REAL)n)) *1.96-0.5+0.5*(REAL)n);
    med_conf_lo = xh[xcrit];
    med_conf_up = xh[n-xcrit-1];
  }
    /*
       printf("crit(50)=%i\n", (int)floor(-(0.5*sqrt(50.0)) *1.96-0.5+0.5*50.0));
       */
  if(!bernhard)
  {
    colorize(ClHeader);
    out_r(_("\nResult general statistical information:\n") );
    colorize(ClDefault);
    out_r(_("Count                          : %i\n"),  n);
    out_r(_("Mean                           : %g  [%g, %g] (95%%)\n"),
	mean, mean-conf, mean+conf);
    out_r(_("Median                         : %g  [%g, %g] (95%%)\n"),
	median, med_conf_lo, med_conf_up);
  /*  out_r("                               :     [%g, %g] (90%%)\n",
      no_l, no_u);
      */
    out_r(_("25%% quartile                   : %g\n"), q_l);
    out_r(_("75%% quartile                   : %g\n"), q_u);
    out_r(_("Standard deviation             : %g\n"), sdv);
    out_r(_("Variation coefficient          : %f %%\n"), var_coef);
    out_r(_("Standard error of mean         : %g\n"), sdv_mean);
    out_r(_("Minimum                        : %g\n"), get_min(x,n));
    out_r(_("Maximum                        : %g\n\n"), get_max(x,n));
  }else
  {
    colorize(ClHeader);
    out_r(_("#Result general statistical information in a table\n") );
    colorize(ClDefault);
    out_r("#n\tmean\tm-conf\tm+conf\tmedian\tme_c_lo\tme_c_up"
	"\tquar_lo\tquar_up\tsdv\tvarc(%%)\tsdv_err\tmin\tmax\n");
    out_r("%i\t",  n);
    out_r("%g\t%g\t%g\t", mean, mean-conf, mean+conf);
    out_r("%g\t%g\t%g\t", median, med_conf_lo, med_conf_up);
    out_r("%g\t", q_l);
    out_r("%g\t", q_u);
    out_r("%g\t", sdv);
    out_r("%f\t", var_coef);
    out_r("%g\t", sdv_mean);
    out_r("%g\t", get_min(x,n));
    out_r("%g\t\n", get_max(x,n));

  }
  out_end();

}

/* =================================================================== */


void rank_order(REAL x[], REAL y[], int n) {
  REAL  rho, z, alpha, a;
  int   i, df;
/*  BOOLEAN infinity = TRUE; */

/* critical raw-Values for significance levels 1%, 2%, 5% and 10% for n <=30
   indexes (two-tailed), is valid only for n >= 5! (from: SACHS, p. 230)
  */

  const float sig[31][4] = {
         {0.0000,0.0000,0.0000,0.0000}, {0.0000,0.0000,0.0000,0.0000},
	 {0.0000,0.0000,0.0000,0.0000}, {0.0000,0.0000,0.0000,0.0000},
         {0.0000,0.0000,0.0000,0.0000}, {0.9001,0.9000,0.9000,0.8000},
	 {0.9429,0.8857,0.8286,0.7714}, {0.8929,0.8571,0.7450,0.6786},
         {0.8571,0.8095,0.6905,0.5952}, {0.8167,0.7667,0.6833,0.5833},
	 {0.7818,0.7333,0.6364,0.5515}, {0.7454,0.7000,0.6091,0.5273},
	 {0.7273,0.6713,0.5804,0.4965}, {0.6978,0.6429,0.5549,0.4780},
	 {0.6474,0.6220,0.5341,0.4593}, {0.6536,0.6000,0.5179,0.4429},
	 {0.6324,0.5824,0.5000,0.4265}, {0.6152,0.5637,0.4853,0.4118},
	 {0.5975,0.5480,0.3994,0.3994}, {0.5825,0.5333,0.4579,0.3895},
	 {0.5684,0.5203,0.4451,0.3789}, {0.5545,0.5078,0.4351,0.3688},
	 {0.5426,0.4963,0.4241,0.3597}, {0.5306,0.4852,0.4150,0.3518},
	 {0.5200,0.4748,0.4061,0.3435}, {0.5100,0.4654,0.3977,0.3362},
	 {0.5002,0.4564,0.3894,0.3299}, {0.4915,0.4481,0.3822,0.3236},
	 {0.4828,0.4401,0.3749,0.3175}, {0.4744,0.4320,0.3685,0.3113},
	 {0.4665,0.4251,0.3620,0.3059}
       };


  if ((rho=get_rank_correlation(x, y, n))==REAL_MIN) {
    return;
  }

  out_start();
  df = n-2;
  colorize(ClHeader);
  out_r(_("\nResult SPEARMAN's Rank-Correlation:\n") );
  colorize(ClDefault);
  out_r(_("rho = %f\n"), rho);
  out_r(_("Degrees of freedom = n-2 = %i\n\n"), df);

  out_r(_("Hypothesis H0: rho=0 versus hypothesis H1: rho#0 (->two-sided)\n"));

  if (n<5) {
    out_r(_("Test not possible since n<5 (too few values!)\n\n") );
    out_end();
    return;
  }
  else if ((n>=5) && (n<=30)) {
    for (i=0; i<4; i++) {
      if (fabs(rho) > (REAL)sig[n][i]) {
	break;
      }
    }
    switch (i) {
    case 0: alpha = 0.01;
      break;
    case 1: alpha = 0.02;
      break;
    case 2: alpha = 0.05;
      break;
    case 3: alpha = 0.10;
      break;
    default: alpha = 1.0;
      break;
    }

    if (alpha < 1.0) {
      out_r(_("H0 rejected at a level of significance of %4.2f\n\n"),
	    alpha);
    }
    else {
      out_r(_("Alpha > 0.1 ==> H0 can not be rejected\n\n") );
    }
  }
  else {
/*
    a = get_norm_int(2.326);
    alpha =  1.- ((1.-a)*2.0);
    out_d("norm(2.326)=%f\n", alpha);
*/
    z = fabs(rho) * sqrt((REAL)n-1.0);
    out_r(_("Significance checked using the normal distribution\n") );
    out_d("z = %f\n", z);
    a = get_norm_int(z);
    alpha =  1.- ((1.-a)*2.0);
    out_r(_("Probability of H0 = %6.4f\n\n"), 1.0-alpha);
  }
  out_end();
}

/* =================================================================== */

void t_test(REAL x[], int nx, REAL y[], int ny) {
  REAL denom, t, x_mean, y_mean;
  REAL x_sum, y_sum, x_qsum, y_qsum, alpha, a;
  REAL x_qmean, y_qmean;
  int  df;
/* t-Test for difference between two means of two samples */

  x_sum = get_sum(x,nx);
  y_sum = get_sum(y,ny);
  x_qsum = get_qsum(x,nx);
  y_qsum = get_qsum(y,ny);

  x_mean = x_sum/(REAL)nx;
  y_mean = y_sum/(REAL)ny;


  DIV(x_qmean, SQR(x_sum), (REAL)nx);
  DIV(y_qmean,  SQR(y_sum), (REAL)ny);
  DIV(denom, (x_qsum - x_qmean + y_qsum - y_qmean), (REAL)(nx + ny - 2));
  denom *= (1./nx + 1./ny);
  DIV(t, fabs(x_mean - y_mean), sqrt(denom));
  df = nx + ny -2;

  out_start();
  colorize(ClHeader);
  out_r(_("\nResult t-Test for independent samples\n") );
  colorize(ClDefault);
  out_r(_("Degrees of freedom = n1 + n2 - 2 = %i\n"), df);
  if (t !=0 ) {
    alpha = get_t_int(fabs(t), df);
    out_r("t = %f\n", t);
    out_r(_("\nHypothesis H0: mu1=mu2 versus hypothesis H1: mu1#mu2 (->two-sided)\n") );
    a = alpha;
    out_r(_("Probability of H0 = %6.4f\n\n"), 1.0-a);
  }
  else {
    out_r(_("t-Test not possible since t = 0!\n") );
  }

/*
  out_r("Hypothese H0: mu1=mu2  gegen  Hypothese H1: mu1>mu2 (->einseitig)\n");
  sgn = get_sgn(t);
  a = 1.- 0.5*(1. + alpha*sgn);
  out_r("Wahrscheinlichkeit fuer H0 = %6.4f\n\n", 1.0-a);

  out_r("Hypothese H0: mu1=mu2  gegen  Hypothese H1: mu1<mu2 (->einseitig)\n");
  a = 1.- 0.5*(1. - alpha*sgn);
  out_r("Wahrscheinlichkeit fuer H0 = %6.4f\n\n", 1.0-a);
*/
  out_end();
}

/* ====================================================================== */



void pair_t_test(REAL x[], REAL y[], int n) {
  REAL t, alpha, a, dif_sdv, dif_mean, *dif;
  int  df, i;
/* t-Test for difference between two means of two paired samples */

  dif = (REAL*)m_calloc(n, sizeof(REAL));
  for (i=0; i<n; i++) {
    dif[i] = x[i]-y[i];
  }
  dif_sdv = get_sdv(dif, n);
  dif_mean = get_mean(dif, n);
  DIV(t, dif_mean*sqrt((REAL)n), sqrt(SQR(dif_sdv)));
  t = fabs(t);
  /*  t = fabs( (dif_mean*sqrt((REAL)n))/(sqrt(SQR(dif_sdv))) ); */
  df = n-1;

  out_start();
  colorize(ClHeader);
  out_r(_("\nResult t-Test for pairwise ordered samples\n") );
  colorize(ClDefault);
  out_r(_("Degrees of freedom n-1 = %i\n"), df);
  if (t !=0 ) {
    alpha = get_t_int(fabs(t), df);
    out_r("t = %f\n", t);
    out_r(_("\nHypothesis H0: mu1=mu2 versus hypothesis H1: mu1#mu2 (->two-sided)\n") );
    a = alpha;
    out_r(_("Probability of H0 = %6.4f\n\n"), 1.0-a);
  }
  else {
    out_r(_("t-Test not possible since t = 0!\n") );
  }

/*
  out_r("Hypothese H0: mu1=mu2  gegen  Hypothese H1: mu1>mu2 (->einseitig)\n");
  sgn = get_sgn(t);
  a = 1.- 0.5*(1. + alpha*sgn);
  out_r("Wahrscheinlichkeit fuer H0 = %6.4f\n\n", 1.0-a);

  out_r("Hypothese H0: mu1=mu2  gegen  Hypothese H1: mu1<mu2 (->einseitig)\n");
  a = 1.- 0.5*(1. - alpha*sgn);
  out_r("Wahrscheinlichkeit fuer H0 = %6.4f\n\n", 1.0-a);
*/
  out_end();
}

/* ====================================================================== */



void u_test(REAL x[], int nx, REAL y[], int ny) {

  SORTREC *vrec;
  int  i, k, j, n;
  REAL m, ux, uy, rx, ry, umin, z, alpha, t, denom, arg;

  n = nx+ny;

  vrec = (SORTREC*)m_calloc(n, sizeof(SORTREC));

  for (i=0; i<nx; i++) {
    vrec[i].val = x[i];
    vrec[i].ind = 0;
  }

  for (i=nx; i<n; i++) {
    vrec[i].val = y[i-nx];
    vrec[i].ind = 1;
  }

  qsort(vrec, n, sizeof(SORTREC), u_rank_compar);

  i = 0;
  k = 0;
  m = 0.;
  t = 0.;
  while (i<n) {
    vrec[i].rank = (REAL)i+1.;
    if ( (i<(n-1)) && (equal(vrec[i].val, vrec[i+1].val)) ) {
      k++;
      m += (REAL)i;
    }
    else {
      if (k!=0) {
        m += (REAL)i;
        k ++;
        /* t += (pow((REAL)k, 3.) - (REAL)k)/12.; */
        t +=  ( (REAL)k * (SQR((REAL)k) -1.) ) /12.;
        m = m/ (REAL)(k);
        for (j=i; j>(i-k); j--) {
          vrec[j].rank = m+1.;
        }
      }
      k=0;
      m=0.;
    }
    i++;
  }

  rx = 0.0;
  ry = 0.0;

  for (i=0; i<n; i++) {
    if (vrec[i].ind == 0) {
      rx += vrec[i].rank;
    }
    else {
      ry += vrec[i].rank;
    }
/*
    out_r("Nr %i  ind %i  Wert %f  Rank %f\n",
           i, vrec[i].ind, vrec[i].val, vrec[i].rank);
*/
  }

  ux = (REAL)nx*(REAL)ny + ( (REAL)nx*(REAL)(nx+1) )/2. - rx;
  uy = (REAL)nx*(REAL)ny + ( (REAL)ny*(REAL)(ny+1) )/2. - ry;

  umin = (ux<uy) ? ux : uy;
/*  z = fabs(umin-((REAL)nx*(REAL)ny)/2.) /
    sqrt( ( (REAL)nx*(REAL)ny*( (REAL)nx+(REAL)ny+1.) )/12.);
*/
  arg = (REAL)nx*(REAL)ny/((REAL)n*(REAL)(n-1))*((REAL)n*(SQR((REAL)n)-1.)/12. -t);
  SQRT(denom, arg);
  DIV(z, fabs(umin-((REAL)nx*(REAL)ny)/2.), denom);
/*
  z = fabs(umin-((REAL)nx*(REAL)ny)/2.) /
    sqrt( ( (REAL)nx*(REAL)ny/((REAL)n*(REAL)(n-1)) *
           ( (REAL)n * (SQR((REAL)n)-1.)/12. -t ) )  );
*/
  out_start();
  colorize(ClHeader);
  out_r(_("\nResult u-Test:\n") );
  colorize(ClDefault);
  out_r("Rx = %f   Ry = %f\n", rx, ry);
  out_r("Ux = %f   Uy = %f\n", ux, uy);
  out_r("nx = %i   ny = %i\n", nx, ny);
  out_r("U = %f\n", umin);

  out_r(_("\nHypothesis H0: x and y originate from the same set versus\n") );
  /* Test on x > y if ux < uy. */
  if( ux < uy ) {
	out_r(_("Hypothesis H1: x stochastically larger than y (-> one-sided test) or\n") );

  } else { /* ux >= uy */
	out_r(_("Hypothesis H1: x stochastically smaller than y (-> one-sided test) or\n") );
  }
  out_r(_("              x is different from y (-> two-sided test)\n\n") );


/*  if ( ((nx>=8) && (ny>=8)) || (nx>20) || (ny>20) ) {  */
  if ( ((nx<7) || (ny<7)) || (nx<=19) || (ny<=19) ) {
    out_r(_("Warning: Only rough approximation of probability possible!\n") );
    out_r(_("Please check exact probability of alpha for h having %i "
	  "degrees of freedom\n"), ncol-1);
    out_r(_("in the literature, e.g. in Table 16/17, pp. 599 in WEBER \n\n") );
  }

  out_r(_("Normally distributed random variable   z = %f\n"), z);
  out_r(_("Correction term for equal ranks t = %f\n"), t);
  alpha = get_norm_int(z);
  out_r(_("Probability of H0 (one-sided) = %6.4f\n"), 1.0-alpha);
  alpha = 1.- ((1.-alpha)*2.);
  out_r(_("Probability of H0 (two-sided) = %6.4f\n\n"), 1.0-alpha);
  out_end();
}


/* ====================================================================== */


void kruskal_test(PREAL x[], int nrow[], int ncol) {
  SORTREC *vrec;
  unsigned int n=0, ir, ic, k=0, i, j;
  REAL *r, nr, t, m, cor, prob, h=0.;

  for (ic=0; ic<ncol; ic++) {
    n += nrow[ic];
  }

  vrec = (SORTREC*)m_calloc(n, sizeof(SORTREC));
  r = (REAL*)m_calloc(ncol, sizeof(REAL));

  for (ic=0; ic<ncol; ic++) {
    r[ic] = 0.;
    for (ir=0; ir<nrow[ic]; ir++) {
      vrec[k].val = x[ic][ir];
      vrec[k].ind = ic;
      k++;
    }
  }
  out_start();
#ifdef DEBUG
  if (n != k) {
    out_r(_("Error: n!=k\n") );
  }
#endif

  qsort(vrec, n, sizeof(SORTREC), u_rank_compar);

  i = 0;
  k = 0;
  m = 0.;
  t = 0.;
  while (i<n) {
    vrec[i].rank = (REAL)i+1.;
    if ( (i<(n-1)) && (equal(vrec[i].val, vrec[i+1].val)) ) {
      k++;
      m += (REAL)i;
    }
    else {
      if (k!=0) {
        m += (REAL)i;
        k ++;
        t += (REAL)k * (SQR((REAL)k)-1.);
        m = m/ (REAL)(k);
        for (j=i; j>(i-k); j--) {
          vrec[j].rank = m+1.;
        }
      }
      k=0;
      m=0.;
    }
    i++;
  }

  for (i=0; i<n; i++) {
    r[vrec[i].ind] += vrec[i].rank;
  }


  nr = (REAL)n;
  cor = 1. - t/(SQR(nr)*(nr-1.));

  for (ic=0; ic<ncol; ic++) {
    h += (SQR(r[ic])/(REAL)nrow[ic]);
    /* out_r("r[%u] = %g\n", ic, r[ic]);  */
  }

  h = (12./(nr*(nr+1.)) * h - 3.*(nr+1.))/cor;

  colorize(ClHeader);
  out_r(_("\nResult Kruskal-Wallis-Test:\n") );
  colorize(ClDefault);
  out_r(_("Chi^2-distributed random variable H = %g\n"), h);
  out_r(_("Correction term for equal ranks (ties) = %g\n"), cor);
  out_r(_("Degrees of freedom = %i\n"), ncol-1);

  out_r(_("\nHypothesis H0: Samples originate from the same set versus\n") );
  out_r(_("Hypothesis H1: Samples do not originate from the same set\n") );

 /* if ((ncol>=3) && (nrow[0]>=5) && (nrow[1]>=5) && nrow[2]>=5) {  */
  if ((ncol<2) || (nrow[0]<4) || (nrow[1]<4) || nrow[2]<4) {
    out_r(_("Warning: Only rough approximation of probability possible!\n") );
    out_r(_("Please check exact probability of alpha for h having %i degrees "
	  "of freedom\n"), ncol-1);
    out_r(_("in the literature, e.g. in Table 16/17, pp. 599 in WEBER \n\n") );
  }
  if (h > 0.0) {
    prob = get_chi_int(h, (ncol-1));
    out_r(_("Probability alpha for H0 = %6.4f\n\n"), prob);
  }
  else {
    out_err(MAT, ERR_FILE, ERR_LINE,
    	_("Calculation of Chi^2-distribution not possible since h<0!") );
  }
  out_end();
}


/* ====================================================================== */

void vierfeld_test(REAL x[], REAL y[], int n) {
   const REAL  pi=3.14159265358979323846;
   long unsigned int a=0, b=0, c=0, d=0, i, xi, yi;
   REAL chi, r, T, prob, alpha, delta, beta, gamma;

  out_start();
   if (n!=2) {
     out_r(_("Characteristics are counted (1='existent', 0='not existent')\n\n"));
     for (i=0; i<n; i++) {
       xi = get_round(x[i]);
       yi = get_round(y[i]);
       if ( ((xi!=1) && (xi!=0)) || ((yi!=1) && (yi!=0)) ) {
         out_err(ERR, ERR_FILE, ERR_LINE,
	     _("Columns must contain only 0's and 1's!") );
	   out_end();
         return;
       }
       if ( xi &&  yi) { a++; }
       if (!xi &&  yi) { b++; }
       if ( xi && !yi) { c++; }
       if (!xi && !yi) { d++; }
     }
   }
   else {
     out_r(_("Values are being interpreted as fourfold-table:\n\n") );
     a = (unsigned int) x[0];
     b = (unsigned int) y[0];
     c = (unsigned int) x[1];
     d = (unsigned int) y[1];
     n = a+b+c+d;
   }

   out_r( _("Fourfold-table:\n\n") );
   out_r( _(" |                   |  A existent   |  A not existent     |\n"));
   out_r(" +-------------------+---------------+---------------------+\n");
   out_r( _(" | B existent        |      a        |         b           |\n"));
   out_r( _(" | B not existent    |      c        |         d           |\n"));
   out_r(" +-------------------+---------------+---------------------+\n\n");

   alpha = ((REAL)(a+b) * (REAL)(a+c))/(REAL)n;
   delta = ((REAL)(d+b) * (REAL)(d+c))/(REAL)n;
   beta  = ((REAL)(a+b) * (REAL)(b+d))/(REAL)n;
   gamma = ((REAL)(c+d) * (REAL)(a+c))/(REAL)n;

   if ((alpha>=5.) && (delta>=5.) && (beta>=5.) && (gamma>=5.)) {
     chi = (SQR( (REAL)a*(REAL)d - (REAL)b*(REAL)c  )  *(REAL)n) /
       ( (REAL)(a+b) * (REAL)(c+d) * (REAL)(a+c) * (REAL)(b+d) );
   }
   else {
    out_r(_("Using according to YATES corrected Chi^2 value, since frequencies <= 5 expected!\n\n") );

    chi = (SQR(fabs( (REAL)a*(REAL)d - (REAL)b*(REAL)c) -0.5*(REAL)n)*(REAL)n)/
      ((REAL)(a+b) * (REAL)(c+d) * (REAL)(a+c) * (REAL)(b+d));
  }
  T = (((REAL)a*(REAL)d - (REAL)b*(REAL)c))/
    (sqrt( (REAL)(a+b) * (REAL)(c+d) * (REAL)(a+c) * (REAL)(b+d) ));
  r = sin(T*(pi/4.));

  out_r(_("observed: a=%lu,  b=%lu,  c=%lu,  d=%lu,  n=%i\n"),
		a, b, c, d, n);
  out_r(_("expected: a=%4.2f,  b=%4.2f,  c=%4.2f,  d=%4.2f,  n=%i\n"),
         alpha, beta, gamma, delta, n);
  out_r( _("Chi^2 = %f\n"), chi);
  out_r(_("Contingency coefficient r (according to HAMMING) = %f\n\n"), r);
  out_r(_("Chi^2-test:\n") );
  out_r(_("Hypothesis H0: Both items are independent (uncorrelated)\n") );
  out_r(_("versus H1: Both items are dependent (correlated)\n") );
  prob = get_chi_int(chi, 1);
  out_r(_("Probability of H0: %6.4f\n\n"), prob);
  out_end();
}


/* ====================================================================== */


void tafel_test(PREAL xx[], int nrow, int ncol) {
  REAL fieldsum=0.0, sum=0.0, *colsum, *rowsum, lncols=0.0, lnrows=0.0,
       prob, g, fsum, chi=0.0;
  int i, k, df, max_err = 100, zero_cols=0, zero_rows=0;
  PREAL *e;
  BOOLEAN freq_too_small=FALSE, printall=TRUE;

  colsum = (REAL*)m_calloc(ncol, sizeof(REAL));
  rowsum = (REAL*)m_calloc(nrow, sizeof(REAL));
  e = (PREAL*) m_calloc(ncol, sizeof(PREAL));
  for (i=0; i<ncol; i++) {
    e[i] = (REAL*)m_calloc(nrow, sizeof(REAL));
  }

  out_start();
  for (i=0; i<ncol; i++) {
    colsum[i] = 0.0;
    for (k=0; k<nrow; k++) {
      if(xx[i][k] < 0.0){
	out_err(MAT, ERR_FILE, ERR_LINE,
	    _("Column \"%s\", line %i, has value < 0"),
	    get_name(xx[i]), k);
	return;
      }
      fieldsum += xx[i][k] * get_ln_0(xx[i][k]);
      sum += xx[i][k];
      colsum[i] += xx[i][k];
    }
    if (colsum[i] == 0.0) {
      zero_cols ++;
    }
    out_r("\n");
  }

  for (k=0; k<nrow; k++) {
    rowsum[k] = 0.0;
    for (i=0; i<ncol; i++) {
      rowsum[k] += xx[i][k];
    }
    if (rowsum[k] == 0.0) {
      zero_rows ++;
    }
  }

  for (k=0; k<nrow; k++) {
    for (i=0; i<ncol; i++) {
      e[i][k] =(colsum[i]*rowsum[k])/sum;
      if (e[i][k] != 0.0) {
        chi += SQR(xx[i][k]-e[i][k])/e[i][k];
      }
      if (e[i][k] < 5.) {
        freq_too_small = TRUE;
      }
    }
  }


  out_r(_("Analysis of two-items table:\n\n") );

  for (i=0; i<ncol; i++) {
    if(colsum[i] < 0.0){
      out_err(MAT, ERR_FILE, ERR_LINE,
	  _("Sum of values of column \"%s\" is lower than 0"),
	  get_name(xx[i]));
      return;
    }
    lncols += colsum[i] * get_ln_0(colsum[i]);
  }

  for (k=0; k<nrow; k++) {
    if(rowsum[k] < 0.0){
      out_err(MAT, ERR_FILE, ERR_LINE,
	  _("Sum of values of row %i is lower than 0"), k);
      return;
    }
    lnrows += rowsum[k] * get_ln_0(rowsum[k]);
  }

  if(sum < 0.0){
    out_err(MAT, ERR_FILE, ERR_LINE,
	_("Sum of all values of all columns is lower than 0"));
    return;
  }
  fsum = sum * get_ln_0(sum);
  g = 2 * (fieldsum-lncols-lnrows+fsum);
  df = (ncol-1-zero_cols) * (nrow-1-zero_rows);

/*
  out_r("Transformierte Zeilensummen = %f\n", lnrows);
  out_r("Transformierte Spaltensummen = %f\n", lncols);
  out_r("Transformierte Haeufigkeiten = %f\n", fieldsum);
  out_r("Transformierte Gesamtsumme = %f\n", fsum);

*/


  out_r(  _(" Item A   |                   Item B                 \n") );
  out_r("          |");
  for (i=0; i<ncol; i++) {
    out_r(_("Class%3i  |"), (i+1));
  }
  out_r(_("  Sum     |\n") );
  for (i=-2; i<ncol; i++) {
    out_r("----------+");
  }
  out_r("\n");


  for (k=0; k<nrow; k++) {
    if(k == max_err){
      printall = FALSE;
      out_i(_("There are %i rows to be printed yet. "
	    "Do you want to see them? (%s) "), (nrow - max_err), _("y/N"));
      GETNLINE;
      if(!(empty) && (line[0] ==_("y")[0] || line[0] != _("Y")[0])){
	printall = TRUE;
      }
    }
    if(printall){
      out_r(_("Class%3i  |"), k+1);
      for (i=0; i<ncol; i++) {
	out_r("%7u   |", (unsigned int)xx[i][k]);
      }
      out_r("%7i   |\n", (unsigned int)rowsum[k]);
      out_r(_("exp. frq. |") );
      for (i=0; i<ncol; i++) {
	out_r("%7u   |", (unsigned int)e[i][k]);
      }
      out_r("          |\n");
      if (k<(nrow-1)) {
	for (i=-2; i<ncol; i++) {
	  out_r("----------+");
	}
	out_r("\n");
      }
    }
  }

  for (i=-2; i<ncol; i++) {
    out_r("----------+");
  }
  out_r(_("\n Sum      |") );
  for (i=0; i<ncol; i++) {
    out_r("%7i   |", (unsigned int)colsum[i]);
  }
  out_r("%7i   |\n", (unsigned int)sum);

  out_r( _("\nChi^2                            = %f\n"), chi);
  out_r(_("G (check value for Chi^2-Test) = %f\n"), g);
  out_r(_("Degrees of freedom = %i\n\n"), df);
  if (freq_too_small) {
    out_r(_("Warning: Expected frequencies < 5!\n\n") );
  }

  out_r(_("Chi^2-test:\n") );
  out_r(_("Hypothesis H0: Both items are independent (uncorrelated)\n") );
  out_r(_("versus H1: Both items are dependent (correlated)\n") );
  prob = get_chi_int(g, df);
  out_r(_("Probability of H0: %6.4f\n\n"), prob);
  out_end();

}


/* ====================================================================== */

void wilcoxon_test(REAL x[], REAL y[], int n) {
  SORTREC *vrec;
  int i, k, j, nv, alpha, navg, t;
  REAL *diff, *avg;
  REAL m, s_plus, s_minus, s_min, d, nr, nd, prob, z, median, conf_u, conf_l;
  const short int sig[20][3] = {
    {  0, -1, -1 }, {  2,  0, -1 }, {  4,  2,  0 }, {  6,  3,  2 },
    {  8,  5,  3 }, { 11,  7,  5 }, { 14, 10,  7 }, { 17, 13, 10 },
    { 21, 16, 13 }, { 25, 20, 16 }, { 30, 24, 20 }, { 35, 28, 23 },
    { 40, 33, 28 }, { 46, 38, 32 }, { 52, 43, 38 }, { 59, 49, 43 },
    { 66, 56, 49 }, { 73, 62, 55 }, { 81, 69, 61 },
    { 89, 77, 68 }
  };


  diff = (REAL*) m_calloc(n, sizeof(REAL));
  vrec = (SORTREC*)m_calloc(n, sizeof(SORTREC));

  nv = 0;
  for (i=0; i<n; i++) {
    d = x[i]-y[i];
    diff[i] = d;
    if (d != 0.0) {
      vrec[nv].val = d;
      nv ++;
    }
  }

  if (nv == 0) {
    out_err(MWA, ERR_FILE, ERR_LINE,
    	_("All value pairs are equal. WILCOXON-Test thus not possible/has no meaning."));
    return;
  }

  qsort(vrec, nv, sizeof(SORTREC), wilcoxon_rank_compar);

  i = 0;
  k = 0;
  m = 0.;
  while (i<nv) {
    vrec[i].rank = (REAL)i+1.;
    if ( (i<(nv-1)) && (equal(fabs(vrec[i].val), fabs(vrec[i+1].val))) ) {
      k++;
      m += (REAL)i;
    }
    else {
      if (k!=0) {
        m += (REAL)i;
        k ++;
        m = m / (REAL)(k);
        for (j=i; j>(i-k); j--) {
          vrec[j].rank = m+1.;
        }
      }
      k=0;
      m=0.;
    }
    i++;
  }

  s_plus = 0.0;
  s_minus = 0.0;
  for (i=0; i<nv; i++) {
/*    printf("diff=%20.17f rang=%f\n", vrec[i].val, vrec[i].rank);     */
    if (vrec[i].val > 0.0) {
      s_plus += vrec[i].rank;
    }
    else {
      s_minus += vrec[i].rank;
    }
  }

  median = get_median(diff, n);
  navg = n*(n+1)/2;
  avg = (REAL*) m_calloc(navg, sizeof(REAL));
  nr = (REAL)nv;
  nd = (REAL)m;

  k = 0;
  for (i=0; i<n; i++) {
    for (j=i; j<n; j++) {
      avg[k] = (diff[i] + diff[j])/2.0;
      k ++;
    }
  }

  qsort(avg, navg, sizeof(REAL), real_compar_up);
/*  nr = 50; */
  if (n<26) {
    t = sig[n-6][2];
  }
  else {
    z = get_z(0.99);
    t = (int) (nd*(nd+1.0)/4.0 - z*sqrt(nd*(nr+1.0)*(2.0*nd+1.0)/24.0) - 0.5);
  }
  if ( (t<0) || (t>=navg) ) {
    t = 0;
  }

  conf_l = avg[t];
  conf_u = avg[navg-t-1];

/*
  for (i=0; i<navg; i++) {
    printf("%i %f\n", i, avg[i]);
  }

  printf("t= %i, z=%f\n", t, z);

  mean = get_mean(diff, nv);



  if (nv % 2 == 1) {
    i = (nv-1)/2;
    median = vrec[i].val;
  }
  else {
    i = (nv/2) -1;
    median = vrec[i].val;
    i++;
    median = (median + vrec[i].val)/2.;
  }
*/

  s_min = (s_plus<s_minus) ? s_plus : s_minus;
  out_start();
  colorize(ClHeader);
  out_r(_("\nResult Wilcoxon-Test:\n") );
  colorize(ClDefault);
  out_r(_("Positive rank sum S+  : %g\n"), s_plus);
  out_r(_("Negative rank sum S-  : %g\n"), s_minus);
  out_r(_("Smallest rank sum S   : %g\n"), s_min);
  out_r(_("Number of value pairs : %i\n"), n);
  out_r(_("Size of the set       : %i\n"), nv);
  out_r(_("Number 0-differences  : %i\n"), n-nv);
  out_r(_("Mean of differences   : %g +/- %g\n"), get_mean(diff, n),
	get_sdv(diff, n));
  out_r(_("Median of differences : %f [%f, %f] (99%%)\n\n"),
	median, conf_l, conf_u);
/*  out_r("Mittel der Differenzen: %f\n\n", mean);  */

  out_r(_("Hypothesis H0: x and y are 'treated' equally versus\n") );
  out_r(_("Hypothesis H1: x and y are 'treated' unequally (-> two-sided)\n\n"));

  if (nv<6) {
    out_r(_("Calculation of probability not possible if n < 6!\n") );
    out_end();
    return;
  }
  else if (nv<26) {
    i = nv-6;
    out_r(_("Critical values for S (two-sided) from table:\n") );
    out_r("   5%%   2%%   1%%\n");
    out_r(" %4hi %4hi %4hi\n", sig[i][0], sig[i][1], sig[i][2]);

    if (s_min <= (REAL)sig[i][2]) {
      alpha = 1;
    }
    else if (s_min <= (REAL)sig[i][1]) {
      alpha = 2;
    }
    else if (s_min <= (REAL)sig[i][0]) {
      alpha = 5;
    }
    else {
      alpha = -1;
    }

    if (alpha != -1) {
      out_r(_("H0 rejected at level of significance of %i%% (two-sided)\n\n"),
	  alpha);
    }
    else {
      out_r(_("H0 can not be rejected\n\n") );
    }
  }
/*  else {  */

  z = (s_min -  (nr*(nr+1.0))/4.0) /
      sqrt( (nr*(nr+1.0)*(2.0*nr+1.0))/24.0 );

  out_r(_("normally distributed variable z = %f\n"), z);
  prob = get_norm_int(z);
  out_r(_("Probability of H0 (one-sided) = %6.4f\n"), prob);
  out_r(_("Probability of H0 (two-sided) = %6.4f\n\n"), prob*2.0);
/*  }      */
  out_end();
}

/* ====================================================================== */


void outlier(int xx_ind, int n) {
  PREAL x = xx[xx_ind];
  FILE *outfile;
  char analias[80];
  int i, iy, iout, i_l, i_u;
  REAL index, q_l, q_u, w_l, w_u, m_l, m_u, max, min;
  REAL median, mean, no_l, no_u;
  REAL *xs;

  if (get_min(x,n)==get_max(x,n)) {
    out_err(MAT, ERR_FILE, ERR_LINE, _("All values are equal!"));
    return;
  }
  xs = (REAL*) m_calloc(n, sizeof(REAL));
  for (i=0; i<n; i++) {
    xs[i] = x[i];
  }
  qsort(xs, n, sizeof(REAL), real_compar_up);

  if (n % 2 == 1) {
    i = (n-1)/2;
    median = xs[i];
  }
  else {
    i = (n/2) -1;
    median = xs[i];
    i++;
    median = (median + xs[i])/2.;
  }

  mean = get_mean(xs, n);
  max = xs[n-1];
  min = xs[0];
  index = (REAL)n*0.25;
  if (index==floor(index)) {
    i_l = (int)index-1;
    i_u = (int)index;
  }
  else {
    i_l = (int) floor(index-1.);
    i_u = (int) ceil(index-1.);
  }

  q_l = 0.5 * (xs[i_l]+xs[i_u]);
  q_u = 0.5 * (xs[n-i_l-1]+xs[n-i_u-1]);

  m_l = q_l - 1.5*(q_u-q_l);
  m_u = q_u + 1.5*(q_u-q_l);

  w_l = max;
  w_u = min;

  for (i=0; i<n; i++) {
    if ((xs[i]>w_u) && (xs[i]<=m_u)) {
      w_u = xs[i];
    }
    if ((xs[i]<w_l) && (xs[i]>=m_l)) {
      w_l = xs[i];
    }
  }

  no_u = median + 1.58*(q_u-q_l)/sqrt((REAL)n);
  no_l = median - 1.58*(q_u-q_l)/sqrt((REAL)n);

  if (!noplot) {
    plot_box(x, n, median, mean, q_l, q_u, w_l, w_u, no_l, no_u, get_label(x));
  }

  out_start();
  iout = 0;
  for (i=0; i<n; i++) {
    if ((x[i]>w_u) || (x[i]<w_l)) {
      iout ++;
      out_r(_("Outliers:  x[%i]=%f\n"), i+1, x[i]);
    }
  }
  out_r(_("\n%i possible outliers found\n\n"), iout);

  if (iout == 0) {
    out_end();
    return;
  }
 #ifndef STATIST_X
  out_i(_("Delete outliers? (%s) "), _("y/N") );
  GETNLINE;
  if(!(empty) && (line[0]!=_("y")[0] || line[0] != _("Y")[0])){
      out_end();
      return;
  }

  iy = 0;
  strcpy(analias, "out_");
  strncat(analias, alias[acol[0]], 79-strlen(analias));
  make_new_col(analias, nn[xx_ind]);
  outfile = tmpptr[ncol - 1];
  alloc_cols(1, 3);
  x = xx[xx_ind]; /* We redeclare `x' because alloc_cols delete all columns from
		     memory */
  for (i=0; i<nn[xx_ind]; i++) {
    if (!((x[i]>w_u) || (x[i]<w_l))) {
     iy  ++;
     FWRITE(&x[i], sizeof(REAL), 1, outfile);
    }
    else {
     iy  ++;
     FWRITE(&SYSMIS, sizeof(REAL), 1, outfile);
    }
  }

  out_r(_("%i possible outliers deleted\n"), iout);
				/* ncol is the number, but counting starts at 0
				so we have to substract one the the right one*/
  out_r(_("Created new column %s having %i values!\n\n"), alias[ncol-1], iy);
#endif // ndef STATIST_X
  out_end();
}

/* =================================================================== */

char *center(char result[80], char *str, int width){
  int i, j, l;
  width += strlen(str) - stringLen(str);
  for(i = 0; i < 80; i++)
    result[i] = 0;
  l = strlen(str);
  i = (width - l) / 2;
  for(j = 0; j < i; j++)
    result[j] = ' ';
  for(i = 0; i < l; i++){
    result[j] = str[i];
    j++;
  }
  while(j < width){
    result[j] = ' ';
    j++;
  }
  return result;
}

void freq_table(){
  typedef struct { 
    char *c[6]; /* the value, the frequency, and the percentages */
    void *next;
  } FreqRow;

  FreqRow *row = NULL, *first_row = NULL;
  char b[20];
  char h[256], l[256];
  PREAL y = xx[acol[0]];
  REAL x, p, vp, ap = 0.0, vap = 0.0; /* value, global %, valid %, 
					 cumulated %, valid cumulated % */
  int i, n, m = 0, w[6], r; /* r = 3 because of the table header */
  BOOLEAN truncated = FALSE;
  qsort(y, nn[acol[0]], sizeof(REAL), real_compar_up);

  /* Setting default column widths */
  w[0] = stringLen(_("Value")) + 4;
  w[1] = stringLen(_("N")) + 2;
  w[2] = 8;
  w[3] = stringLen(_("Valid %")) + 2;
  w[4] = stringLen(_("Cum. %")) + 2;
  w[5] = stringLen(_("V. Cum. %")) + 2;
  if(w[3] < 8)
    w[3] = 8;
  if(w[4] < 8)
    w[4] = 8;
  if(w[5] < 8)
    w[5] = 8;

  /* Counting frequencies */
  i = 0;
  do{
    x = y[i];
    n = 0;
    while(i < nn[acol[0]] && y[i] == x){
      n++;
      i++;
    }
    if(row == NULL){
      first_row = (FreqRow*)m_calloc(1, sizeof(FreqRow));
      row = first_row;
    } else{
      row->next = (FreqRow*)m_calloc(1, sizeof(FreqRow));
      row = row->next;
    }
    for(r = 2; r < 6; r++)
      row->c[r] = (char*)m_calloc(7, sizeof(char));
    sprintf(b, "%i", n);
    row->c[1] = (char*)m_calloc((strlen(b) + 1), sizeof(char));
    strcpy(row->c[1], b);
    p = 100 * (REAL)n / (REAL)nn[acol[0]];
    sprintf(row->c[2], "%6.2f", p);
    ap += p;
    sprintf(row->c[4], "%6.2f", ap);
    if(x == SYSMIS){
      m = n;
      row->c[0] = (char*)m_calloc((strlen(_("Missing")) + 1), sizeof(char));
      strcpy(row->c[0], _("Missing"));
      strcpy(row->c[3], "--- ");
      strcpy(row->c[5], "--- ");
    } else{
      if(names[acol[0]]){
	for(r = 0; r <  names[acol[0]]->n; r++){
	  if(x == names[acol[0]]->v[r]){
	    row->c[0] = names[acol[0]]->l[r];
	    break;
	  }
	}
	if(row->c[0] == NULL){
	  sprintf(b, "%6g", x);
	  row->c[0] = (char*)m_calloc((strlen(b) + 1), sizeof(char));
	  strcpy(row->c[0], b);
	}
      } else{
	  sprintf(b, "%6g", x);
	  row->c[0] = (char*)m_calloc((strlen(b) + 1), sizeof(char));
	  strcpy(row->c[0], b);
      }
      vp = 100 * (REAL)n / (REAL)(nn[acol[0]] - m);
      sprintf(row->c[3], "%6.2f", vp);
      vap += vp;
      sprintf(row->c[5], "%6.2f", vap);
    }
  } while(i < nn[acol[0]]);
  if(first_row == NULL){
    out_err(ERR, ERR_FILE, ERR_LINE, "first_row == NULL!");
    return;
  }

  /* Calculating column widths for "value" and "N"*/
  row = first_row;
  while(row){
    for(i = 0; i < 2; i++)
      if((stringLen(row->c[i]) + 2) > w[i])
	w[i] = stringLen(row->c[i]) + 2;
    row = row->next;
  }

  /* Printing header */
  memset(l, 0, 256);
  if(names[acol[0]] && names[acol[0]]->ctitle)
    strncpy(l, names[acol[0]]->ctitle, 255);
  else
    strncpy(l, alias[acol[0]], 255);
  colorize(ClHeader);
  out_r(_("Frequencies: %s\n"), l);
  colorize(ClDefault);
  strcpy(b, _("Value"));
  sprintf(h, " %s ", center(l, b, w[0]));
  strcpy(b, _("N"));
  strcat(h, center(l, b, w[1]));
  strcat(h, center(l, "%", w[2]));
  strcpy(b, _("Valid %"));
  strcat(h, center(l, b, w[3]));
  strcpy(b, _("Cum. %"));
  strcat(h, center(l, b, w[4]));
  strcpy(b, _("V. Cum. %"));
  strcat(h, center(l, b, w[5]));
  r = stringLen(h) + 1;
  for(i = 0; i < r; i++)
    l[i] = '=';
  l[i] = 0;
  out_r("%s\n", l);
  colorize(ClHeader);
  out_r("%s\n", h);
  colorize(ClDefault);
  for(i = 0; i < r; i++)
    l[i] = '-';
  l[i] = 0;
  out_r("%s\n", l);
  for(i = 0; i < r; i++)
    l[i] = '=';
  l[i] = 0;

  /* Printing rows */
  row = first_row;
  r = 4; /* table title + table header = 4 lines */
  int diff;
  while(row){
    diff = strlen(row->c[0]) - stringLen(row->c[0]);
    if(names[acol[0]])
      sprintf(b, " %%-%is", w[0]+diff);
    else
      sprintf(b, " %%%is", w[0]+diff);
    out_r(b, row->c[0]);
    for(i = 1; i < 6; i++){
      diff = strlen(row->c[i]) - stringLen(row->c[i]);
      sprintf(b, "%%%is", w[i]+diff);
      out_r(b, row->c[i]);
    }
    out_r("\n");
    if ((r+1) % (SCRLINES - 1) == 0){
      mywait();
      if(!empty){
	truncated = TRUE;
	break;
      }
    }
    r++;
    row = row->next;
  }
  
  if(truncated){
    out_r("\n");
    out_err(MWA, ERR_FILE, ERR_LINE,
	_("         ---  TABLE TRUNCATED  ---"));
  } else{
    out_r("%s\n", l);
  }
}

typedef struct { 
  REAL v; /* the same as xx[0][j] */
  int i;  /* the index: j */
} Mirror;

int mirror_compare_up(const void *a, const void *b) {
  if (((Mirror*)a)->v > ((Mirror*)b)->v){
    return  1;
  } else{
    if (((Mirror*)a)->v < ((Mirror*)b)->v){
      return -1;
    } else{
      return  0;
    }
  }
}

int whatwidth_i(int i){
  char s[20];
  sprintf(s, "%i", i);
  return(strlen(s));
}

int whatwidth_r(REAL r){
  char s[20];
  int i, l;
  if(r == SYSMIS)
    return 1;
  i = floor(r);
  sprintf(s, "%i", i);
  l = strlen(s);
  if(i >= 100)
    return(l + 3);		/*  100.00 */
  else
    if(i <= -100)
      return(l + 4);		/* -100.00 */
    else
      if(i >=0 && i < 100)	/*  1.0000 */
	return(l + 5);
      else			/* -1.0000 */
	return(l + 6);
}

/* The calculus and the printing of the Table of Means are truncated if there
 * are more than MRESULT. I think that such a table is only meaniful when there
 * is a small number of classes of frequency, and I suppose that the user has
 * chosen this analysis by mistake in these cases */
void compare_means(int c){
  typedef struct { 
    PREAL r; /* the class of frequency, and the means */
    PREAL s; /* the standard deviation */
    int *n;  /* the frequency, and n. of valid values used to calc. the mean */
    void *next;
  } MeanRow;

  int i, j, k, p, n, nc, w;
  int *nw, *rw, *sw; /* lists of widths for columns */ 
  char b[80], *l, *L; 
  MeanRow *row, *firstrow; 
  BOOLEAN truncated = FALSE, vlabelfound = FALSE;
  Mirror *y = (Mirror*)m_calloc(nn[acol[0]], sizeof(Mirror)); 
  REAL *z = (PREAL)m_calloc(nn[acol[0]], sizeof(REAL));

  for(i = 0; i < nn[acol[0]]; i++){
    y[i].v = xx[acol[0]][i];
    y[i].i = i;
  }

  /* Sort "y" to calculate frequencies, as in freq_table(), but keep track of
   * original indexes to calculate means. */
  qsort(y, nn[acol[0]], sizeof(Mirror), mirror_compare_up);

  /* Calculating means, and filling the table of means */
  j = 0;
  k = 0;
  nc = 0;
  row = (MeanRow*)m_calloc(1, sizeof(MeanRow));
  firstrow = row;
  do{
    row->r = (PREAL)m_calloc(c, sizeof(REAL));
    row->s = (PREAL)m_calloc(c, sizeof(REAL));
    row->n = (int*)m_calloc(c, sizeof(int));
    row->next = NULL;
    row->r[0] = y[j].v;
    n = 0;
    i = j;
    while(j < nn[acol[0]] && y[j].v == row->r[0]){ /* Counting frequency */
      n++;
      j++;
    }
    row->n[0] = n;
    for(p = 1; p < c; p++){ /* Filling "z[]" to calculate the mean */
      n = 0;
      for(k = i; k < j; k++){
	if(xx[acol[p]][y[k].i] != SYSMIS){
	  z[n] = xx[acol[p]][y[k].i];
	  n++;
	}
      }
      if(n > 0)
	row->r[p] = get_mean(z, n);
      else
	row->r[p] = 10.0;           /* This value will not be printed,   */
      if(n > 1)                     /* but will be used to calculate the */
	row->s[p] = get_sdv(z, n);  /* column width: 0.0000 is broader   */
      else			    /* than 10.00 */
	row->s[p] = 10.0;
      row->n[p] = n;
    }
    nc++;
    if(nc > MRESULT)
      break;
    if(j < nn[acol[0]]){
      row->next = (MeanRow*)m_calloc(1, sizeof(MeanRow));
      row = row->next;
    }
  } while(j < nn[acol[0]]);

  /* Calculating the necessary widths for columns */
  rw = (int*)m_calloc(c, sizeof(int));
  sw = (int*)m_calloc(c, sizeof(int));
  nw = (int*)m_calloc(c, sizeof(int));
  for(p = 0; p < c; p++){
    if(p == 0){
      rw[p] = stringLen(_("value"));
      sw[p] = 0;
    } else{
      rw[p] = stringLen(_("mean"));
      sw[p] = stringLen(_("sdev"));
    }
    nw[p] = 2;
    row = firstrow;
    do{
      i = whatwidth_i(row->n[p]);
      if(i > nw[p])
        nw[p] = i;
      i = whatwidth_r(row->r[p]);
      if(i > rw[p])
	rw[p] = i;
      if(p > 0){
	i = whatwidth_r(row->s[p]);
	if(i > sw[p])
	  sw[p] = i;
      }
      if(row->next)
	row = row->next;
    }while(row->next);
  }
  if(names[acol[0]]){
    k = rw[0];
    for(i = 0; i < names[acol[0]]->n; i++)
      if((2 + stringLen(names[acol[0]]->l[i])) > k)
	k = stringLen(names[acol[0]]->l[i]) + 2;
    rw[0] = k;
  }

  /* calculating the total width of the table */
  w = 0;
  sw[0] = 0;
  for(p = 0; p < c; p++){
    if(p != 0){
      rw[p] += 2;		/* width += blank space between columns */
      sw[p] += 2;
    }
    nw[p] += 2;
    w += rw[p] + sw[p] + nw[p];
  }

  /* Printing the table of means */
  if(names[acol[0]] && names[acol[0]]->ctitle)
    strncpy(b, names[acol[0]]->ctitle, 79);
  else
    strncpy(b, alias[acol[0]], 79);
  colorize(ClHeader);
  out_r(_("Comparison of means: %s\n"), b);
  colorize(ClDefault);
  l = (char*)m_calloc((w+3), sizeof(char));
  L = (char*)m_calloc((w+3), sizeof(char));
  for(i = 0; i <= w; i++){
    l[i] = '-';
    L[i] = '=';
  }
  l[i] = '\n';
  L[i] = '\n';
  out_r(L);
  colorize(ClHeader);
  for(p = 0; p < c; p++){
    i = rw[p] + sw[p] + nw[p];
    out_r("%s", center(b, alias[acol[p]], i));
  }
  colorize(ClDefault);
  out_r("\n");
  for(p = 0; p < c; p++){
    strcpy(b, " ");
    for(i = 0; i < (rw[p] + sw[p] + nw[p] - 1); i++)
      strcat(b, "-");
    out_r(b);
  }
  out_r("\n");
  colorize(ClHeader);
  out_r("%s ", center(b, _("value"), rw[0]));
  out_r("%s", center(b, _("N"), nw[0]));
  for(p = 1; p < c; p++){
    out_r("%s", center(b, _("mean"), rw[p]));
    out_r("%s", center(b, _("sdev"), rw[p]));
    out_r("%s", center(b, _("N"), nw[p]));
  }
  colorize(ClDefault);
  out_r("\n%s", l);

  i = 5; /* header = 5 lines */
  row = firstrow;
  int diff;
  while(row){
    for(p = 0; p <  c; p++){
      if(p == 0){
	if(names[acol[0]] && names[acol[0]]->n > 0){  /* print value label */
	  if(row->r[p] == SYSMIS){
	    sprintf(b, " %%-%is%%%ii", rw[p], nw[p]);
	    out_r(b, NODATA, row->n[p]);
	  } else{
	    for(k = 0; k <  names[acol[0]]->n; k++){
	      if(row->r[p] == names[acol[0]]->v[k]){
		diff = strlen(names[acol[0]]->l[k]) - stringLen(names[acol[0]]->l[k]);
		sprintf(b, " %%-%is%%%ii", rw[p]+diff, nw[p]);
		out_r(b, names[acol[0]]->l[k], row->n[p]);
		vlabelfound = TRUE;
		break;
	      }
	    }
	    if(vlabelfound)
	      vlabelfound = FALSE;
	    else{
	      sprintf(b, " %%-%i.%if%%%ii", rw[p], 4, nw[p]);
	      out_r(b, row->r[p], row->n[p]);
	    }
	  }
	} else{ /* no label for values: print value */
	  if(row->r[p] == SYSMIS){
	    sprintf(b, "%%%is%%%ii", rw[p], nw[p]);
	    out_r(b, NODATA, row->n[p]);
	  } else{
	    sprintf(b, "%%%i.%if%%%ii", rw[p], 4, nw[p]);
	    out_r(b, row->r[p], row->n[p]);
	  }
	}
      } else{
	if(row->n[p] > 1){
	  if((row->r[p] >= 100.0 || row->r[p] <= -100.0) &&
	      (row->s[p] >= 100.0 || row->s[p] <= -100.0))
	    sprintf(b, "%%%i.%if%%%i.%if%%%ii", rw[p], 2, sw[p], 2, nw[p]);
	  else if((row->r[p] >= 100.0 || row->r[p] <= -100.0) &&
	      (row->s[p] < 100.0 && row->s[p] > -100.0))
	    sprintf(b, "%%%i.%if%%%i.%if%%%ii", rw[p], 2, sw[p], 4, nw[p]);
	  else if((row->r[p] < 100.0 && row->r[p] > -100.0) &&
	      (row->s[p] >= 100.0 || row->s[p] <= -100.0))
	    sprintf(b, "%%%i.%if%%%i.%if%%%ii", rw[p], 4, sw[p], 2, nw[p]);
	  else
	    sprintf(b, "%%%i.%if%%%i.%if%%%ii", rw[p], 4, sw[p], 4, nw[p]);
	  out_r(b, row->r[p], row->s[p], row->n[p]);
	} else{
	  if(row->n[p] == 1){
	    if(row->r[p] >= 100.0 || row->r[p] <= -100.0) 
	      sprintf(b, "%%%i.%if%%%is%%%ii", rw[p], 2, sw[p], nw[p]);
	    else
	      sprintf(b, "%%%i.%if%%%is%%%ii", rw[p], 4, sw[p], nw[p]);
	    out_r(b, row->r[p], "--", row->n[p]);
	  } else{
	    sprintf(b, "%%%is%%%is%%%ii", rw[p], sw[p], nw[p]);
	    out_r(b, "--", "--", 0);
	  }
	}
      }
    }
    out_r("\n");
    i++;
    if ((i+1) % (SCRLINES - 1) == 0){
      mywait();
      if(!empty){
	truncated = TRUE;
	break;
      }
    }
    row = row->next;
  }
  if(truncated){
    out_r("\n");
    out_err(MWA, ERR_FILE, ERR_LINE,
	_("         ---  TABLE TRUNCATED  ---"));
  } else
    if(nc > MRESULT){
      out_r("\n");
      out_err(MWA, ERR_FILE, ERR_LINE,
	_("         ---  TABLE TRUNCATED  ---"));
    } else
      out_r("%s\n", L);
}

/* =================================================================== */


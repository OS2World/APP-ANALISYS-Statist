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

/* funcs.h fuer STATIST */

typedef struct {
  short unsigned int *a, n;  
} TUPEL;


extern REAL  get_sum(REAL x[], int n);
extern REAL  get_qsum(REAL x[], int n);
extern REAL  get_xysum(REAL x[], REAL y[], int n);
extern REAL  get_sdv(REAL x[], int n);
extern REAL  get_mean(REAL x[], int n);
extern REAL  get_max(REAL x[], int n);
extern int   get_maxint(int x[], int n);
extern REAL  get_min(REAL x[], int n);
extern REAL  get_median(REAL x[], int n);
extern REAL  get_rank_correlation(REAL x[], REAL y[], int n);
extern int   get_round(REAL x);
extern REAL  get_sgn(REAL x);
extern REAL  get_norm_int(REAL x);
extern REAL  get_norm_ord(REAL x);
extern REAL  get_t_int(REAL t, int f);
extern REAL  get_chi_int(REAL chi, int f);
extern REAL  get_f_int(REAL f, int f1, int f2);
extern REAL  rise(REAL x, int exp);
extern REAL  get_z(REAL alpha);
extern REAL  get_t(REAL alpha, int df);
extern REAL  get_ln_0(REAL x);

extern int  rank_compar(const void *x, const void *y);
extern int  real_compar_down(const void *x, const void *y);
extern int  real_compar_up(const void *x, const void *y);
extern int  u_rank_compar(const void *x, const void *y);
extern int  wilcoxon_rank_compar(const void *x, const void *y);
extern BOOLEAN equal(REAL x, REAL y);
extern int pks(REAL d, int n);
extern REAL get_multiple_reg(REAL y[], PREAL x[], int nrow, int ncol,
			      REAL b[], REAL *sdv, REAL *f_calc);
extern void copy_tupel(TUPEL* dest, const TUPEL* src);
extern void print_tupel(TUPEL atupel);
extern TUPEL* create_tupel(TUPEL *atupel, int ndata);
extern BOOLEAN equal_tupel(TUPEL t1, TUPEL t2);
extern REAL get_cross_validate(REAL y[], PREAL x[], int nrow, int ncol, 
			REAL ypred[]);



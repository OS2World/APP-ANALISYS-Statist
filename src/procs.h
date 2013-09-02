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

/* procs.h for STATIST */



void histogram(REAL x[], int n, int mclass, REAL min, REAL max);
void standard(REAL x[], int n, int mclass, REAL min, REAL max);
extern void lin_reg(REAL x[], REAL y[], int n);

extern void rank_order(REAL x[], REAL y[], int n);
extern void t_test(REAL x[], int nx, REAL y[], int ny);
extern void pair_t_test(REAL x[], REAL y[], int n);
extern void kolmo_test(REAL x[], int n);
extern void percentiles(REAL x[], int n);
extern void probit(REAL dose[], REAL num[], REAL effect[], int n);
extern  void correl_matrix(PREAL xx[], int nrow, int ncol);
extern void rank_matrix(PREAL xx[], int nrow, int ncol);

extern void multiple_reg(REAL y[], PREAL x[], int nrow, int ncol);
extern void poly_reg(REAL x[], REAL y[], int n, int m);
extern void part_corr(PREAL xx[], int nrow, int ncol);
extern void point_biserial_reg(REAL x[], REAL y[], int n);
extern void vierfeld_test(REAL x[], REAL y[], int n);
extern void tafel_test(PREAL xx[], int nrow, int ncol);
extern void u_test(REAL x[], int nx, REAL y[], int ny);
extern void kruskal_test(PREAL x[], int nrow[], int ncol);
extern void ks(REAL data1[], int n1, REAL data2[], int n2, REAL *d, REAL *prob);
extern void equal_freq(REAL x[], int n);
extern void compare_freq(REAL x[], int nx, REAL y[], int ny);
extern void wilcoxon_test(REAL x[], REAL y[], int n);
extern void outlier(int xx_ind, int n);

extern void random_tupel(REAL y[], PREAL x[], int nrow, int ncol, int ntup);
extern void cross_validate(REAL y[], PREAL x[], int nrow, int ncol);
extern void freq_table();
extern void compare_means(int n);

/*
  extern void compare_dist(REAL x[], int nx, REAL y[], int ny);
*/




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

/* plot.h fuer STATIST */

extern BOOLEAN init_gnuplot();
extern BOOLEAN plot_pair(REAL x[], REAL y[], int n, REAL a0, REAL a1,
		  char *xlab, char *ylab);
extern BOOLEAN plot_tripl(REAL x[], REAL y[], REAL z[], int n,
                   REAL a0, REAL a1, REAL a2,
		   char *xlab, char *ylab, char *zlab);
extern BOOLEAN plot_histo(REAL x[], int y[], int n, REAL step, REAL data[],
			  char *datalab);
extern BOOLEAN plot_histo_discrete(REAL x[], int y[], int n, REAL step, REAL data[],
			  char *datalab);
extern BOOLEAN plot_probit(REAL dose[], REAL num[], REAL effect[], int n,
                           REAL a0, REAL a1, REAL dose0, REAL dose1,
			   char *doselab, char *effectlab);
extern BOOLEAN plot_poly(REAL x[], REAL y[], int n, REAL a[], int npoly,
			 char *xlab, char *ylab);
extern BOOLEAN plot_box(REAL x[], int n, REAL median, REAL mean, REAL q_l,
                        REAL q_u, REAL w_l, REAL w_u, REAL no_l,
			REAL no_u, char *xlab);
BOOLEAN plot_cdf_ks(REAL z[], int n, REAL z0, REAL fn1, REAL fn2, REAL x[],
		 char *xlab);
BOOLEAN plot_cdf(REAL z[], int n, char *xlab);

#ifndef MSDOS
extern void save_png();
extern BOOLEAN plot_command();
#endif

extern char GPL_DAT[];
extern int histsize; /* number of last gnuplot commands to memorize */
extern char **plothist;
extern char *gnuplot_charset;
extern BOOLEAN has_graph;

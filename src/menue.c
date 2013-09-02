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
**  some changes by Bernhard Reiter  http://www.usf.Uni-Osnabrueck.DE/~breiter
**  $Id: menue.c,v 1.28 2009/12/10 13:17:03 jakson Exp $
***************************************************************/

/* menue.c for STATIST */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "statist.h"
#include "procs.h"
#include "funcs.h"
#include "data.h"
#include "plot.h"
#include "menue.h"

#include "gettext.h"

int getint() {
   int anint=-1;
   if (!empty) {
      status = sscanf(line, "%i", &anint);
      if (status != 1) {
         out_err(ERR, ERR_FILE, ERR_LINE, _("No valid number!") );
         empty = TRUE;
      }
   }
   return anint;
}

REAL getreal() {
   REAL areal= -1.0;
   if (!empty) {
      status = sscanf(line, "%lf", &areal);
      if (status != 1) {
         out_err(ERR, ERR_FILE, ERR_LINE, _("Invalid number!") );
         empty = TRUE;
      }
   }
   return areal;
}


BOOLEAN equal_rows(int nm) {
/* Tests whether columns have different lengths */
  int i;
  for (i=1; i<nm; i++) {
    if (nn[acol[i]] != nn[acol[0]]) {
      out_err(ERR, ERR_FILE, ERR_LINE,
	  _("Columns have different number of entries!") );
      return FALSE;
    }
  }
  return TRUE;
}

void printline(){
  colorize(ClMenuSep);
  out_d("================================================\n\n");
  colorize(ClDefault);
}

/* ====================================================================== */

void main_menue() {
   int choice = 99;

   while (choice != QUIT) {
     out_d( _("MAIN MENU: \n\n") );
     out_d( _("   0 = Quit\n") );
     out_d( _("   1 = Data management\n") );
     out_d( _("   2 = Regressions and correlations\n") );
     out_d( _("   3 = Tests\n"));
     out_d( _("   4 = Miscellaneous\n") );
     out_d( _("   5 = Data manipulation\n") );
     out_d( _("   6 = Preferences\n") );

     out_d( _("\n  Your choice: ") );
     GETNLINE;
     status = sscanf(line,"%i", &choice);
     if ((status==0) || (empty)) {
       choice = 99;
     }
     out_d("\n\n");

     switch(choice) {
       case 0 : ;
	       break;
       case 1: {
		 printline();
		 data_menu();
	       }
	       break;
       case 2: {
		 printline();
		 regress_menu();
	       }
	       break;
       case 3: {
		 printline();
		 test_menu();
	       }
	       break;
       case 4: {
		 printline();
		 misc_menu();
	       }
	       break;
       case 5: {
		 printline();
		 manipulate_menu();
	       }
	       break;
       case 6: {
		 printline();
		 prefs_menu();
	       }
	       break;
       default:
	       out_err(ERR, ERR_FILE, ERR_LINE,
		   _("Illegal instruction!") );
	       break;
     }
     printline();
   }
}


/* ====================================================================== */
void misc_menu() {
   int choice = 99, i;
   REAL min, max;

   while (choice != QUIT) {
     out_d(_("MISCELLANEOUS: \n\n") );
     out_d(_("   0 = Main menu\n") );
     out_d(_("   1 = Standard deviation, mean, median, etc. \n") );
     out_d(_("   2 = Probit analysis\n") );
     out_d(_("   3 = Outliers & Box-Whisker-plot\n") );
     out_d(_("   4 = Percentiles\n") );
     out_d(_("   5 = Frequency table\n") );
     out_d(_("   6 = Compare means\n") );
#ifndef MSDOS
     if(gnupl_open == TRUE && has_graph == TRUE)
       out_d(_("   8 = Save last gnuplot graphic as png\n") );
     if (!noplot)
       out_d(_("   9 = Enter gnuplot commands\n") );
#endif

     out_d(_("\n  Your choice: ") );
     GETRLINE;
     status = sscanf(line,"%i", &choice);
     if ((status==0) || (empty)) {
       choice = 99;
     }
     out_d("\n\n");

     switch(choice) {
       case 0: ;
               break;
       case 1:
               i = getcols(1, 1, TRUE);
	       if(i == 0)
		 return;
	       out_i(_("Number of classes (0 or RETURN for auto): ") );
	       GETNLINE;
	       i = getint();
	       if (empty) {
		 i = 0;
	       }
	       if (i>MCLASS) {
		 out_err(ERR, ERR_FILE, ERR_LINE,
		     _("More than %i classes not allowed!"), MCLASS);
		 break;
	       }
	       else if ((i<0)||(i==1)) {
		 out_err(ERR, ERR_FILE, ERR_LINE,
		     _("Please enter a meaningful number of classes!") );
		 break;
	       }
	       else {
		 out_i("Minimum (%f): ", get_min(xx[acol[0]], vn[acol[0]]));
		 GETNLINE;
		 min = getreal();
		 if (empty) {
		   min = get_min(xx[acol[0]], vn[acol[0]]);
		 }
		 out_i("Maximum (%f): ", get_max(xx[acol[0]], vn[acol[0]]));
		 GETNLINE;
		 max = getreal();
		 if (empty) {
		   max = get_max(xx[acol[0]], vn[acol[0]]);
		 }
	       }
	       standard(xx[acol[0]], vn[acol[0]], i, min, max);
               break;
       case 2:
               out_i(_("Set number to 100 (--> percent? ) %s "), _("y/N") );
               GETNLINE;
               if (!(empty) && (line[0] == _("y")[0] || line[0] == _("Y")[0])){
                 out_i(_("Please select columns containing dose and "
		       "effect data\n") );
                 i = getcols(2, 2, TRUE);
                 if (i != 0) {
                   tempcol = (REAL*)m_calloc(vn[acol[0]], sizeof(REAL));
                   for (i=0; i<vn[acol[0]]; i++ ) {
                     tempcol[i] = 100.0;
                   }
                   probit(xx[acol[0]], tempcol, xx[acol[1]], vn[acol[0]]);
                 }
               }
               else {
                 out_i(_("Please select columns containing dose, number, "
		       "effect data\n") );
                 i = getcols(3, 3, TRUE);
                 if (i != 0) {
                   probit(xx[acol[0]], xx[acol[1]], xx[acol[2]], vn[acol[0]]);
                 }
               }
               break;
       case 3:
               i = getcols(1, 1, TRUE);
               if(i != 0)
		 outlier(acol[0], vn[acol[0]]);
	       break;
       case 4:
               i = getcols(1, 1, TRUE);
               if(i != 0)
                 percentiles(xx[acol[0]], vn[acol[0]]);
               break;
       case 5:
	       i = getcols(1, 1, 3);
	       if(i != 0)
		 freq_table();
	       break;
       case 6:
               out_i(_("Please select columns to compare means,\n") );
               out_i(_("(the first one will be taken as y-value)\n") );
	       i = getcols(2, ncol, 3);
	       if(i == 0)
		 break;
               if (equal_rows(i))
                 compare_means(i);
	       else
		 out_err(ERR, ERR_FILE, ERR_LINE, _("The columns must have "
		       "the same number of data points for this analysis!"));
               break;
#ifndef MSDOS
       case 8:
	       save_png();
	       break;
       case 9:
               if (!noplot) {
                 plot_command();
               }
               break;
#endif
       default:
	       out_err(ERR, ERR_FILE, ERR_LINE,	_("Illegal instruction!") );
               break;
     }
     m_freeall(); /* deallocate used memory for tmp variables */
     if (choice != 0) {
       mywait();
     }
   }
   return;
} /* misc_menu() */

/* =================================================================== */

void test_menu()  {
   int i, k, choice = 99;

   while (choice != QUIT) {
     out_d( _("TESTS:\n\n"));
     out_d(_("   0 = Main menu\n") );
     out_d(_("   1 = t-test for comparison of two means of two samples\n") );
     out_d(_("   2 = t-test for comparison of pairwise ascertained samples\n"));
     out_d(_("   3 = Test of normal distribution (KS-Lilliefors-Test)\n") );
     out_d(_("   4 = Chi^2-fourfold-test\n") );
     out_d(_("   5 = Chi^2 two-items-test\n") );
     out_d(_("   6 = u-test (Test of independence of two samples)\n") );
     out_d(_("   7 = H-test (Kruskal-Wallis) for k independent samples\n") );
     out_d(_("   8 = Wilcoxon-Rank-test for pairwise ascertained samples\n") );
     out_d(_("   9 = Chi^2-test of equal frequency\n") );
     out_d(_("  10 = Chi^2-test of correspondence between measured and "
	   "theoretical frequency\n") );

     out_d(_("\n  Your choice: ") );
     GETRLINE;
     status = sscanf(line,"%i", &choice);
     if ((status==0) || (empty)) {
       choice = 99;
     }
     out_d("\n\n");

     switch(choice) {
       case 0:
               break;
       case 1:
               i = getcols(2, 2, FALSE);
               if(i != 0)
                 t_test(xx[acol[0]],vn[acol[0]], xx[acol[1]],vn[acol[1]]);
               break;
       case 2:
               i = getcols(2, 2, TRUE);
               if(i != 0)
                 pair_t_test(xx[acol[0]], xx[acol[1]],vn[acol[0]]);
               break;
       case 3:
               i = getcols(1, 1, TRUE);
               if(i != 0)
                 kolmo_test(xx[acol[0]], vn[acol[0]]);
               break;
      case 4:
               i = getcols(2, 2, TRUE);
               if(i != 0)
                 vierfeld_test(xx[acol[0]], xx[acol[1]], vn[acol[0]]);
               break;
      case 5:
	       i = getcols(2, ncol, TRUE);
               if (i != 0){
                 yy = (REAL**)m_calloc(i, sizeof(PREAL));
                 for (k=0; k<i; k++) {
                   yy[k] = xx[acol[k]];
                 }
                 tafel_test(yy, vn[acol[0]], i);
               }
               break;
       case 6:
               i = getcols(2, 2, FALSE);
               if (i != 0)
                 u_test(xx[acol[0]], vn[acol[0]], xx[acol[1]], vn[acol[1]]);
               break;
       case 7:
	       i = getcols(3, ncol, FALSE);
               if (i != 0)  {
                 yy = (REAL**)m_calloc(i, sizeof(PREAL));
                 ny = (int*)m_calloc(i, sizeof(int));
                 for (k=0; k<i; k++) {
                   yy[k] = xx[acol[k]];
                   ny[k] = vn[acol[k]];
                 }
                 kruskal_test(yy, ny, i);
               }
               break;
       case 8:
	       i = getcols(2, 2, TRUE);
	       if (i != 0) {
		 wilcoxon_test(xx[acol[0]], xx[acol[1]], vn[acol[0]]);
	       }
	       break;
       case 9:
               i = getcols(1, 1, TRUE);
               if (i != 0) {
                 equal_freq(xx[acol[0]], vn[acol[0]]);
               }
               break;
       case 10:
	       out_i(_("Variable 1 (x) = measured distribution\n") );
	       out_i(_("Variable 2 (y) = expected (theoretical) distribution\n\n") );
               i = getcols(2, 2, FALSE);
               if (i != 0) {
                 compare_freq(xx[acol[0]], vn[acol[0]], xx[acol[1]], vn[acol[1]]);
               }
               break;
       default:
	       out_err(ERR, ERR_FILE, ERR_LINE,
		   _("Illegal instruction!") );
               break;
     }
     m_freeall(); /* deallocate used memory for tmp variables */
     if (choice != 0) {
       mywait();
     }
   }
   return;
} /* test_menu() */

/* =================================================================== */

void regress_menu()  {
   int i, k, choice = 99;

   while (choice != QUIT) {
     out_d(_("REGRESSION & CORRELATION: \n\n") );
     out_d(_("   0 = Main menu\n") );
     out_d(_("   1 = Linear regression and correlation\n") );
     out_d(_("   2 = SPEARMAN rank-correlation-coefficient\n") );
     out_d(_("   3 = Multiple linear correlation\n") );
     out_d(_("   4 = Partial linear correlation (maximum: 5 variables)\n") );
     out_d(_("   5 = Polynomial regression\n") );
     out_d(_("   6 = Matrix of the linear correlation coefficients\n") );
     out_d(_("   7 = Matrix of SPEARMAN correlation coefficients\n") );
     out_d(_("   8 = Point-biserial (linear) correlation\n") );
     out_d(_("   9 = Cross-validation of multiple linear regression\n") );
     out_d(_("  10 = Randomization of multiple linear regression\n") );
     out_d(_("\n  Your choice: ") );
     GETRLINE;
     status = sscanf(line,"%i", &choice);
     if ((status==0) || (empty)) {
       choice = 99;
     }
     out_d("\n\n");

     switch(choice) {
       case 0:
               break;
       case 1:
               i = getcols(2, 2, TRUE);
               if (i != 0) {
                 lin_reg(xx[acol[0]], xx[acol[1]], vn[acol[0]]);
               }
               break;
       case 2:
               i = getcols(2, 2, TRUE);
               if (i != 0) {
                 rank_order(xx[acol[0]], xx[acol[1]], vn[acol[0]]);
               }
               break;
       case 3:
               out_i(_("First selected column will be taken as y-value!") );
	       i = getcols(2, ncol, TRUE);
               if (i != 0) {
                 yy = (REAL**)m_calloc((i-1), sizeof(PREAL));
                 for (k=0; k<(i-1); k++) {
                   yy[k] = xx[acol[k+1]];
                 }
                 out_r(_("y = Column %s\n"),
		     			(alias[acol[0]]));
                 for (k=1; k<i; k++) {
                   out_r(_(" x[%i] = Column %s\n"),
			 		k, (alias[acol[k]]));
                 }
                 multiple_reg(xx[acol[0]], yy, vn[acol[0]], i-1);
               }
               break;
       case 4:
	       i = getcols(3, 5, TRUE);
	       if (i != 0) {
		 yy = (REAL**)m_calloc(i, sizeof(PREAL));
		 for (k=0; k<i; k++) {
		   yy[k] = xx[acol[k]];
		 }
		 part_corr(yy, vn[acol[0]], i);
	       }
	       break;
       case 5:
               out_i(_("Please select columns for x- and y-values "
		     "(variable 1 = x, variable 2 = y)\n\n") );
               i = getcols(2, 2, TRUE);
	       if (i != 0) {
		 out_i(_("Please enter order of regression polynom (max. %i): ")
		     , MPOLY);
                 GETBLINE;
                 i = getint();
                 if ( (i>0) && (i<=MPOLY) && (vn[acol[0]]==vn[acol[1]]) ) {
                   poly_reg(xx[acol[0]], xx[acol[1]], vn[acol[0]], i);
                 }
                 else {
                   out_err(ERR, ERR_FILE, ERR_LINE,
		       _("Illegal order or x- and y-columns have different "
			 "number of values!") );
                 }
               }
               break;
       case 6:
	       i = getcols(2, ncol, TRUE);
               if (i != 0) {
                 yy = (REAL**)m_calloc(i, sizeof(PREAL));
                 for (k=0; k<i; k++) {
                  yy[k] = xx[acol[k]];
                 }
		 correl_matrix(yy, vn[acol[0]], i);
               }
               break;
       case 7:
	       i = getcols(2, ncol, TRUE);
               if (i != 0) {
                 yy = (REAL**)m_calloc(i, sizeof(PREAL));
                 for (k=0; k<i; k++) {
                   yy[k] = xx[acol[k]];
                 }
		 rank_matrix(yy, vn[acol[0]], i);
               }
               break;
       case 8:
               out_i(_("First column must contain only 0's and 1's!"));
	       out_d("\n");
               i = getcols(2, 2, TRUE);
               if (i != 0) {
                 point_biserial_reg(xx[acol[0]], xx[acol[1]], vn[acol[0]]);
               }
               break;
       case 9:
               out_i( _("First selected column will be taken as y-value!"));
	       i = getcols(2, ncol, TRUE);
               if (i != 0) {
                 yy = (REAL**)m_calloc((i-1), sizeof(PREAL));
                 for (k=0; k<(i-1); k++) {
                   yy[k] = xx[acol[k+1]];
                 }
                 out_r( _("y = column %s\n"), (alias[acol[0]]));
                 for (k=1; k<i; k++) {
                   out_r( _("x[%i] = column %s\n"),
                       k, (alias[acol[k]]));
                 }
                 cross_validate(xx[acol[0]], yy, vn[acol[0]], i-1);
               }
               break;

       case 10:
               out_i(_("First selected column will be taken as y-value!"));
	       i = getcols(2, ncol, TRUE);
               if (i != 0) {
                 yy = (REAL**)m_calloc((i-1), sizeof(PREAL));
                 for (k=0; k<(i-1); k++) {
                   yy[k] = xx[acol[k+1]];
                 }
                 out_r(_("y = Column %s\n"), (alias[acol[0]]));
                 for (k=1; k<i; k++) {
                   out_r(_("x[%i] = Column %s\n"), k, (alias[acol[k]]));
                 }
		 out_i(_("Please enter number of randomizations: ") );
                 GETBLINE;
                 k = getint();
                 random_tupel(xx[acol[0]], yy, vn[acol[0]], i-1, k);
               }
               break;

       default:
	       out_err(ERR, ERR_FILE, ERR_LINE,	_("Illegal instruction!") );
               break;
     }
     m_freeall(); /* deallocate used memory for tmp variables */
     if (choice != 0) {
       mywait();
     }
   }
   return;
} /* regress_menu() */

/* =================================================================== */

void data_menu() {
   int i, j, k, n_max, choice=99;
   FILE *ascii_file;
   char filename[MLINE], label[MLINE];
   BOOLEAN found, print;

   while (choice != QUIT) {
     out_d(_("DATA MANAGEMENT: \n\n") );
     out_d(_("   0 = Main menu\n") );
     out_d(_("   1 = List data of columns\n") );
     out_d(_("   2 = Read another file\n") );
     out_d(_("   3 = List names of columns\n") );
     out_d(_("   4 = Rename column\n") );
     out_d(_("   5 = Read column from terminal\n") );
     out_d(_("   6 = Export columns as ASCII-data\n") );
     out_d(_("   7 = Export data base as fixed width data file\n") );
     out_d(_("   8 = File format options\n") );
     out_d(_("\n  Your choice: ") );
     GETRLINE;
     status = sscanf(line,"%i", &choice);
     if ((status==0) || (empty)) {
       choice = 99;
     }
     out_d("\n\n");

     switch(choice) {
       case 0:
               break;
       case 1:
	       printcols();
	       break;
       case 2:
               newsourcefile();
               break;
       case 3:
	       if(first_labels){
		 out_i(_("Include value labels? (%s) "), _("y/N"));
		 GETNLINE;
		 if (!(empty) && (line[0] == _("y")[0] || line[0] == _("Y")[0])) /* Use out_r because */
		   print = TRUE;	   /* this function is useful to */
		 else			   /* create a sub-set of the labels */
		   print = FALSE;	   /* file */
		 for(i = 0; i < ncol; i++){
		   if(names[i] && names[i]->ctitle){
		     out_r(_("Column %2i: "), i+1);
		     if(print)
		       out_r("\n");
		     colorize(ClHeader);
		     out_r("%-10s", alias[i]);
		     colorize(ClDefault);
		     out_r(" %s\n", names[i]->ctitle);
		   }
		   else
		     out_r(_("Column %2i = %s\n"), (i+1), alias[i]);
		   if(print && names[i]){
		     for(j = 0; j < names[i]->n; j++)
		       out_r("      %6g %s\n", names[i]->v[j], names[i]->l[j]);
		     out_r("\n");
		   }
		 }
	       } else
		 for(i = 0; i < ncol; i++)
		   out_d(_("Column %2i = %s\n"), (i+1), alias[i]);
	       out_d("\n");
	       break;
       case 4:
	       out_d(_("Columns: ") );
	       for (k=0; k<ncol; k++) {
		 out_d("%s ", alias[k]);
	       }
	       out_d("\n\n");

               out_i(_("Please enter name of column: ") );
               GETBLINE;
	       found = FALSE;
	       sscanf(line, "%s", label);
	       for (i=0; i<ncol; i++) {
		 if (strcmp(line, alias[i])==0) {
		   found = TRUE;
		   out_i(_("Please enter name for the column: ") );
		   GETBLINE;
		   myfree(alias[i]);
		   sscanf(line, "%s", label);
		   alias[i] = mymalloc(strlen(label)+1);
		   strcpy(alias[i], label);
		   break;
		 }
	       }
	       if (!found) {
		 out_err(ERR, ERR_FILE, ERR_LINE,
		     _("Column name \"%s\" doesn't exist!"), label);
	       }
               break;
       case 5:
               readcol_from_term();
               break;

       case 6:
	       i = getcols(1, ncol, 3);
	       n_max=0;
	       for (k=0; k<i; k++) {
		 if (nn[acol[k]] > n_max) {
		   n_max = nn[acol[k]];
		 }
	       }
	       if(i != 0) {
		 out_i(_("Please enter name of the export file: ") );
		 GETBLINE;
		 sscanf(line, "%s", filename);
		 FOPEN(filename, "wt", ascii_file);
#ifndef NO_GETTEXT
		 SET_C_LOCALE;
#endif
		 /* begin writing column labels */
		 if(has_header || detect_header)
		   fprintf(ascii_file, "%s", alias[acol[0]]);
		 else
		   fprintf(ascii_file, "#%%%s", alias[acol[0]]);
		 for(j=1; j<i; j++)
		   fprintf(ascii_file, " %s", alias[acol[j]]);
		 fprintf(ascii_file, "\n"); /* end writing column labels */
		 if(int_as_int){
		   for (k=0; k<n_max; k++) {
		     for (j=0; j<i; j++) {
		       if (nn[acol[j]]>k) {
			 if(xx[acol[j]][k] == SYSMIS)
			   fprintf(ascii_file, " %s", NODATA);
			 else
			   fprintf(ascii_file, " %g", xx[acol[j]][k]);
		       }
		       else {
			 fprintf(ascii_file, " %s", NODATA);
		       }
		     }
		     fprintf(ascii_file, "\n");
		   }
		 } else{
		   for (j=0; j<i; j++) {
		     if (nn[acol[j]]>k) {
		       if(xx[acol[j]][k] == SYSMIS)
			 fprintf(ascii_file, " %s", NODATA);
		       else
			 fprintf(ascii_file, " %10.5e", xx[acol[j]][k]);
		     }
		     else {
		       fprintf(ascii_file, " %s", NODATA);
		     }
		   }
		   fprintf(ascii_file, "\n");
		 }
		 FCLOSE(ascii_file);
#ifndef NO_GETTEXT
		 RESET_LOCALE;
#endif
		 out_d(_("Created file %s with %i columns!\n"), filename, i);
	       }
	       break;
       case 7:
	       exp_fwdf();
	       break;
       case 8:
	       if(sourcename)
		 show_file_head(sourcename);
	       set_fileformat();
	       break;
       default:
	       out_err(ERR, ERR_FILE, ERR_LINE, _("Illegal instruction!") );
         break;
       }
     m_freeall(); /* deallocate used memory for tmp variables */
     if (choice != 0) {
       mywait();
     }
   }
   return;
} /* data_menu() */



/* =================================================================== */

void manipulate_menu() {
   int i, j, k, choice=99;

   while (choice != QUIT) {
     out_d(_("DATA MANIPULATION: \n\n") );
     out_d(_("   0 = Main menu\n") );
     out_d(_("   1 = Log-transformation (base 10)\n") );
     out_d(_("   2 = Invert values (1/x)\n") );
     out_d(_("   3 = z-transformation [(x-mu)/sdv]\n") );
     out_d(_("   4 = Sort\n") );
     out_d(_("   5 = Join columns\n") );
     out_d(_("   6 = Exponentiation to base 10\n") );
     out_d(_("   7 = Create columns for weighted mean\n") );
     out_d(_("   8 = Log-transformation (natural logarithm)\n") );
     out_d(_("   9 = Exponentiation to base 'e'\n") );

     out_d(_("\n  Your choice: ") );
     GETRLINE;
     status = sscanf(line,"%i", &choice);
     if ((status==0) || (empty)) {
       choice = 99;
     }
     out_d("\n\n");

     switch(choice) {
       case 0:
               break;
       case 1:
               log_transform();
               break;
       case 2:
               inv_transform();
               break;
       case 3:
               z_transform();
               break;
       case 4:
	       sort_col();
	       break;
       case 5:
	       i = getcols(2, ncol, FALSE);
               if (i != 0) {
		 create_columns(1);
                 nn[ncol - 1]= 0;
                 for (j=0; j<i; j++) {
                   out_d(_("Number of values in column %s: %i\n"),
		       alias[acol[j]], nn[acol[j]]);
                   nn[ncol - 1] += nn[acol[j]];
                   for (k=0; k<nn[acol[j]]; k++) {
                     FWRITE(&(xx[acol[j]][k]), sizeof(REAL), 1, tmpptr[ncol -1]);
                   }
                 }
                 out_d(_("\nCreated column %s with %i values!\n"),
		     alias[ncol - 1], nn[ncol - 1]);
               }
	       break;

	case 6:
	       power_10_transform();
	       break;
        case 7:
               out_d(_("Columns: ") );
	       for (k=0; k<ncol; k++) {
		 out_d("%s ", alias[k]);
	       }
	       out_d("\n\n");
	       out_i(_("Please select column with values and column with "
		     "factors:\n") );
	       i = getcols(2, 2, TRUE);
	       if (i != 0) {
		 create_columns(1);
		 nn[ncol - 1] = 0;
		 for (i=0; i<nn[acol[1]]; i++) {
		   for (j=0; j<(int)xx[acol[1]][i]; j++) {
		     FWRITE(&(xx[acol[0]][i]), sizeof(REAL), 1, tmpptr[ncol - 1]);
		   }
		   nn[ncol - 1] += (int)xx[acol[1]][i];
		 }
		 out_d(_("\nCreated column %s with %i values!\n"),
		     alias[ncol - 1], nn[ncol - 1]);
	       }
	 break;
	case 8:
		ln_transform();
		break;

	case 9:
		power_e_transform();
		break;

       default:
		out_err(ERR, ERR_FILE, ERR_LINE, _("Illegal instruction!") );
         break;
     }
     m_freeall(); /* deallocate used memory for tmp variables */
     if (choice != 0) {
       mywait();
     }
   }
   return;
} /* manipulate_menu() */


/* =================================================================== */

void prefs_menu(){
  int i, choice = 99;
  char filename[MLINE];

  while (choice != QUIT) {
    set_winsize();
    out_d(_("PREFERENCES: \n\n") );
    out_d(_("   0 = Main menu\n") );
    out_d(_("   1 = Save preferences\n"));
    out_d(_("   2 = Verbose"));
    if(verbose)
      out_d(_(" [yes]\n"));
    else
      out_d(_(" [no]\n"));
    out_d(_("   3 = Gnuplot graphics") );
    if(noplot)
      out_d(_(" [no]\n"));
    else
      out_d(_(" [yes]\n"));
    out_d(_("   4 = Beep at errors and warnings") );
    if(nobell)
      out_d(_(" [no]\n"));
    else
      out_d(_(" [yes]\n"));
    out_d(_("   5 = Histogram as text graphic instead of gnuplot-graphic") );
    if(thist)
      out_d(_(" [yes]\n"));
    else
      out_d(_(" [no]\n"));
    out_d(_("   6 = Special output changes from Bernhard") );
    if(bernhard)
      out_d(_(" [yes]\n"));
    else
      out_d(_(" [no]\n"));
    out_d(_("   7 = Use system command \"%s\""), ls_cmd);
    if(system_ls)
      out_d(_(" [yes]\n"));
    else
      out_d(_(" [no]\n"));
    out_d(_("   8 = Use value labels") );
    if(first_labels)
      out_d(_(" [yes]\n"));
    else
      out_d(_(" [no]\n"));
    out_d(_("   9 = Maximum number of rows before aborting table printing"
	  " [%i]\n"), MRESULT);
    out_d(_("  10 = Screen number of columns [%i]\n"), SCRCOLS);
    out_d(_("  11 = Screen number of lines [%i]\n"), SCRLINES);
 
    out_d(_("\n  Your choice: ") );
    GETRLINE;
    status = sscanf(line,"%i", &choice);
    if ((status==0) || (empty)) {
      choice = 99;
    }
    out_d("\n\n");

    switch(choice) {
      case 0:
	break;
      case 1:
	save_prefs();
	break;
      case 2:
	if(verbose)
	  verbose = FALSE;
	else
	  verbose = TRUE;
	break;
      case 3:
	if(noplot)
	  noplot = FALSE;
	else
	  noplot = TRUE;
	break;
      case 4:
	if(nobell)
	  nobell = FALSE;
	else
	  nobell = TRUE;
	break;
      case 5:
	if(thist)
	  thist = FALSE;
	else
	  thist = TRUE;
	break;
      case 6:
	if(bernhard)
	  bernhard = FALSE;
	else
	  bernhard = TRUE;
	break;
      case 7:
	if(system_ls)
	  system_ls = FALSE;
	else
	  system_ls = TRUE;
	break;
      case 8:
	if(first_labels){
	  old_first_labels = first_labels;
	  first_labels = NULL;
	  for(i = 0; i < ncol; i++)
	    names[i] = NULL;
	} else{
	  if(old_first_labels)
	    first_labels = old_first_labels;
	  else{
	    out_i(_("No file with labels was loaded yet.\n"
		  "Do you want to load one now? (%s) "), _("y/N"));
	    GETNLINE;
	    if (!(empty) && (line[0] == _("y")[0] || line[0] == _("Y")[0])){
	      ls();
	      out_i(_("Name of file with labels: "));
	      GETBLINE;
	      sscanf(line, "%s", filename);
	      if(myexist(filename))
		read_labels(filename);
	      else
		out_err(ERR, ERR_FILE, ERR_LINE,
		    _("File \"%s\" not found!"), filename);
	    }
	  }
	  attach_labels_to_columns();
	}
	break;
      case 9:
	out_i(_("Please, choose a new value: "));
	GETBLINE;
	i = getint();
	if(i == -1)
	  break;
	if(i < 5){
	  out_err(ERR, ERR_FILE, ERR_LINE,
	      _("The number must be bigger than %i\n"), 10);
	  mywait();
	  break;
	}
	MRESULT = i;
	break;
      case 10:
	out_i(_("Please, choose a new value: "));
	GETBLINE;
	i = getint();
	if(i == -1)
	  break;
	if(i < 15){
	  out_err(ERR, ERR_FILE, ERR_LINE,
	      _("The number must be bigger than %i\n"), 15);
	  mywait();
	  break;
	}
	rcols = i;
	set_winsize();
	break;
      case 11:
	out_i(_("Please, choose a new value: "));
	GETBLINE;
	i = getint();
	if(i == -1)
	  break;
	if(i < 10){
	  out_err(ERR, ERR_FILE, ERR_LINE,
	      _("The number must be bigger than %i\n"), 10);
	  mywait();
	  break;
	}
	rlines = i;
	set_winsize();
	break;
      default:
	out_err(ERR, ERR_FILE, ERR_LINE, _("Illegal instruction!"));
	break;
    }
  }
  return;
} /* prefs_menu() */

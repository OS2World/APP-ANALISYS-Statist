/* This file is part of statist
**
** It is distributed under the GNU General Public License.
** See the file COPYING for details.
**
**
** (c) 1997 Dirk Melcher
**  Doerper Damm 4
**  49134 Wallenhorst
**  GERMANY
**  Tel. 05407/7636
**  email: Dirk.Melcher@usf.Uni-Osnabrueck.DE
**
**  some changes by Bernhard Reiter  http://www.usf.Uni-Osnabrueck.DE/~breiter
**  $Id: data.h,v 1.16 2006/09/09 18:12:38 jakson Exp $
***************************************************************/

/* data.h for STATIST */
#include <stdio.h>

extern int   getcols(int min, int max, BOOLEAN eraserow);
extern void  alloc_cols(int n_alloc, BOOLEAN eraserow);
extern PREAL readcol(int i);
extern void  readsourcefile();
extern void  read_labels(char *labelsfile);
extern void  attach_labels_to_columns();
extern void  newsourcefile();
extern void  create_columns(int amount);
extern void delete_last_columns(int i);
extern void  erasetempfiles();
extern int   parsecomment(const char *comment, BOOLEAN is_comment);
extern void  printcol(REAL x[], int n);
extern void  printcols();
extern BOOLEAN make_new_col(char *alias, int n);
extern void  log_transform();
extern void  power_10_transform();
extern void  ln_transform();
extern void  power_e_transform();

extern void  inv_transform();
extern void  z_transform();
extern void  sort_col();
extern void  readcol_from_term();
extern BOOLEAN str_in_str(const char *s1, const char *s2);
extern char *get_label(PREAL x);
extern char *get_name(PREAL x);
extern BOOLEAN emptyline(const char *s);
int col_exist(char *analias, BOOLEAN is_error);
extern void set_fileformat();
extern void show_file_head(char *fn);
void extract_sample(int argc, char * argv[]);
void extract_cols(int argc, char * argv[]);
void exp_fwdf();


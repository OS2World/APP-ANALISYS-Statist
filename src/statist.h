/* This file is part of statist
**
** It is distributed under the GNU General Public License.
** See the file COPYING for details.
**
** (c) 1997 Dirk Melcher
** old email address: Dirk.Melcher@usf.Uni-Osnabrueck.DE
**
** adaptions for StatistX: Andreas Beyer, 1999, abeyer@usf.uni-osnabrueck.de
**
** published by Bernhard Reiter  http://www.usf.Uni-Osnabrueck.DE/~breiter
** $Id: statist.h,v 1.41 2009/12/10 13:17:03 jakson Exp $
***************************************************************/

#include <stdio.h>

#include "gettext.h"
#include "memory_handling.h"


/* next two defines are important for gettext */
#define PACKAGE "statist"
/* LOCALEDIR should be set by Makefile, this is a fallback	*/
#ifndef LOCALEDIR
#  define LOCALEDIR "/usr/share/locale"
#endif


#define _STATIST_VERSION_NUMBER "1.4.2"
#define VERSION_INFO N_("\nSTATIST --- Version %s\n " \
"(C) Dirk Melcher 1997-1999\n \
(C) Bernhard Reiter 1998-2009\n \
(C) Jakson Aquino 2005-2009\n \
STATIST comes with ABSOLUTELY NO WARRANTY; it is Free Software under GNU GPL.\n\
 Read the file COPYING for details.\n\n")

#define STATIST_VERSION "statist -- version " _STATIST_VERSION_NUMBER

#ifdef MSDOS
  /* change the MSDOS_FIXED_TMP_FILE, if you want statist to create
	the temporary files at a different place on MSDOS
	"%i"  must be in the string once, to replaced by an (integer) number
   */
#ifndef MSDOS_FIXED_TMP_FILE
#define MSDOS_FIXED_TMP_FILE "c:\\tmp\\stat%i.tmp"
#endif
#endif


/* others may define TRUE & FALSE, too! */
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define CRASH -1
#define CHAR_OFFSET 97  /* 'a' position in ASCII-Code */
#define QUIT 0
#define COMMENT '#'  /* Symbol that indicates commentary in data files */
#define _ALL_ N_("all")  /* alias for: "all columns" */
#define MCLASS 60    /* Max. number of classes on normal distribution test */
#define MPOLY 20     /* Maximum degree of polynomial */
#define MOPT  24     /* Total number of options */
#define MLINE 255    /* Maximum number of chars for standard line variable */

#define WAR 0	/* Error - Warning */
#define ERR 1	/* Normal Error */
#define FAT 2	/* Fatal Error */
#define MAT 3	/* Math Error */
#define MWA 4	/* Math Warning */

#define SQR(x) ((x) * (x))

#define GETRLINE  fgets((line), (MLINE-1), (stdin));	\
		  if (strlen(line) <=1) {		\
		    empty = TRUE;			\
		    return;	       			\
		  }					\
		  else {				\
		    line[strlen(line)-1]='\0';		\
		    empty = FALSE;			\
		  }

#define GETNLINE  fgets((line), (MLINE-1), (stdin));	\
		  if (strlen(line) <=1) {		\
		    empty = TRUE;			\
		  }					\
		  else {				\
		    line[strlen(line)-1]='\0';		\
		    empty = FALSE;			\
		  }

#define GETBLINE fgets((line), (MLINE-1), (stdin));	\
		  if (strlen(line) <=1) {		\
		    empty = TRUE;			\
		    break;				\
		  }					\
		  else {				\
		    line[strlen(line)-1]='\0';		\
		    empty = FALSE;			\
		  }

#ifdef SUNOS
extern char *sys_errlist[];
#define STRERROR(errn) sys_errlist[errn]
#else
#define STRERROR(errn) strerror(errn)
#endif

#ifdef DEBUG
#define ERR_LINE __LINE__
#define ERR_FILE __FILE__
#else
#define ERR_LINE 0
#define ERR_FILE ""
#endif


#define FOPEN(name, mode, fp)					\
  if ((fp = fopen((name), (mode)))==NULL) {			\
    out_err(FAT, ERR_FILE, ERR_LINE,				\
    _(" System reports error while opening file %s:\n   %s"),	\
       name, STRERROR(errno));					\
  }

#define FWRITE(ptr, size, nmemb, fp)				\
  if (fwrite((ptr), (size), (nmemb), (fp)) != (nmemb)) {	\
    out_err(FAT,ERR_FILE, ERR_LINE,				\
    _(" System reports error while writing with fwrite:\n %s"),	\
	      STRERROR(errno));					\
  }

#define FREAD(ptr, size, nmemb, fp)					\
  if (fread((ptr), (size), (nmemb), (fp)) != (nmemb)) {			\
    if (feof(fp)) {							\
      out_err(FAT, ERR_FILE, ERR_LINE,					\
    _(" Error while reading with fread: Unexpected end of file") );	\
    }									\
    else {								\
      out_err(FAT, ERR_FILE, ERR_LINE,					\
	  _("System reports error while reading with fread:\n %s"),	\
	      STRERROR(errno));						\
    }									\
  }

#define FGETS(line, size, fp)						\
  if (fgets((line), (size), (fp))==NULL) {				\
    if (feof(fp)) {							\
      out_err(FAT, ERR_FILE, ERR_LINE,					\
	_("   Error while reading with fgets: Unexpected end of file\n") );\
    }									\
    else {								\
      out_err(FAT, ERR_FILE, ERR_LINE,					\
	_(" System reports error while reading with fgets:\n   %s"),	\
	      STRERROR(errno));						\
    }									\
  }

#define FCLOSE(fp)							\
  if (fclose(fp) != 0) {						\
    out_err(ERR, ERR_FILE, ERR_LINE,					\
    _("System reports error while attempting to close file:\n  %s"),	\
	    STRERROR(errno));						\
  }

#define DIV(r, c, d)						\
   if ((d)==0.0) {						\
     out_err(MAT, ERR_FILE, ERR_LINE, _("Division by 0!") );	\
     return;							\
   }								\
   else {							\
     r = (c)/(d);						\
   }

#define RDIV(r, c, d)						\
   if ((d)==0.0) {						\
     out_err(MAT, ERR_FILE, ERR_LINE, _("Division by 0!") );	\
     return REAL_MIN;						\
   }								\
   else {							\
     r = (c)/(d);						\
   }

#define SQRT(y, x)						\
  if ((x) < 0.0) {						\
     out_err(MAT, ERR_FILE, ERR_LINE,				\
     	_("Square root with negative argument!") );		\
     return;							\
  }								\
  else {							\
     y = sqrt(x);						\
  }

#define RSQRT(y, x)						\
  if ((x) < 0.0) {						\
     out_err(MAT, ERR_FILE, ERR_LINE,				\
     	_("Square root with negative argument!") );		\
     return REAL_MIN;						\
  }								\
  else {							\
     y = sqrt(x);						\
  }


#define REAL_EPSILON DBL_EPSILON
#define REAL_MAX     DBL_MAX
#define REAL_MIN    (-DBL_MAX)

typedef short int BOOLEAN;
typedef double REAL;
typedef REAL* PREAL;

typedef struct {
  REAL val, rank;
  int ind;
} SORTREC;


#include <errno.h>
extern void mywait();
extern BOOLEAN myexist(char *name);
#ifdef STATIST_X
extern void out_start(); /* start results output  */
extern void out_end();   /* end results output    */
#else
#	define out_start()
#	define out_end()
#endif
extern void out_d(char *fmt, ...);
extern void out_i(char *fmt, ...);
extern void out_r(char *fmt, ...);
extern void out_err(int errn, char *modul, int lno, char *fmt, ...);
extern void print_histo(REAL x[], int y[], int n);
extern void finish();
extern BOOLEAN ls();
extern void save_prefs();

extern void set_winsize();
extern size_t stringLen(char *s);
extern int SCRLINES, SCRCOLS, rcols, rlines, MRESULT;
extern REAL SYSMIS; /* 'M's are put in tempfiles with this value */
extern char *NODATA;   /* Symbol that indicates missing values in data files */
extern char *sourcename;
extern char sep; /* Field separator in csv files */
extern char dec; /* Decimal delimiter in data files */
extern PREAL  *yy, tempcol;
extern int MCOL;  /* Max. number of columns that can be loaded before
		     reallocating memory */
extern int *ny;

typedef struct { 
  char *clabel; /* The column label as it appears in the datafile */
  char *ctitle; /* A longer name for the column (can contain spaces) */
  int n;	/* Number of value labels for this column */
  REAL *v;	/* List of values with labels */
  char **l;	/* List of value labels */
  void *next;
} Labels;

extern Labels **names; /* Pointer to labels for specific columns */
extern Labels *first_labels, /* The first one in a linked list of labels */
       *old_first_labels; /* Will store first_labels if the user chooses do not
			     show labels */

/* "struct" COLUMN */
extern PREAL  *xx;      /* column data */
extern char **alias;    /* column label */
extern FILE **tmpptr;  /* pointers to temporary files */
extern int *nn,	 /* number of data points in tempfile (total) */
    *vn;	 /* number of (valid) data points in ram */
extern BOOLEAN *x_read; /* is the column already in ram? */

extern BOOLEAN empty, gnupl_open, verbose;
extern int     ncol, /* number of columns in the entire database */
	       *acol, /* list of columns currently in ram */
	       status;
extern char    line[MLINE];
extern FILE    *logfile;
extern BOOLEAN silent, noplot, help, log_set, nofile, nobell, thist, bernhard,
       has_header, noheader, color, quoted;

extern FILE   *pipef;

extern char *ls_cmd, *gtprefix, *gplt_default_term, *gplt_png;
extern BOOLEAN system_ls, format_columns_out, detect_header, ask_fileformat, int_as_int;
extern BOOLEAN is_utf8; /* Is the user locale set to UTF-8 ?*/
#ifdef MSDOS
extern char msdos_temp_dir[255];
#endif
/* These macros are used because Gnuplot has problems to
   deal with other locales than the standard C-locale. */
#ifndef NO_GETTEXT
extern char *saved_locale, *old_locale;
#define SET_C_LOCALE	\
	old_locale = setlocale (LC_NUMERIC, NULL); \
	saved_locale = mymalloc(strlen(old_locale)+1); \
	strcpy(saved_locale, old_locale); \
	setlocale(LC_NUMERIC, "C");

#define RESET_LOCALE \
	setlocale (LC_NUMERIC, saved_locale); \
	free(saved_locale); saved_locale = NULL;
#endif


void colorize(int cl);
#define ClDefault 0
#define ClError 1
#define ClErrorD 2
#define ClInstr 3
#define ClHeader 4
#define ClLineNum 5
#define ClMenuSep 6


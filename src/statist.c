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
**  $Id: statist.c,v 1.41 2009/12/21 17:08:51 jakson Exp $
***************************************************************/

/* statist.c */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <signal.h>
#include <dirent.h>
#ifndef NO_VALUES_H
#include <values.h>
#endif
#ifndef NO_IOCTIL_H
#include <sys/ioctl.h> /* terminal number of columns and lines */
#endif
#ifndef MSDOS
#include <unistd.h> /* getlogin */
#include <pwd.h> /* getpwnam */
#endif

#include "statist.h"
#include "menue.h"
#include "data.h"
#include "funcs.h"
#include "plot.h"
#include "gettext.h"

typedef struct {
  char *optstr;
  BOOLEAN *optvar;
} OPTREC;

BOOLEAN get_options(char **args, int options[], int nopt, OPTREC optrec[]);
BOOLEAN read_options();
void help_msg();
void set_prefs_filename();

#ifndef MSDOS
void handle_pipe_interrupt(int sig);
char * user_dir();
#endif
void handle_fpe_interrupt(int sig);
void handle_term_interrupt(int sig);
void set_default_colors();

Labels **names;
Labels *first_labels, *old_first_labels;
PREAL   *xx, *yy, tempcol;
#ifndef NO_VALUES_H
REAL SYSMIS = (-DBL_MAX);
#else
REAL SYSMIS = (1.79769313486231570e-308);
#endif
char  *sourcename;
char *NODATA = NULL;
char sep = 0;
char dec = '.';
int     *nn, *vn, *acol, *ny;
int MRESULT = 200;

int SCRCOLS, SCRLINES; /* terminal number of cols and lines */
int rcols = 0, rlines = 0; /* values for SCRCOLS and SCRLINES read from
			      statistrc (or set in Preferences menu) */

BOOLEAN is_utf8;
int MCOL = 0;
BOOLEAN *x_read, empty, gnupl_open=FALSE;
int     ncol=0, status;
char    line[MLINE];
char    **alias; /* Column labels*/
FILE    **tmpptr; /* Pointers to temporary files */
FILE    *logfile;

BOOLEAN silent=FALSE, noplot=FALSE, help=FALSE, log_set=FALSE, nofile=FALSE,
	nobell=FALSE, thist=FALSE, bernhard=FALSE, xsample=FALSE, xcols=FALSE,
	color=FALSE, version=FALSE, noheader=FALSE, has_header=FALSE;
BOOLEAN system_ls = FALSE, verbose = TRUE, format_columns_out = TRUE,
	detect_header = TRUE, ask_fileformat = FALSE, int_as_int = FALSE;

FILE    *pipef;

char msdos_temp_dir[255];
char *ls_cmd; /* Command send to system to list content of current directory */
char *saved_locale, *old_locale, *prefs_filename, *gtprefix, *gplt_default_term,
     *gplt_png;
char *clError, *clInstr, *clHeader, *clLineNum,  *clMenuSep;

/* ================================================================== */

#ifndef NO_GETTEXT
BOOLEAN get_utf8_status(){
  char *s;
  s = setlocale(LC_CTYPE, NULL);
  if(strstr(s, "UTF-8") != NULL)
    return TRUE;
  else
    return FALSE;
}

size_t stringLen(char *s){
  wchar_t dest[256];
  if(is_utf8)
    return mbstowcs(dest, s, 255);
  else
    return strlen(s);
}
#else
size_t stringLen(char *s){
  return strlen(s);
}
#endif


int main(int argc, char *argv[]) {
   int   i, labelsfile = 0;
   BOOLEAN source_set=FALSE;
   BOOLEAN has_strc;
   int options[MOPT], nopt=0;
#ifdef MSDOS
   char *s;
   s = getenv("TMP");
   if(strlen(s))
     strcpy(msdos_temp_dir, s);
   else
     strcpy(msdos_temp_dir, "C:\\tmp");
   set_prefs_filename(argv[0]);
   has_strc = FALSE;
#endif

   NODATA = (char*)mymalloc(2 * sizeof(char));
   strcpy(NODATA, "M");

   OPTREC optrec[MOPT] = {
      {"-h", &help}, {"-help", &help}, {"--help", &help}, {"-?", &help},
      {"-silent", &silent}, {"--silent", &silent},
      {"-noplot", &noplot}, {"--noplot", &noplot},
      {"-log", &log_set}, {"--log", &log_set},
      {"-nofile", &nofile}, {"--nofile", &nofile},
      {"-nobell", &nobell}, {"--nobell", &nobell},
      {"-thist", &thist}, {"--thist", &thist},
      {"--bernhard",&bernhard},
      {"--header", &has_header},
      {"--noheader", &noheader},
      {"--xsample", &xsample},
      {"--xcols", &xcols},
      {"--color", &color},
      {"--quiet",&silent},
      {"--version", &version} };

#ifndef MSDOS
   prefs_filename = (char*)mymalloc(15 * sizeof(char));
   strcpy(prefs_filename, "/etc/statistrc");
   has_strc = read_options();
   myfree(prefs_filename);
   prefs_filename = NULL;
   set_prefs_filename();
   signal(SIGPIPE, handle_pipe_interrupt);
#endif
   signal(SIGFPE, handle_fpe_interrupt);
   signal(SIGINT, handle_term_interrupt);
   signal(SIGTERM, handle_term_interrupt);

#ifndef NO_GETTEXT
   setlocale (LC_ALL, "");
#ifdef MSDOS
   char locale_dir[255];
   strcpy(locale_dir, prefs_filename);
   int k = strlen(locale_dir);
   for(i = (k-1); i > 0; i--)
     if(locale_dir[i] == '/'){
       locale_dir[i] = 0;
       break;
     }
   strcat(locale_dir, "/locale");
   bindtextdomain (PACKAGE, locale_dir);
#else
   bindtextdomain (PACKAGE, LOCALEDIR);
#endif
   textdomain (PACKAGE);
   is_utf8 = get_utf8_status();
#endif

   if(has_strc)
     read_options();
   else
     has_strc = read_options();

   set_default_colors();
   for (i=1; i<argc; i++) {
     if (argv[i][0]=='-') {
       options[nopt] = i;
       nopt ++;
       if (nopt > MOPT) {
	 out_err(FAT, ERR_FILE, ERR_LINE, _("Too many options") );
       }
     }
   }

   if (!get_options(argv, options, nopt, optrec)) {
     exit(1);
   }
   out_d(_(VERSION_INFO), _STATIST_VERSION_NUMBER);
#ifdef DEBUG
   out_d(_(" Compiled with DEBUG.\n"));
#ifndef NO_GETTEXT
  #ifdef MSDOS
    out_d(" LOCALEDIR=%s\n", locale_dir);
  #else
    out_d(" LOCALEDIR=%s\n", LOCALEDIR);
  #endif
#endif
#endif

  if(version){
    exit(0);
  }
  if(silent)
    verbose = FALSE;

  if(noheader)
    detect_header = FALSE;

  if(verbose){
    out_d(_(" Hint: You can often answer \"%s\" for the number of columns.\n"),
	_(_ALL_) );
#ifndef MSDOS
    if(has_strc == FALSE)
      out_d(_("\n You can create an ~/.statistrc file to increase statist "
	    "functionality.\n There is an example of statistrc at %s.\n"),
	  DOCDIR);
#endif
  }

   out_d("\n");
   if (help) {
     help_msg();
     exit(0);
   }

   if(xsample){
     extract_sample(argc, argv);
     finish();
     exit(0);
   }

   for (i=1; i<argc; i++) {
     if (!(argv[i][0]=='-')){
       if(strcmp(argv[i-1], "--labels") == 0){
	 labelsfile = i;
       } else
	 if(strcmp(argv[i-1], "--na-string") == 0){
	   NODATA = (char*)mymalloc((strlen(argv[i]) + 1) * sizeof(char));
	   strcpy(NODATA, argv[i]);
	 } else
	   if(strcmp(argv[i-1], "--sep") == 0){
	     if(strcmp(argv[i], "\\t") == 0)
	       sep = '\t';
	     else
	       sep = argv[i][0];
	   } else
	     if(strcmp(argv[i-1], "--dec") == 0){
	       dec = argv[i][0];
	     } else
	       if(!xcols){
		 if (!source_set) {
		   sourcename = (char*) mycalloc(strlen(argv[i]) + 1, sizeof(char));
		   strcpy(sourcename, argv[i]);
		   source_set = TRUE;
		 }
		 else {
		   if(silent)
		     out_err(FAT, ERR_FILE, ERR_LINE,
			 _(" More than one data file given!"));
		   else{
		     out_err(ERR, ERR_FILE, ERR_LINE,
			 _(" More than one data file given!\n"
			   "   Only the first one will be read."));
		     if(verbose)
		       out_d(_(
			     "   If you really want to load more than one file, please\n"
			     "   choose the menu option:\n"
			     "   'Data management | Read another file'\n"));
		     mywait();
		   }
		 }
	       }
     }
   }

   if(xcols){
     extract_cols(argc, argv);
     finish();
     exit(0);
   }

   if ( (!source_set) && (!nofile) ) {
     ls();
     out_i(_("Name of data file: ") );
     GETNLINE
       if (empty) {
	 out_d(_("No data file given. "));
	 nofile = TRUE;
       } else {
	 sourcename = (char*) mycalloc(strlen(line) + 1, sizeof(char));
	 sscanf(line, "%s", sourcename);
	 if(!(myexist(sourcename))){
	   nofile = TRUE;
	   out_err(ERR, ERR_FILE, ERR_LINE,
	       _("File \"%s\" not found!"), sourcename);
	   myfree(sourcename);
	   sourcename = NULL;
	 }
       }
   }

   if (log_set) {
     FOPEN("statist.log", "wt", logfile);
     if(!nofile)
       out_r(_("Source file = %s\n\n"), sourcename);
   }

   if(dec != '.' && dec != ','){
     if(silent){
       out_err(WAR, ERR_FILE, ERR_LINE,
	   _("Invalid decimal delimiter: '%c'. Please, choose either ',' or '.'"), dec);
     } else{
       out_err(ERR, ERR_FILE, ERR_LINE,
	   _("Invalid decimal delimiter: '%c'. Please, choose either ',' or '.'"), dec);
       if(!nofile)
	 show_file_head(sourcename);
       set_fileformat();
     }
   }

   if(sep && !(sep == ' '  || sep == ',' || sep == ';' || sep == '\t')){
       if(silent){
	 out_err(WAR, ERR_FILE, ERR_LINE, _("Invalid field separator: '%c'"), sep);
       } else{
	 out_err(ERR, ERR_FILE, ERR_LINE, _("Invalid field separator: '%c'"), sep);
	 if(!nofile)
	   show_file_head(sourcename);
	 set_fileformat();
       }
     }

   if(nofile){
     out_d(_("Load file using menu option\n"
	   "Data management | Read another file\n"));
   } else{
     readsourcefile(sourcename);
     mywait();
   }
   out_d("\n");

   if(labelsfile){
     if(myexist(argv[labelsfile]))
       read_labels(argv[labelsfile]);
     else{
       out_err(ERR, ERR_FILE, ERR_LINE,
	   _("File \"%s\" not found!"),
	   argv[labelsfile]);
     }
   }

   main_menue();
   finish();
   return 0;
}

/* =================================================================== */

void finish() {
  int i;
  Labels *ptr;
   erasetempfiles();
   if (log_set) {
     FCLOSE(logfile);
   }
#ifndef NOPIPE
   if (gnupl_open == TRUE) {
     fprintf(pipef, "quit\n");
     pclose(pipef);
   }
   if (myexist(GPL_DAT)) {
     remove(GPL_DAT);
   }
#endif
   if(prefs_filename)
     myfree(prefs_filename);
   if(sourcename)
     myfree(sourcename);
   if(ls_cmd)
     myfree(ls_cmd);
   if(gnuplot_charset)
     myfree(gnuplot_charset);
   if(old_first_labels)
     first_labels = old_first_labels;
   while(first_labels){
     ptr = first_labels->next;
     myfree(first_labels->clabel);
     myfree(first_labels->ctitle);
     for(i = 0; i < first_labels->n; i++){
       myfree(first_labels->l[i]);
     }
     myfree(first_labels->v);
     myfree(first_labels->l);
     myfree(first_labels);
     first_labels = ptr;
   }
   if(plothist){
     for(i = 0; i < histsize; i++)
       if(plothist[i])
	 myfree(plothist[i]);
     myfree(plothist);
   }
   if(gplt_png)
     myfree(gplt_png);
   if(clError)
     myfree(clError);
   if(clInstr)
     myfree(clInstr);
   if(clHeader)
     myfree(clHeader);
   if(clLineNum)
     myfree(clLineNum);
   if(clMenuSep)
     myfree(clMenuSep);
   if(gtprefix)
     myfree(gtprefix);
   myfree(NODATA);
}



void help_msg() {
out_d(_("Usage:\n"
      "  statist [ options ] data_file\n"));
out_d(_("Options:\n"));
out_d(_("--help, -h, -?  : print this help message and exit\n"));
out_d(_("--silent, --quiet\n"
        "                : don't print menu etc. (for batch/script usage)\n"));
out_d(_("--log           : write results to log file `statist.log'\n"));
out_d(_("--nofile        : don't read a data file when starting the program\n"));
out_d(_("--nobell        : no beep at errors and warnings\n"));
out_d(_("--thist         : histogram as text graphic instead of gnuplot-graphic\n"));
out_d(_("--noplot        : no gnuplot-graphic\n"));
#ifndef MSDOS
out_d(_("--color         : colorized output\n"));
#endif
out_d(_("--labels <file> : read column titles and value labels from file\n"));
out_d(_("--header        : the file has column names in the first line\n"));
out_d(_("--noheader      : the file does not have column names\n"));
out_d(_("--sep <char>    : field separator character\n"));
out_d(_("--dec <char>    : decimal delimiter character (default: '.')\n"));
out_d(_("--na-string <string>\n"
"                : indicator of missing values (default: \"M\")\n"));
out_d(_("--xcols <conf_file> <orig_datafile> <new_data_file>\n"
"                : extract columns from a fixed width data file\n"));
out_d(_("--xsample <percentage> <database> <dest_file>\n"
"                : extract a sample of rows from a file\n"));
out_d(_("--bernhard      : special output changes from Bernhard, i.e.:\n"
"                  - table output at Miscellaneous/Standard deviation\n"
"                  - if --noplot defined no text histogram at\n"
"                       Miscellaneous/Standard deviation\n"));
out_d("\n");
out_d(_("Missing values must be indicated by '%s'\n"), NODATA);
#ifndef MSDOS
#ifdef DOCDIR
out_d("\n");
out_d(_("Statist documentation can be found at: %s\n"), DOCDIR);
#endif
#endif
}

/* =================================================================== */

BOOLEAN get_options(char **args, int options[], int nopt, OPTREC optrec[]) {
  int i, k;
  BOOLEAN found;

  for (i=0; i<nopt; i++) {
    found = FALSE;
    for (k=0; k<MOPT; k++) {
      if (strcmp(args[options[i]], optrec[k].optstr) == 0) {
	*(optrec[k].optvar) = TRUE;
	found = TRUE;
	if ((strlen(optrec[k].optstr) > 2) && optrec[k].optstr[0] == '-'
	    && optrec[k].optstr[1] != '-')
	  out_err(WAR, ERR_FILE, ERR_LINE,
	      _("The option short format %s is deprecated\n"
		"  and might be removed some time in the future.\n"
		"  Use -%s instead."), optrec[k].optstr, optrec[k].optstr);
	break;
      }
    }
    if (!found && !(strcmp(args[options[i]], "--labels") == 0 || 
	  strcmp(args[options[i]], "--na-string") == 0 || 
	  strcmp(args[options[i]], "--dec") == 0 || 
	  strcmp(args[options[i]], "--sep") == 0)) {
      out_err(ERR, ERR_FILE, ERR_LINE,
		_("Illegal option: '%s'. Try -h for help."), args[options[i]]);
      return FALSE;
    }
  }
  return TRUE;
}

/* =================================================================== */

void mywait() {
  if (!silent) {
    out_i(_("  --- Please, continue with <RETURN> ... ---") );
    GETRLINE
  }
}



BOOLEAN myexist(char *name) {
  FILE *F;
  F = fopen(name,"r");
  if (F) {
    fclose(F);
    return TRUE;
  }
  else {
    return FALSE;
  }
}


/* This function is not very good, but it's possible to use the system "ls".
 * It will be used only if for some reason the user preferred it.  */
BOOLEAN ls(){
  struct dirent* entry;
  DIR* dir;
  char s[10];
  int i, l = 1, nc;
  if(silent)
    return(TRUE);
  if(system_ls && ls_cmd != NULL){
    out_d("\n");
    out_d(_("Content of current directory:"));
    out_d("\n\n");
    i = system(ls_cmd);
    out_d("\n");
    if(i == 0)
      return(TRUE);
    else
      out_err(WAR, ERR_FILE, ERR_LINE,
	  _("Error while running the command \"%s\"."), ls_cmd);
  }
  dir = opendir(".");
  if (!dir) {
    perror("opendir");
    return FALSE;
  } else{
    set_winsize();
    out_d("\n");
    out_d(_("Content of current directory:"));
    out_d("\n");
    if(format_columns_out){
      while((entry = readdir(dir)) != NULL){
	if(entry->d_name[0] != '.'){
	  i = strlen(entry->d_name) + 1;
	  if(i > l)
	    l = i;
	}
      }
      sprintf(s, "%%-%is ", l); /* formatting the column width */
      nc = SCRCOLS / l; /* number of columns */
      rewinddir(dir);
      i = 0;
      while((entry = readdir(dir)) != NULL){
	if(entry->d_name[0] != '.'){
	  if((i % nc) == 0)
	    out_d("\n");
	  out_d(s, entry->d_name);
	  i++;
	}
      }
    } else{
      while((entry = readdir(dir)) != NULL)
	if(entry->d_name[0] != '.')
	  out_d("%s  ", entry->d_name);
    }
    if (closedir(dir) == -1) {
      perror("closedir");
      return FALSE;
    }
    out_d("\n\n");
  }
  return TRUE;
}

void out_r(char *fmt, ...) {
  va_list argptr;

  va_start(argptr, fmt);
  vprintf(fmt, argptr);
  va_end(argptr);
  if (log_set) {
    va_start(argptr, fmt);
    vfprintf(logfile, fmt, argptr);
    va_end(argptr);
  }
}


void out_d(char *fmt, ...) {
  va_list argptr;

  if (!silent) {
    va_start(argptr, fmt);
    vprintf(fmt, argptr);
    va_end(argptr);
  }
}

void out_i(char *fmt, ...) {
  va_list argptr;

  if (!silent) {
    va_start(argptr, fmt);
    colorize(ClInstr);
    vprintf(fmt, argptr);
    colorize(ClDefault);
    va_end(argptr);
  }
}



void out_err(int errn, char *modulname, int lno, char *fmt, ...) {
  va_list argptr;
  char where[128];

  va_start(argptr, fmt);
  switch (errn) {
    case WAR: strcpy(line, _("> statist-warning ") );
      break;
    case ERR: strcpy(line, _("> statist-error ") );
      break;
    case FAT: strcpy(line, _("\n> statist-fatal error ") );
      break;
    case MAT: strcpy(line, _("> statist-numerical error ") );
      break;
    case MWA: strcpy(line, _("> statist-warning ") );
      break;
  }

  if (lno != 0) {
    sprintf(where, _("(Module %15s, line %i):\n  "), modulname, lno);
    strncat(line, where, 254-strlen(line));
  }
  else {
    strncat(line, ":\n  ", 254-strlen(line));
  }
  strncat(line, fmt, 254-strlen(line));
  strncat(line, "\n\n", 254-strlen(line));

  if(!nobell)
    fprintf(stderr, "\a");

  colorize(ClError);
  if(errn == FAT && !(silent)){
    strncat(line, "> ",	254-strlen(line));
    strncat(line, _("Please, report bugs to: statist-list@intevation.de\n\n"),
	254-strlen(line));
  }
  vfprintf(stderr, line, argptr);
  if ( ((errn==MAT) || (errn==MWA)) && (log_set) ) {
    vfprintf(logfile, line, argptr);
  }
  colorize(ClErrorD);
  va_end(argptr);

  fflush(stderr);
  if (errn == FAT) {
    finish();
    exit(1);
  }
}



/* =================================================================== */

void print_histo(REAL x[], int y[], int n) {
  int maxclass;
  int   k, cols, i;
  char show[80];

  maxclass = get_maxint(y, n);
  for (i=0; i<n; i++) {
    cols = y[i]*60/maxclass;
    strcpy(show,"");
    for (k=0; k<cols; k++) {
      strcat(show, "*");
    }
    out_r("%9f %s %i\n", x[i], show, y[i]);
    if (((i+1) % (SCRLINES - 1)) == 0) {
      mywait();
    }
  }
}


/* =================================================================== */

#ifndef MSDOS
void handle_pipe_interrupt(int sig) {
  gnupl_open = FALSE;
  out_err(WAR, ERR_FILE, ERR_LINE,
	_("Broken pipe to gnuplot. Re-connecting.") );
  pclose(pipef);
}
#endif

void handle_fpe_interrupt(int sig) {
  out_err(FAT, ERR_FILE, ERR_LINE,
	_("Signal from operating system: floating point exception\n\
		Division by zero? Terminating!") );
}

void handle_term_interrupt(int sig) {
  out_d("\n\nstatist: ");
  out_d(_("Received signal SIGTERM (Ctrl-C?). Terminating!") );
  out_d("\nstatist: ");
  out_d(_("Please, report bugs to: statist-list@intevation.de\n\n"));
  finish();
  exit(1);
}

#ifdef MSDOS
void set_prefs_filename(char *p){
  char *s;
  int i, l;
  s = (char*)mymalloc(300 * sizeof(char));
   strcpy(s, p);
   l = strlen(s);
   i = 0;
   for(i = 0; i < l; i++)
     if(s[i] == '\\')
       s[i] = '/';
   for(i = (l-1); i > 0; i--)
     if(s[i] == '/'){
       s[i] = 0;
       break;
     }
   strcat(s, "/statistrc.txt");
#else
void set_prefs_filename(){
  char *s;
  struct passwd *uinf = NULL;
  s = getlogin();
  if(s == NULL)
    return;
  uinf = getpwnam(s);
  if(uinf == NULL)
    return;
  s = (char*)mymalloc((12 + strlen(uinf->pw_dir)) * sizeof(char));
  strcpy(s, uinf->pw_dir);
  strcat(s, "/.statistrc");
#endif
  prefs_filename = (char*)mymalloc((strlen(s) + 1) * sizeof(char));
  strcpy(prefs_filename, s);
  myfree(s);
}

void set_color(char *opt, char *value){
  char s[20], *v, *t;
  s[0] = 0;
  
  if(strstr(value, "underln") == value){
    sprintf(s, "\033[3;");
    v = value + 7;
  } else{
    if(strstr(value, "bright") == value){
      sprintf(s, "\033[1;");
      v = value + 6;
    } else{
      if(strstr(value, "blink") == value){
	sprintf(s, "\033[5;");
	v = value + 5;
      } else{
	if(strstr(value, "dim") == value){
	  sprintf(s, "\033[2;");
	  v = value + 3;
	} else{
	  sprintf(s, "\033[");
	  v = value;
	}
      }
    }
  }
  if(strcmp(v, "magenta") == 0){
    strcat(s, "35m");
  } else 
    if(strcmp(v, "yellow") == 0){
      strcat(s, "33m");
    } else 
      if(strcmp(v, "green") == 0){
	strcat(s, "32m");
      } else 
	if(strcmp(v, "cyan") == 0){
	  strcat(s, "36m");
	} else 
	  if(strcmp(v, "blue") == 0){
	    strcat(s, "34m");
	  } else 
	    if(strcmp(v, "white") == 0){
	      strcat(s, "37m");
	    } else 
	      if(strcmp(v, "black") == 0){
		strcat(s, "30m");
	      } else 
		if(strcmp(v, "red") == 0){
		  strcat(s, "31m");
		}

  if(strlen(s) < 5){
    out_err(WAR, ERR_FILE, ERR_LINE,
	_("Unknown color name \"%s\" in file \"%s\"."), value, prefs_filename);
    return;
  }
  t = (char*)mymalloc((strlen(s) + 1) * sizeof(char));
  strcpy(t, s);
  if((strcmp(opt, "cl_menu_separator") == 0)){
    if(clMenuSep)
      myfree(clMenuSep);
    clMenuSep = t;
  } else
    if((strcmp(opt, "cl_header") == 0)){
      if(clHeader)
	myfree(clHeader);
      clHeader = t;
    } else
      if((strcmp(opt, "cl_line_num") == 0)){
	if(clLineNum)
	  myfree(clLineNum);
	clLineNum = t;
      } else
	if((strcmp(opt, "cl_error") == 0)){
	  if(clError)
	    myfree(clError);
	  clError = t;
	} else
	  if((strcmp(opt, "cl_instructions") == 0)){
	    if(clInstr)
	      myfree(clInstr);
	    clInstr = t;
	  } else{
	    myfree(t);
	  }
}

void set_default_colors(){
  if(clError == NULL){
    clError = (char*)mymalloc(8 * sizeof(char));
    strcpy(clError, "\033[1;31m");
  }
  if(clInstr == NULL){
    clInstr = (char*)mymalloc(6 * sizeof(char));
    strcpy(clInstr, "\033[33m");
  }
  if(clHeader == NULL){
    clHeader = (char*)mymalloc(8 * sizeof(char));
    strcpy(clHeader, "\033[1;37m");
  }
  if(clLineNum == NULL){
    clLineNum = (char*)mymalloc(6 * sizeof(char));
    strcpy(clLineNum, "\033[35m");
  }
  if(clMenuSep == NULL){
    clMenuSep = (char*)mymalloc(6 * sizeof(char));
    strcpy(clMenuSep, "\033[34m");
  }
}

BOOLEAN read_options(){
  char *s, b[255], opt[128], value[128];
  FILE *F;
  int i, k;

  if(!ls_cmd){
#ifndef MSDOS
    ls_cmd = (char*)mymalloc(3 * sizeof(char));
    strcpy(ls_cmd, "ls");
#else
    ls_cmd = (char*)mymalloc(4 * sizeof(char));
    strcpy(ls_cmd, "DIR");
#endif
  }
  if(!(prefs_filename && myexist(prefs_filename)))
    return FALSE;

  FOPEN(prefs_filename, "r", F);
  while(fgets(b, 254, F)){
    s = b;
    while(s[0] == ' ' || s[0] == '\t')
      s++;
    if(strlen(s) < 4)
      continue;
    if(s[0] == '#')
      continue;
    memset(opt, 0, 128);
    memset(value, 0, 128);
    i = 0;
    while(i < 128 && ((s[0] >= 'a' && s[0] <= 'z') ||
	  (s[0] >= '0' && s[0] <= '9') || s[0] == '_')){
      opt[i] = s[0];
      i++;
      s++;
    }

    /* Don't remove \" from neither gnuplot_default_term nor gnuplot_png_font */
    if(strcmp(opt, "gnuplot_png_font") == 0 || 
	strcmp(opt, "gnuplot_default_term") == 0){
      while(s[0] == ' ' || s[0] == '\t' || s[0] == '=')
	s++;
    } else{
      while(s[0] == ' ' || s[0] == '\t' || s[0] == '=' || 
	  s[0] == '\"' || s[0] == '\'')
	s++;
    }
    i = 0;
    while(s[i] != '\n'){
      value[i] = s[i];
      i++;
    }
    i--;
    if(strcmp(opt, "gnuplot_png_font") == 0 || 
	strcmp(opt, "gnuplot_default_term") == 0)
      while(value[i] == ' ' || value[i] == '\t'){
	value[i] = 0;
	i--;
      }
    else
      while(value[i] == ' ' || value[i] == '\t' || value[i] == '\"'
	  || s[0] == '\''){
	value[i] = 0;
	i--;
      }

    if(strcmp(opt, "gnuplot_png_font") == 0){
      i = strlen(value) + 1;
      if(i > 3){
	if(gplt_png)
	  myfree(gplt_png);
	gplt_png = (char*)mymalloc(i * sizeof(char));
	strcpy(gplt_png, value);
      }
    }
    if(strcmp(opt, "gnuplot_default_term") == 0){
      i = strlen(value) + 1;
      if(i > 3){
	if(gplt_default_term)
	  myfree(gplt_default_term);
	gplt_default_term = (char*)mymalloc(i * sizeof(char));
	strcpy(gplt_default_term, value);
      }
    }
    if(strcmp(opt, "use_gnuplot") == 0)
      if(strcmp(value, "no") == 0)
	noplot = TRUE;
    if(strcmp(opt, "graphs_title_prefix") == 0){
      i = strlen(value) + 1;
      if(i > 1){
	if(gtprefix)
	  myfree(gtprefix);
	gtprefix = (char*)mymalloc(i * sizeof(char));
	strcpy(gtprefix, value);
      }
    }
    if(strcmp(opt, "color") == 0)
      if(strcmp(value, "yes") == 0)
	color = TRUE;
    if(opt[0] == 'c' && opt[1] == 'l' && opt[2] == '_')
      set_color(opt, value);
    if(strcmp(opt, "log") == 0)
      if(strcmp(value, "yes") == 0)
	log_set = TRUE;
    if(strcmp(opt, "bell") == 0)
      if(strcmp(value, "no") == 0)
	nobell = TRUE;
    if(strcmp(opt, "text_histogram") == 0)
      if(strcmp(value, "yes") == 0)
	thist = TRUE;
    if(strcmp(opt, "bernhard") == 0)
      if(strcmp(value, "yes") == 0)
	bernhard = TRUE;
    if(strcmp(opt, "verbose") == 0)
      if(strcmp(value, "no") == 0)
	verbose = FALSE;
    if(strcmp(opt, "max_results") == 0){
      k = atoi(value);
      if(k > 5 && k < 1000000)
	MRESULT = k;
    }
    if(strcmp(opt, "screen_lines") == 0){
      k = atoi(value);
      if(k > 10 && k < 100)
	rlines = k;
    }
    if(strcmp(opt, "screen_columns") == 0){
      k = atoi(value);
      if(k > 10 && k < 200)
	rcols = k;
    }
    if(strcmp(opt, "format_columns_out") == 0)
      if(strcmp(value, "no") == 0)
	format_columns_out = FALSE;
    if(strcmp(opt, "use_system_ls") == 0 && strcmp(value, "yes") == 0)
	system_ls = TRUE;
    if(strcmp(opt, "field_separator") == 0 && strlen(value) > 0)
      sep = value[0];
    if(strcmp(opt, "autodetect_header") == 0 && strcmp(value, "no") == 0)
	detect_header = FALSE;
    if(strcmp(opt, "int_as_int_in_ascii_files") == 0 && strcmp(value, "yes") == 0)
	int_as_int = TRUE;
    if(strcmp(opt, "ask_fileformat") == 0 && strcmp(value, "yes") == 0)
	ask_fileformat = TRUE;
    if(strcmp(opt, "na_string") == 0 && (strlen(value)) > 1){
      if(NODATA)
	myfree(NODATA);
      NODATA = (char*)mymalloc(sizeof(char) * (strlen(value) + 1));
      strcpy(NODATA, value);
    }
    if(strcmp(opt, "system_ls_command") == 0 && (strlen(value)) > 1){
      if(ls_cmd)
	myfree(ls_cmd);
      ls_cmd = (char*)mymalloc(sizeof(char) * (strlen(value) + 1));
      strcpy(ls_cmd, value);
    }
    if(strcmp(opt, "gnuplot_charset") == 0 && (strlen(value)) > 1){
      if(gnuplot_charset)
	myfree(gnuplot_charset);
      gnuplot_charset = (char*)mymalloc(sizeof(char) * (strlen(value) + 1));
      strcpy(gnuplot_charset, value);
    }
  }
  FCLOSE(F);
  return TRUE;
}

char * STRCAT(char *d, char *o, int *max){
  if((strlen(d) + strlen(o)) >= *max){
    *max += 1000;
    d = (char*)myrealloc(d, *max * sizeof(char));
  }
  strcat(d, o);
  return d;
}

void save_prefs(){
  char *s, b[255], *b2, value[255], opt[255];
  FILE *F;
  int max = 1000, i;
  BOOLEAN found = TRUE, o1 = TRUE, o2 = TRUE, o3 = TRUE, o4 = TRUE,
	  o5 = TRUE, o6 = TRUE, o7 = TRUE, o8 = TRUE, o9 = TRUE;

  if(prefs_filename == NULL){
    out_err(ERR, ERR_FILE, ERR_LINE,
	_("Could not build the name of the preferences file!"));
    return;
  }
  b2 = (char*)mymalloc(max * sizeof(char));
  b2[0] = 0;

  /* Reading the statistrc file, and replacing option values (in the buffer)*/
  if(myexist(prefs_filename)){
    FOPEN(prefs_filename, "r", F);
    while(fgets(b, 254, F)){
      s = b;
      while(s[0] == ' ' || s[0] == '\t')
	s++;
      if(strlen(s) < 4){
	b2 = STRCAT(b2, b, &max);
	continue;
      }
      if(s[0] == '#'){
	b2 = STRCAT(b2, b, &max);
	continue;
      }
      for(i = 0; i < 20; i++){
	opt[i] = 0;
	value[i] = 0;
      }
      i = 0;
      while((s[0] >= 'a' && s[0] <= 'z') || s[0] == '_'){
	opt[i] = s[0];
	i++;
	s++;
      }
      while(s[0] == ' ' || s[0] == '\t' || s[0] == '=')
	s++;
      i = 0;
      while((s[i] >= 'a' && s[i] <= 'z') || (s[i] >= '0' && s[i] <= '9')
	  || s[i] == '-' || s[i] == ' '){
	value[i] = s[i];
	i++;
      }
      i--;
      while(value[i] == ' ' || value[i] == '\t'){
	value[i] = 0;
	i--;
      }
      if(strcmp(opt, "verbose") == 0){
	found = TRUE;
	o1 = FALSE;
	if(verbose)
	  strcpy(value, "verbose = yes\n");
	else
	  strcpy(value, "verbose = no\n");
      }
      if(strcmp(opt, "use_gnuplot") == 0){
	found = TRUE;
	o2 = FALSE;
	if(noplot)
	  strcpy(value, "use_gnuplot = no\n");
	else
	  strcpy(value, "use_gnuplot = yes\n");
      }
      if(strcmp(opt, "bell") == 0){
	found = TRUE;
	o3 = FALSE;
	if(nobell)
	  strcpy(value, "bell = no\n");
	else
	  strcpy(value, "bell = yes\n");
      }
      if(strcmp(opt, "text_histogram") == 0){
	found = TRUE;
	o4 = FALSE;
	if(thist)
	  strcpy(value, "text_histogram = yes\n");
	else
	  strcpy(value, "text_histogram = no\n");
      }
      if(strcmp(opt, "bernhard") == 0){
	found = TRUE;
	o5 = FALSE;
	if(bernhard)
	  strcpy(value, "bernhard = yes\n");
	else
	  strcpy(value, "bernhard = no\n");
      }
      if(strcmp(opt, "use_system_ls") == 0){
	found = TRUE;
	o6 = FALSE;
	if(system_ls)
	  sprintf(value, "use_system_ls = yes\n");
	else
	  sprintf(value, "use_system_ls = no\n");
      }
      if(strcmp(opt, "max_results") == 0){
	found = TRUE;
	o7 = FALSE;
	sprintf(value, "max_results = %i\n", MRESULT);
      }
      if(strcmp(opt, "screen_lines") == 0){
	found = TRUE;
	o8 = FALSE;
	sprintf(value, "screen_lines = %i\n", rlines);
      }
      if(strcmp(opt, "screen_columns") == 0){
	found = TRUE;
	o9 = FALSE;
	sprintf(value, "screen_columns = %i\n", rcols);
      }
      if(found)
	b2 = STRCAT(b2, value, &max);
      else
	b2 = STRCAT(b2, b, &max);
      found = FALSE;
    }
    FCLOSE(F);
  }
  FOPEN(prefs_filename, "w", F);

  /* Saving preferences if they weren't found in statistrc yet, and if they
   * don't have the default values */
  if(o1 && !(verbose))
    STRCAT(b2, "verbose = no\n", &max);
  if(o2 && noplot)
    STRCAT(b2, "use_gnuplot = no\n", &max);
  if(o3 && nobell)
    STRCAT(b2, "bell = no\n", &max);
  if(o4 && thist)
    STRCAT(b2, "text_histogram = yes\n", &max);
  if(o5 && bernhard)
    STRCAT(b2, "bernhard = yes\n", &max);
  if(o6 && system_ls)
    STRCAT(b2, "use_system_ls = yes\n", &max);
  if(o7 && MRESULT != 200){
    sprintf(value, "max_results = %i\n", MRESULT);
    STRCAT(b2, value, &max);
  }
#ifdef MSDOS
  i = 24;
#else
  i = 40;
#endif
  if(o8 && rlines){
    sprintf(value, "screen_lines = %i\n", rlines);
    STRCAT(b2, value, &max);
  }
  if(o9 && rcols){
    sprintf(value, "screen_columns = %i\n", rcols);
    STRCAT(b2, value, &max);
  }
  FWRITE(b2, 1, strlen(b2), F);
  FCLOSE(F);

  out_d(_("File \"%s\" saved!"), prefs_filename);
  out_d("\n\n");
  mywait();
  myfree(b2);
}

void set_winsize(){
  SCRLINES = 0;
  SCRCOLS = 0;
#ifndef NO_IOCTIL_H
  struct winsize ws;
  if(ioctl(1, TIOCGWINSZ, &ws) != -1) {
    SCRCOLS = ws.ws_col;
    SCRLINES = ws.ws_row;
  }
#endif
  if(rcols)
    SCRCOLS = rcols;
  if(rlines)
    SCRLINES = rlines;
  if(SCRLINES < 10 || SCRLINES > 400){ /* xterm with tiny font: 1020x359 */
#ifdef MSDOS
    SCRLINES = 24;
#else
    SCRLINES = 40;
#endif
  }
  if(SCRCOLS < 10 || SCRCOLS > 2000)
    SCRCOLS = 80;
}

void colorize(int cl){
  if(!color)
    return;
#ifdef MSDOS
  return;
#endif
  switch(cl){
    case ClDefault :
      out_d("\033[0m");
      break;
    case ClError :
      fprintf(stderr, "%s", clError);
      break;
    case ClErrorD :
      fprintf(stderr, "\033[0m");
      break;
    case ClInstr :
      out_d(clInstr);
      break;
    case ClHeader :
      out_d(clHeader);
      break;
    case ClLineNum :
      out_d(clLineNum);
      break;
    case ClMenuSep :
      out_d(clMenuSep);
      break;
  }
}

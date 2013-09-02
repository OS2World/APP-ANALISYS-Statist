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
**  $Id: plot.c,v 1.20 2006/09/30 21:10:38 jakson Exp $
***************************************************************/
/* plot.c for STATIST */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#ifndef MSDOS
#include <iconv.h>
#endif

#include "statist.h"
#include "funcs.h"
#include "data.h"
#include "plot.h"
#include "gettext.h"

#ifndef NOPIPE
#define PLOT_MSG   out_d(_("Creating gnuplot-graphic ...\n\n"))
#define PLOT_CLOSE fflush(pipef)
#else
#define PLOT_MSG   out_d(_("Creating gnuplot command file '%s'\n\n"), GPL_COM)
#define PLOT_CLOSE FCLOSE(pipef)
#endif

char GPL_DAT[300];
char GPL_COM[16];
char **plothist;
int histsize = 30;
int lastcmd = -1;
BOOLEAN has_graph = FALSE;

char *gnuplot_charset;

/* Charset conversion function */
#ifndef MSDOS
char * g(char *s){
  if(gnuplot_charset == NULL || is_utf8 == FALSE || 
      strcmp(gnuplot_charset, "UTF-8") == 0)
    return s;
  char *b1, *b2;
  b1 = s;
  b2 = s;
  int error;
  char *inbuf, *outbuf;
  size_t t, inbytesleft, outbytesleft;
  iconv_t cd = iconv_open(gnuplot_charset, "UTF-8");
  if(cd == (iconv_t)-1){
    if(errno == EINVAL)
      out_err(ERR, ERR_FILE, ERR_LINE,
	  _("Invalid encoding specification: \"%s\""), gnuplot_charset);
    else
      out_err(ERR, ERR_FILE, ERR_LINE, "iconv_open() error");
  }
  inbuf = b1;
  inbytesleft = strlen(b1);
  outbytesleft = inbytesleft;
  outbuf = (char*)m_calloc((outbytesleft + 1), sizeof(char));
  b2 = outbuf;
  t = iconv(cd, &inbuf, &inbytesleft, &outbuf, &outbytesleft);
  error = errno;
  if (iconv_close(cd) != 0) {
    out_err(ERR, ERR_FILE, ERR_LINE, "iconv_close() error");
  }
  if(t == -1){
    char *errmsg = (char*)m_calloc((strlen(s) + 255), sizeof(char));
    char where[80];
    int i = 0, l;
    l = strlen(b2);
    switch(error){
      case EILSEQ :
	while(i < 77 && i < l){
	  where[i] = ' ';
	  i++;
	}
	where[i] = '^';
	where[i+1] = 0;
	sprintf(errmsg, "%s\n  %s  %s", _("An invalid multibyte sequence was "
	      "encountered in the input:"), s, where);
	out_err(WAR, ERR_FILE, ERR_LINE, errmsg);
	break;
      case EINVAL :
	out_err(WAR, ERR_FILE, ERR_LINE, 
	    "An  incomplete multibyte sequence was encountered in the "
	    "input, and the input byte sequence terminates after it");
	break;
    }
    b2 = b1;
  }
  return b2;
}
#else
char * g(char *s){
  return(s);
}
#endif

char * fixquote(char *s){
  char * fixedlabel;
  int i = 0, j = 0;
  int l = strlen(s);
  fixedlabel = (char*)m_calloc((2 * l), sizeof(char));
  while(i < l){
    if(s[i] == '\"' && (i == 0 || s[i-1] != '\\')){
      fixedlabel[j] = '\\';
      fixedlabel[j+1] = '\"';
      j += 1;
    } else{
      fixedlabel[j] = s[i];
    }
    i++;
    j++;
  }
  fixedlabel[j] = 0;
  return g(fixedlabel);
}

void gnuplotcmd(char *fmt, ...){
  va_list argptr;
  char s[255];
  va_start(argptr, fmt);
  vsprintf(s, fmt, argptr);
  fprintf(pipef, "%s\n", s);
  if(histsize > 0){
    lastcmd++;
    if(lastcmd == histsize)
      lastcmd = 0;
    plothist[lastcmd] = (char*)mycalloc((strlen(s)+1), sizeof(char));
    strcpy(plothist[lastcmd], s);
  }
  va_end(argptr);
}

void clear_history(){
  int i;
  for(i = 0; i < histsize; i++)
    if(plothist[i] != NULL){
      myfree(plothist[i]);
      plothist[i] = NULL;
    }
}

void set_default(){
  clear_history();
  gnuplotcmd("unset key");
  gnuplotcmd("unset parametric");
  gnuplotcmd("unset log x");
  gnuplotcmd("unset grid");
  gnuplotcmd("unset label");
  gnuplotcmd("unset ylabel");
  gnuplotcmd("set xtics autofreq");
  gnuplotcmd("set ytics autofreq");
  gnuplotcmd("linetype=1");
  has_graph = TRUE;
}

BOOLEAN init_gnuplot() {
  int i = 1;
  if(plothist == NULL)
    plothist = (char**)mycalloc(histsize, sizeof(char*));

#ifndef NOPIPE
  switch (gnupl_open) {
    case CRASH:
      return FALSE;
      break;
    case TRUE:
      return TRUE;
      break;
    case FALSE:
#ifdef MSDOS
      if ((pipef = popen("pgnuplot --version", "w"))==NULL) {
#else
      if ((pipef = popen("gnuplot --version", "w"))==NULL) {
#endif
	out_err(ERR, ERR_FILE, ERR_LINE, _("gnuplot-initialization failed!"));
	gnupl_open = CRASH;
	return FALSE;
      }
      else {
	/* popen() has changed. Now we have to close the pipe to know whether
	 * gnuplot is in the path. Adopting plotdrop solution: */
	int status = pclose(pipef);
	pipef = NULL;
	if(status){
	  out_err(ERR, ERR_FILE, ERR_LINE, _("gnuplot-initialization failed!"));
	  gnupl_open = CRASH;
	  pipef = NULL;
	  return FALSE;
	} else{
#ifdef MSDOS
	  pipef = popen("pgnuplot", "w");
#else
	  pipef = popen("gnuplot -geometry 450x300", "w");
#endif
	  if(pipef){
	    gnupl_open = TRUE;
#ifdef MSDOS
	    sprintf(GPL_DAT, "stat_gpl.dat");
#else
	    sprintf(GPL_DAT, "%.256s/stat_gpl.dat", getenv("HOME"));
	    while(myexist(GPL_DAT)){
	      i++;
	      sprintf(GPL_DAT, "%.256s/stat_gpl_%i.dat", getenv("HOME"), i);
	    }
	    if(gplt_default_term)
	      gnuplotcmd("set term %s", gplt_default_term);
#endif
	    return TRUE;
	  }
	  return FALSE;
	}
      }
    default:
      return FALSE;
  }

#else
  strcpy(GPL_COM, "stat_gpl.plt");
  FOPEN(GPL_COM, "w", pipef);
  gnupl_open = TRUE;
  sprintf(GPL_DAT, "stat_gpl.dat");
  return TRUE;
#endif
}

void save_png(){
  if(gnupl_open != TRUE)
    return;
  char fn[MLINE];
  out_i(_("Please enter a name for the png file: ") );
  GETRLINE;
  sscanf(line, "%s", fn);

  /* Add .png filename extension if necessary */
  int i = strlen(fn);
  if(!(fn[i-4] == '.' && fn[i-3] == 'p' && fn[i-2] == 'n' && fn[i-1] == 'g'))
    strncat(fn, ".png", MLINE - i - 5);

  /* I added out_d() lines because I suppose that users expect some feedback,
   * but I can't guarantee that the picture will be created. */
  out_d("\n");
  if(gplt_png){
    gnuplotcmd("set term png font %s", gplt_png);
    out_d("gnuplot> set term png font %s\n", gplt_png);
  } else{
    gnuplotcmd("set term png");
    out_d("gnuplot> set term png\n");
  }
  gnuplotcmd("set output \"%s\"", fn);
  out_d("gnuplot> set output \"%s\"\n", fn);
  gnuplotcmd("replot");
  out_d("gnuplot> replot\n");
  out_d("\n");
  if(gplt_default_term)
    gnuplotcmd("set term %s", gplt_default_term);
  else
    gnuplotcmd("set term x11");
  has_graph = FALSE;
}

BOOLEAN plot_pair(REAL x[], REAL y[], int n, REAL a0, REAL a1,
		  char *xlab, char *ylab) {
  FILE *tempf;
  int i;
  char aline[160];

  if (!init_gnuplot()) {
    return FALSE;
  }

#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  FOPEN(GPL_DAT, "wt", tempf);
  for (i=0; i<n; i++) {
    fprintf(tempf, "%g %g\n", x[i], y[i]);
  }
  FCLOSE(tempf);
  set_default();
  sprintf(aline, "f(x)=%g + %g*x", a0, a1);
  gnuplotcmd("%s", aline);
  gnuplotcmd("set xlabel \"%s\"", fixquote(xlab));
  gnuplotcmd("set ylabel \"%s\"", fixquote(ylab));
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  if(gtprefix)
    sprintf(aline, "%s%s", fixquote(gtprefix), g(_("Linear Regression")));
  else
    sprintf(aline, "%s", g(_("Linear Regression")));
  gnuplotcmd("set title \"%s\"", aline);
  gnuplotcmd("plot '%s', f(x)", GPL_DAT);
  PLOT_CLOSE;
  PLOT_MSG;
  return TRUE;
}


BOOLEAN plot_tripl(REAL x[], REAL y[], REAL z[], int n,
                   REAL a0, REAL a1, REAL a2,
		   char *xlab, char *ylab, char *zlab) {
  FILE *tempf;
  int i;
  char aline[160];

  if (!init_gnuplot()) {
    return FALSE;
  }

#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  FOPEN(GPL_DAT, "wt", tempf);
  for (i=0; i<n; i++) {
    fprintf(tempf, "%g %g %g\n", x[i], y[i], z[i]);
  }
  FCLOSE(tempf);
  sprintf(aline,  "f(u,v)=%g + %g*u + %g*v\n", a0, a1, a2);
  set_default();
  gnuplotcmd("set parametric");
  gnuplotcmd("%s", aline);
  gnuplotcmd("set xlabel \"%s\"", fixquote(xlab));
  gnuplotcmd("set ylabel \"%s\"", fixquote(ylab));
  gnuplotcmd("set zlabel \"%s\"", fixquote(zlab));
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  if(gtprefix)
    sprintf(aline, "%s%s", fixquote(gtprefix), g(_("Multiple Linear Regression")));
  else
    sprintf(aline, "%s", g(_("Multiple Linear Regression")));
  gnuplotcmd("set title \"%s\"", aline);
#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  gnuplotcmd("splot [%g:%g][%g:%g][] '%s', '%s' with impulses, u,v,f(u,v)",
      get_min(x,n), get_max(x,n),get_min(y,n),get_max(y,n), GPL_DAT, GPL_DAT);
  PLOT_CLOSE;
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  PLOT_MSG;
  return TRUE;
}


BOOLEAN plot_cdf_ks(REAL z[], int n, REAL z0, REAL fn1, REAL fn2,
		 REAL x[], char *xlab) {
  const REAL zmin=(-3.0), zmax=3.0;
  FILE *tempf;
  REAL cdf;
  int i;
  char aline[160];

  if (!init_gnuplot()) {
    return FALSE;
  }

#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  FOPEN(GPL_DAT, "wt", tempf);
  fprintf(tempf, "%g %g\n", zmin, 0.0);
  fprintf(tempf, "%g %g\n", z[0], 0.0);
  for (i=0; i<(n-1); i++) {
    cdf = (REAL)(i+1)/(REAL)n;
    fprintf(tempf, "%g %g\n", z[i], cdf);
    fprintf(tempf, "%g %g\n", z[i+1], cdf);
  }
  fprintf(tempf, "%g %g\n", z[n-1], 1.0);
  fprintf(tempf, "%g %g\n", zmax, 1.0);
/*
  fprintf(tempf, "\n%g %g\n%g %g\n", zmin, fn1, zmax, fn1);
  fprintf(tempf, "\n%g %g\n%g %g\n", zmin, fn2, zmax, fn2);
*/
  FCLOSE(tempf);

  set_default();
  gnuplotcmd("set xlabel \"N(%s)\"", fixquote(xlab));
  gnuplotcmd("set ylabel \"S(%s)\"", fixquote(xlab));
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  if(gtprefix)
    sprintf(aline, "%s%s", fixquote(gtprefix), g(_("KS-Lilliefors-Test")));
  else
    sprintf(aline, "%s", g(_("KS-Lilliefors-Test")));
  gnuplotcmd("set title \"%s\"", aline);
#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  gnuplotcmd("u(x)=%g", fn1);
  gnuplotcmd("l(x)=%g", fn2);
  gnuplotcmd("plot [%g:%g] [0:1.02] '%s' with lines, norm(x), u(x), l(x) with lines 3",
	  zmin, zmax, GPL_DAT);
  PLOT_CLOSE;
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  PLOT_MSG;
  return TRUE;
}

BOOLEAN plot_cdf(REAL z[], int n, char *xlab) {
  const REAL zmin= z[0], zmax = z[n-1];
  FILE *tempf;
  REAL cdf;
  int i;
  char aline[160];

  if (!init_gnuplot()) {
    return FALSE;
  }

#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  FOPEN(GPL_DAT, "wt", tempf);
/*  fprintf(tempf, "%g %g\n", z[0], 0.0);  */
  for (i=0; i<(n-1); i++) {
    cdf = (REAL)(i+1)/(REAL)n;
    fprintf(tempf, "%g %g\n", z[i], cdf);
    fprintf(tempf, "%g %g\n", z[i+1], cdf);
  }
  fprintf(tempf, "%g %g\n", zmax, 1.0);
  FCLOSE(tempf);

  set_default();
  gnuplotcmd("set xlabel \"%s\"", fixquote(xlab));
  gnuplotcmd("set ylabel \"%%(%s)\"", fixquote(xlab));
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  if(gtprefix)
    sprintf(aline, "%s%s", fixquote(gtprefix), g(_("Cumulative distribution")));
  else
    sprintf(aline, "%s", g(_("Cumulative distribution")));
  gnuplotcmd("set title \"%s\"", aline);
#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  gnuplotcmd("plot [%g:%g] [0:1.02] '%s' with lines", zmin, zmax, GPL_DAT);
  PLOT_CLOSE;
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  PLOT_MSG;
  return TRUE;
}


BOOLEAN plot_histo(REAL x[], int y[], int n, REAL step, REAL data[],
		   char *datalab) {
  FILE *tempf;
  unsigned int i, miny, maxy;
  REAL minx, maxx;
  char aline[160];

  if (!init_gnuplot()) {
    return FALSE;
  }

#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  FOPEN(GPL_DAT, "wt", tempf);
  for (i=0; i<n; i++) {
    fprintf(tempf, "%g %i\n", x[i], y[i]);
  }
  miny = 0;
  maxy = (int)(get_maxint(y, n)*1.05) + 1;
  minx = get_min(x, n);
  maxx = get_max(x, n);
  FCLOSE(tempf);
  set_default();
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  if(gtprefix)
    sprintf(aline, "%s%s", fixquote(gtprefix), g(_("Histogram")));
  else
    sprintf(aline, "%s", g(_("Histogram")));
  gnuplotcmd("set title \"%s\"", aline);
  gnuplotcmd("set ylabel '%s'", g(_("Frequency")));
#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  gnuplotcmd("set xlabel \"%s\"", fixquote(datalab));
  if(((maxx+step/2.)-(minx-step/2.)) /step >10)
  {
  	gnuplotcmd("set xtics autofreq");
  } else
  	gnuplotcmd("set xtics %g, %g, %g",
		(minx-step/2.), step, (maxx+step/2.));
  gnuplotcmd("plot [%g:%g][%i:%i] '%s' with boxes linetype linetype",
          minx-step/2., maxx+step/2., miny, maxy, GPL_DAT);
  PLOT_CLOSE;
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  PLOT_MSG;
  return TRUE;
}


BOOLEAN plot_histo_discrete(REAL x[], int y[], int n, REAL step, REAL data[],
		   char *datalab) {
  FILE *tempf;
  unsigned int i, miny, maxy;
  REAL minx, maxx;
  char aline[160];

  if (!init_gnuplot()) {
    return FALSE;
  }

#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  FOPEN(GPL_DAT, "wt", tempf);
  for (i=0; i<n; i++) {
    fprintf(tempf, "%g %i\n", x[i], y[i]);
  }
  miny = 0;
  maxy = get_maxint(y, n) + 1;
  minx = get_min(x, n) - (step/2.);
  maxx = get_max(x, n) + (step/2.);
  FCLOSE(tempf);
  set_default();
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  if(gtprefix)
    sprintf(aline, "%s%s", fixquote(gtprefix), g(_("Histogram")));
  else
    sprintf(aline, "%s", g(_("Histogram")));
  gnuplotcmd("set title \"%s\"", aline);
  gnuplotcmd("set ylabel '%s'", g(_("Frequency")));
  gnuplotcmd("set xlabel \"%s\"", fixquote(datalab));
#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  gnuplotcmd("set xtics %g, %g, %g", (minx+step/2.), step, (maxx-step/2.));
/*  gnuplotcmd("set boxwidth %g\n", step/4. ); */
/*  gnuplotcmd("plot [%g:%g][%i:%i] '%s' with boxes\n", */
  gnuplotcmd("plot [%g:%g][%i:%i] '%s' with impulses",
          minx, maxx, miny, maxy, GPL_DAT);
  PLOT_CLOSE;
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  PLOT_MSG;
  return TRUE;
}

BOOLEAN plot_probit(REAL dose[], REAL num[], REAL effect[], int n,
                    REAL a0, REAL a1, REAL dose0, REAL dose1,
		    char *doselab, char *effectlab) {
  FILE *tempf;
  int i;
  char aline[160];
#ifndef STATIST_X
  out_r("plot_probit: doselab=|%s|, effectlab=|%s|\n",
	 doselab, effectlab);
#endif
  if (!init_gnuplot()) {
    return FALSE;
  }

#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  FOPEN(GPL_DAT, "wt", tempf);
  for (i=0; i<n; i++) {
    fprintf(tempf, "%g %g\n", dose[i], (effect[i]/num[i])*100.);
  }
  FCLOSE(tempf);
  set_default();
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  sprintf(aline, g(_("Dose %s")), fixquote(doselab));
  gnuplotcmd("set xlabel \"%s\"", aline);
  sprintf(aline, g(_("Effect %s [%%]")), fixquote(effectlab));
  gnuplotcmd("set ylabel \"%s\"", aline);
  if(gtprefix)
    sprintf(aline, "%s%s", fixquote(gtprefix), g(_("Probit analysis")));
  else
    sprintf(aline, "%s", g(_("Probit analysis")));
  gnuplotcmd("set title \"%s\"", aline);
#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  gnuplotcmd("set grid");
  gnuplotcmd("set log x");
  gnuplotcmd("a0=%g", a0);
  gnuplotcmd("a1=%g", a1);
  gnuplotcmd("f(x)=norm((log10(x)*a1+a0)-5)*100");
  gnuplotcmd("plot [%g:%g] '%s', f(x)",
          dose0, dose1, GPL_DAT);
  PLOT_CLOSE;
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  PLOT_MSG;
  return TRUE;
}



BOOLEAN plot_box(REAL x[], int n, REAL median, REAL mean, REAL q_l,
                 REAL q_u, REAL w_l, REAL w_u, REAL no_l,
		 REAL no_u, char *xlab) {
  FILE *tempf;
  int i;
  REAL min, max, off;
  char aline[160];

  if (!init_gnuplot()) {
    return FALSE;
  }

  min = get_min(x, n);
  max = get_max(x, n);
  off = (max-min)/20.;

#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  FOPEN(GPL_DAT, "wt", tempf);
  fprintf(tempf, "%g 10\n", q_l);
  fprintf(tempf, "%g 12\n", q_l);
  fprintf(tempf, "%g 12\n", q_u);
  fprintf(tempf, "%g 10\n", q_u);
  fprintf(tempf, "%g 10\n\n", q_l);
  fprintf(tempf, "%g 10\n", median);
  fprintf(tempf, "%g 12\n\n", median);
  fprintf(tempf, "%g 11\n", w_l);
  fprintf(tempf, "%g 11\n\n", q_l);
  fprintf(tempf, "%g 10.7\n", w_l);
  fprintf(tempf, "%g 11.3\n\n", w_l);
  fprintf(tempf, "%g 11\n", q_u);
  fprintf(tempf, "%g 11\n\n", w_u);
  fprintf(tempf, "%g 10.7\n", w_u);
  fprintf(tempf, "%g 11.3\n\n", w_u);
  fprintf(tempf, "%g 11.1\n", no_l+off*0.2);
  fprintf(tempf, "%g 11.1\n", no_l);
  fprintf(tempf, "%g 10.9\n", no_l);
  fprintf(tempf, "%g 10.9\n\n", no_l+off*0.2);
  fprintf(tempf, "%g 11.1\n", no_u-off*0.2);
  fprintf(tempf, "%g 11.1\n", no_u);
  fprintf(tempf, "%g 10.9\n", no_u);
  fprintf(tempf, "%g 10.9\n\n", no_u-off*0.2);
  FCLOSE(tempf);

  set_default();

  for (i=0; i<n; i++) {
    if ((x[i]>w_u) || (x[i]<w_l)) {
      gnuplotcmd("set label '*' at %g, 11 center", x[i]);
    }
  }

  gnuplotcmd("unset ytics");
  gnuplotcmd("set label");
  gnuplotcmd("set label 'o' at %g, 11 center", mean);
  gnuplotcmd("set label 'n=%i' at %g, 12.7 right", n, max);
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  if(gtprefix)
    sprintf(aline, "%s%s", fixquote(gtprefix), g(_("Box-and-Whisker Plot")));
  else
    sprintf(aline, "%s", g(_("Box-and-Whisker Plot")));
  gnuplotcmd("set title \"%s\"", aline);
#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  gnuplotcmd("set xlabel \"%s\"", fixquote(xlab));
  gnuplotcmd("plot [%g:%g][9:13] '%s' with lines",
          min-off, max+off, GPL_DAT);
  PLOT_CLOSE;
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  PLOT_MSG;
  return TRUE;
}




BOOLEAN  plot_poly(REAL x[], REAL y[], int n, REAL a[], int npoly,
		   char *xlab, char *ylab) {
  FILE *tempf;
  int i;
  char term[255], aline[255], grtitle[255];

  if (!init_gnuplot()) {
    return FALSE;
  }

#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  FOPEN(GPL_DAT, "wt", tempf);
  for (i=0; i<n; i++) {
    fprintf(tempf, "%g %g\n", x[i], y[i]);
  }
  FCLOSE(tempf);
  set_default();
  strcpy(aline, "f(x) = ");
  for (i=0; i<=npoly; i++) {
    gnuplotcmd("a%i=%g", i, a[i]);
    sprintf(term, "a%i*x**%i+", i, i);
    strncat(aline, term, 255-strlen(aline));
  }
  strcat(aline, "0\n");
  gnuplotcmd("%s", aline);
  gnuplotcmd("set xlabel \"%s\"", fixquote(xlab));
  gnuplotcmd("set ylabel \"%s\"", fixquote(ylab));
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  sprintf(grtitle, g(_("Polynomial Regression of Order %i")), npoly);
  if(gtprefix)
    sprintf(aline, "%s%s", fixquote(gtprefix), grtitle);
  else
    sprintf(aline, "%s", grtitle);
  gnuplotcmd("set title \"%s\"", aline);
  gnuplotcmd("plot '%s', f(x)", GPL_DAT);
  PLOT_CLOSE;
  PLOT_MSG;
  return TRUE;
}


#ifndef MSDOS
BOOLEAN plot_command() {
  char aline[80];
  char *s;
  int i, j;

  if (!init_gnuplot()) {
    return FALSE;
  }

  if(plothist && lastcmd > -1){
    out_d(_("Last commands sent to gnuplot:\n\n"));
    j = lastcmd;
    for(i = 0; i < histsize; i++){
      j++;
      if(j == histsize)
	j = 0;
      if(plothist[j])
	out_d("  %s\n", plothist[j]);
    }
    out_d("\n");
  }
  out_i(_("Terminate gnuplot commands with '.'\n\n") );
  
  do {
    out_i("gnuplot> ");
    fgets(aline, 79, stdin);
    s = aline;
    while(s[0] == ' ' || s[0] == '\t')
      s++;
    if (s[0] == '.')
      break;
    if (s[0] == 'q' || (s[0] == 'e' && s[1] == 'x')){
      pclose(pipef);
      gnupl_open = FALSE;
      has_graph = FALSE;
      clear_history();
      if(myexist(GPL_DAT))
	remove(GPL_DAT);
      break;
    }
    gnuplotcmd(g(s));
    fflush(pipef);
  } while (1);
  return TRUE;
}
#endif

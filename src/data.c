/* This file is part of statist
**
** It is distributed under the GNU General Public License.
** See the file COPYING for details.
**
** (c) 1997 Dirk Melcher
** old email address: Dirk.Melcher@usf.Uni-Osnabrueck.DE
**
** adapted for statistX: Andreas Beyer, 1999, abeyer@usf.uni-osnabrueck.de
**
** The function get_line() was adapted from GNU coretutils 5.2.1 getndelim2.c
** Copyright (C) 1993, 1996, 1997, 1998, 2000, 2003 Free Software
** Foundation, Inc.
**
**  published by Bernhard Reiter  http://www.usf.Uni-Osnabrueck.DE/~breiter
**  $Id: data.c,v 1.44 2006/10/21 23:15:22 jakson Exp $
***************************************************************/

/* data.c for STATIST */
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "statist.h"
#include "data.h"
#include "funcs.h"
#include "menue.h"

#include "gettext.h"

/* ==================================================================== */

static short int *labelcol;
static short int n_lab;


/* ==================================================================== */

void inflate_MCOL(){
  int i, oldmax = MCOL;
  MCOL += 64;
  xx = (PREAL*)myrealloc(xx, (MCOL * sizeof(PREAL)));
  nn = (int*)myrealloc(nn, (MCOL * sizeof(int)));
  vn = (int*)myrealloc(vn, (MCOL * sizeof(int)));
  acol = (int*)myrealloc(acol, (MCOL * sizeof(int)));
  x_read = (BOOLEAN*)myrealloc(x_read, (MCOL * sizeof(BOOLEAN)));
  alias = (char**)myrealloc(alias, (MCOL * sizeof(char*)));
  tmpptr = (FILE**)myrealloc(tmpptr, (MCOL * sizeof(FILE*)));
  labelcol = (short int*)myrealloc(labelcol, (MCOL * sizeof(short int)));
  names = (Labels**)myrealloc(names, (MCOL * sizeof(Labels*)));

   for (i = oldmax; i < MCOL; i++){
     xx[i] = NULL;
     nn[i] = 0;
     vn[i] = 0;
     acol[i] = 0;
     x_read[i] = FALSE;
     alias[i] = NULL;
     tmpptr[i] = NULL;
     alias[i] = NULL;
     names[i] = NULL;
   }
}

char * get_default_label(int i){
  char * newlabel;
  char * tempLabel = (char*)mycalloc(4, sizeof(char));
  if((i+CHAR_OFFSET) <= 'z'){
    tempLabel[0] = (char)(i+CHAR_OFFSET);
  } else{
    tempLabel[0] = (char)(((i + 26) / 26) + CHAR_OFFSET - 2);
    tempLabel[1] = (char)((i  % 26) + CHAR_OFFSET);
  }
  newlabel = (char*)mymalloc(sizeof(char) * (strlen(tempLabel) + 1));
  strcpy(newlabel, tempLabel);
  myfree(tempLabel);
  return(newlabel);
}

void create_columns(int amount){
  int i;
  for(i = 0; i < amount; i++){
    if(ncol == MCOL)
      inflate_MCOL();
    if(alias[ncol] == NULL)
      alias[ncol] = get_default_label(ncol);
    tmpptr[ncol] = tmpfile();
    if(tmpptr[ncol] == NULL){
      out_err(FAT, ERR_FILE, ERR_LINE, 
	  _("System reports error while opening temporary file:\n  \"%s\""),
	  STRERROR(errno));
    }
    ncol++;
  }
}

/* Free allocated memory, but not the temporary file */
void free_column(int i){
  if((x_read[i])){
    myfree(xx[i]);
    xx[i] = NULL;
    x_read[i] = FALSE;
    vn[i] = 0;
  }
}

/* Free allocated memory and erase temporary file */
void delete_column(int i){
  free_column(i);
  if(alias[i]) {
    myfree(alias[i]);
    alias[i] = NULL;
  }
  names[i] = NULL;
  if(tmpptr[i]){
    FCLOSE(tmpptr[i]);
    tmpptr[i] = NULL;
  }
  nn[i] = 0;
  x_read[i] = FALSE;
  labelcol[i] = 0;
  ncol--;
}

void erasetempfiles() {
   int i;
   if(MCOL == 0)
     return;
   out_d(_("Removing temporary files ...") );
   for(i = 0; i < MCOL; i++)
     delete_column(i);
   myfree(xx);
   myfree(alias);
   myfree(nn);
   myfree(vn);
   myfree(acol);
   myfree(x_read);
   myfree(tmpptr);
   myfree(labelcol);
   myfree(names);
   xx = NULL;
   alias = NULL;
   nn = NULL;
   vn = NULL;
   acol = NULL;
   x_read = NULL;
   tmpptr = NULL;
   labelcol = NULL;
   names = NULL;
   out_d(_(" done\n") );
   n_lab = 0;
   ncol = 0;
   MCOL = 0;
 }


/* ==================================================================== */


/* Adapted from GNU coretutils 5.2.1 getndelim2.c */
int get_line(char **lineptr, size_t *linesize, FILE *stream){
  register int c;
  int pos = -1; /* index of last byte read */
  char * line = *lineptr;
  size_t max = *linesize - 1;
  for(;;){
    c = getc (stream);
    if(c == EOF){
      /* Return partial line, if any.  */
      if (pos == -1)
	return -1;
      else
	break;
    }
    if(pos == max){
      max += 64;
      *linesize += 64;
      *lineptr = myrealloc(*lineptr, *linesize);
      line = *lineptr;
    }
    pos++;
    line[pos] = c;
    if (c == '\n')
      /* Return the line.  */
      break;
  }
  pos++;
  line[pos] = '\0';
  pos++;
  return pos;
}

/* =================================================================== */

void attach_labels_to_columns(){
  int i;
  Labels *ptr = first_labels;
  for(i = 0; i < ncol; i++)
    names[i] = NULL;
  while(ptr){
    for(i = 0; i < ncol; i++){
      if(strcmp(alias[i], ptr->clabel) == 0){
	names[i] = ptr;
	break;
      }
    }
    ptr = ptr->next;
  }
}

void delete_labels(Labels *ptr){
  int i;
  Labels *p;
  if(ptr == first_labels){
    first_labels = first_labels->next;
  } else{
    p = first_labels;
    while(p->next != ptr)
      p = p->next;
    p->next = ptr->next;
  }
  if(ptr->clabel)
    myfree(ptr->clabel);
  if(ptr->ctitle)
    myfree(ptr->ctitle);
  if(ptr->n > 0){
    for(i = 0; i < ptr->n; i++)
      myfree(ptr->l[i]);
  }
  if(ptr->v)
    myfree(ptr->v);
  if(ptr->l)
    myfree(ptr->l);
  myfree(ptr);
}

/* Delete problematic Labels */
void check_labels(){
  Labels *next, *ptr = first_labels;
  while(ptr){
    next = ptr->next;
    if(ptr->clabel == NULL || (ptr->ctitle == NULL && ptr->n == 0))
      delete_labels(ptr);
    ptr = next;
  }
}

/* Creates a linked list of "Labels". The list of Labels might have labels for
 * columns that don't exist in the current datafile. */
void read_labels(char *labelsfile){
  BOOLEAN getting_labels = FALSE;
  Labels *ptr = NULL;
  FILE *F;
  int i, max = 0;
  char b[255], *s, t[255];
  FOPEN(labelsfile, "r", F);
  while(fgets(b, 254, F)){
    s = b;
    while(s[0] == ' ' || s[0] == '\t')
      s++;
    if(strlen(s) < 2){
      if(getting_labels){
	getting_labels = FALSE;
      }
      continue;
    }
    if(s[0] == '#')
      continue;
    i = 0;
    if(!getting_labels){
      if(ptr == NULL){
	first_labels = (Labels*)mycalloc(1, sizeof(Labels));
	ptr = first_labels;
      } else{
	ptr->next = (Labels*)mycalloc(1, sizeof(Labels));
	ptr = ptr->next;
      }
      while(!(s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == '\r')){
	t[i] = s[i];
	i++;
      }
      t[i] = 0;
      ptr->clabel = (char*)mymalloc((strlen(t) + 1) * sizeof(char));
      strcpy(ptr->clabel, t);
      s += i;
      while(s[0] == ' ' || s[0] == '\t')
	s++;
      i = 0;
      while(!(s[i] == '\n' || s[i] == '\r')){
	t[i] = s[i];
	i++;
      }
      t[i] = 0;
      if(strlen(t) > 1){
	ptr->ctitle = (char*)mymalloc((strlen(t) + 1) * sizeof(char));
	strcpy(ptr->ctitle, t);
      }
      max = 0;
      ptr->n = 0;
      getting_labels = TRUE;
    } else{
      if(ptr->n == max){
	max += 10;
	ptr->v = (REAL*)myrealloc(ptr->v, (max * sizeof(REAL)));
	ptr->l = (char**)myrealloc(ptr->l, (max * sizeof(char*)));
      }
      while(!(s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == '\r')){
	t[i] = s[i];
	i++;
      }
      t[i] = 0;
      if(sscanf(t, "%lf", &(ptr->v[ptr->n])) == 1){
	s += i;
	while(s[0] == ' ' || s[0] == '\t')
	  s++;
	i = 0;
	while(!(s[i] == '\n' || s[i] == '\r')){
	  t[i] = s[i];
	  i++;
	}
	t[i] = 0;
	i = strlen(t);
	if(i > 1){
	  ptr->l[ptr->n] = (char*)mymalloc((i+1) * sizeof(char));
	  strcpy(ptr->l[ptr->n], t);
	  ptr->n++;
	} else
	  continue;
      } else
	continue;
    }
  }
  FCLOSE(F);
  check_labels();
  if(first_labels == NULL)
    out_err(ERR, ERR_FILE, ERR_LINE,
	_("No labels found in \"%s\"!"), labelsfile);
  else
    attach_labels_to_columns();
}

void set_fileformat(){
  char answer[80];
  int status = 1;

  out_i(_("Does the file contain the column names? (%s) "), _("y/n"));
  GETRLINE;
  status = sscanf(line, "%s", answer);
  if (status == 0)
    return;
  if(answer[0] == _("y")[0] || answer[0] == _("Y")[0]){
    has_header = TRUE;
    noheader = FALSE;
  } else
    /* FIXME: can't translate "N" because this string is already used in
     * Frequency table and Compare means */
    if(answer[0] == _("n")[0] || answer[0] == 'N'){
      has_header = FALSE;
      noheader = TRUE;
    } else
      return;

  do{
  out_i(_("Decimal delimiter [%c]: "), dec);
  GETBLINE;
  status = sscanf(line, "%s\n", answer);
  if (status == 0)
    return;
  if(answer[0] == ',' || answer[0] == '.')
    dec = answer[0];
  else
    out_err(WAR, ERR_FILE, ERR_LINE,
	_("Invalid decimal delimiter: '%c'. Please, choose either ',' or '.'"), answer[0]);
  } while(!(answer[0] == ',' || answer[0] == '.'));

  do{
    if(sep){
      if(sep == '\t')
	strcpy(answer, "\\t");
      else
	sprintf(answer, "%c", sep);
    } else
      sprintf(answer, " ,;\\t");
    out_i(_("Field separator ( \\t,;)[%s]: "), answer);
    GETBLINE;
    status = sscanf(line, "%s\n", answer);
    if(line[0] == ' ')
      strcpy(answer, " ");
    if (status == 0)
      return;
    if(answer[0] == ','  || answer[0] == ';' || answer[0] == ' ')
      sep = answer[0];
    else
      if(strcmp(answer, "\\t") == 0)
	sep = '\t';
      else{
	out_err(WAR, ERR_FILE, ERR_LINE, _("Invalid field separator: '%c'"), answer[0]);
	status = 0;
      }
  } while(status == 0);

  out_i(_("What string indicates missing values? [%s]: "), NODATA);
  GETNLINE;
  if(!empty){
    status = sscanf(line, "%s\n", answer);
    if (status == 0)
      return;
    myfree(NODATA);
    NODATA = (char*) mymalloc((strlen(answer) + 1) * sizeof(char));
    strcpy(NODATA, answer);
  }
}

void show_file_head(char *fn){
  FILE *source;
  char *aline;
  char *b;
  int i, rlen;
  size_t blen=64;

  FOPEN(fn, "rt", source);
  aline = (char*)mymalloc(blen);
  rlen = get_line(&aline, &blen, source);
  i = 0;

  out_i(_("First lines of \"%s\":"), fn);
  out_d("\n\n");
  set_winsize();
  b = (char*)m_calloc(SCRCOLS + 4, sizeof(char));
  while (rlen > -1 && i < 10){
    snprintf(b, SCRCOLS, "%s", aline);
    b[SCRCOLS - 3] = '\n';
    b[SCRCOLS - 2] = 0;
    out_d(" %s", b);
    rlen = get_line(&aline, &blen, source);
    i++;
  }
  out_d("\n\n");
  FCLOSE(source);
  myfree(aline);
}

void remove_quotes(char *s){
  int i = 0, j = 0, l;
  l = strlen(s);
  while(i < l){
    if(s[i] != '"'){
      s[j] = s[i];
      j++;
    }
    i++;
  }
  s[j] = 0;
}

int parsecomment(const char *theline, BOOLEAN is_comment) {
  int new = 0, j, n, is_rpt;
  char *s, *token, *comment;
  char ignore[] = " ,;\"\n\t";
  char *var_id;
  
  if(is_comment){
    var_id = (char*)mymalloc(sizeof(char) * 3);
    strcpy(var_id, "#%");/* this char indicates the line contains labels */
  } else{
    var_id = (char*)m_calloc(1, 1);
    var_id[0] = 0;
  }

  comment = (char*)m_calloc(sizeof(char), (strlen(theline) + 1));
  strcpy(comment, theline);

  if (strstr(comment, var_id)!=comment) { /* no valid Var - comment */
    if(strstr(comment, "#!")==comment && strcmp(var_id,"#!")!=0 )
      out_err(WAR, ERR_FILE, ERR_LINE,
	  _("'#!' is an illegal indicator of a column of labels.") );
    return -1;
  }
  s = comment+strlen(var_id);  /* jump over var_id */
  n = ncol;
  while ((token = strtok(s, ignore))!= NULL) {
    s = NULL;
    if (token[0]=='$') {
      if(labelcol == NULL)
	inflate_MCOL();
      labelcol[n_lab] = n;
      n_lab++;
      out_d(_("Label in column %i='%s'\n"), (n+1), token);
    }
    else {
      if(n >= MCOL)
	inflate_MCOL();
      if(alias[n] != NULL)
	myfree(alias[n]);
      is_rpt = FALSE; /* avoiding repeated labels */
      for (j=0; j <n; j++)
	if (strcmp(alias[j], token) == 0)
	  is_rpt = TRUE;
      if (is_rpt){
	alias[n] = (char*) mymalloc((strlen(token)+2) * sizeof(char));
	strcpy(alias[n], token);
	alias[n][strlen(token)] = '_';
	alias[n][strlen(token)+1] = 0;
      }
      else{
	alias[n] = (char*) mymalloc((strlen(token)+1) * sizeof(char));
	strcpy(alias[n], token);
      }
      n++;
      new++;
    }
  }

  if (new == 0 && is_comment) {
    if(silent)
      out_err(FAT, ERR_FILE, ERR_LINE, _("No variables found in comment!"));
    out_err(ERR, ERR_FILE, ERR_LINE, _("No variables found in comment!"));
  }

  return new; /* number of new labels (not counting number of $labels) */
}

/* ==================================================================== */

void delete_last_columns(int i){
  int j;
  for(j = 0; j < i; j++){
    delete_column(ncol - 1);
  }
}

void put_dots(char * s, char c){
  while(*s){
    if(*s == c)
      *s = '.';
    s++;
  }
}

void clean_the_line(char *s){
  char *b = s;
  int i = 0;
  while(*s){
    if(*s == '"')
      i = !i;
    if(i && *s == ',')
      *s = '.';
    s++;
  }

  /* If (dec == ',' && sep == ','), the data is quoted and, thus, there is no
   * need of replacing commas with dots because this was already done. */
  if(dec == ',' && sep != ',')
    put_dots(b, dec);
  /* If sep was defined, and the missing values are quoted, the token will be
   * "\"M\"" */
  if(sep)
    remove_quotes(b);
}

void finish_readingsource(FILE *F, char *s, BOOLEAN try_again){
  FCLOSE(F);
  myfree(s);
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif

  /* If ask_fileformat == TRUE, the user was already asked about the file format */
  char answer[80];
  static int attempts = 0;
  if(try_again && !silent && !ask_fileformat && attempts < 3){
    out_d(_("Statist failed to open the file. Perhaps it\n"
	  "didn't detect the file format correctly.\n"));
    out_i(_("Would you like to set the file format? (%s) "), _("Y/n"));
    GETNLINE;
    out_d("\n");
    if(!empty){
      sscanf(line, "%s", answer);
      if(answer[0] == _("n")[0] || answer[0] == 'N'){
	attempts = 0;
	return;
      }
    }
    attempts++;
    show_file_head(sourcename);
    set_fileformat();
    readsourcefile();
  }
  attempts = 0;
}

void readsourcefile(){
  FILE *source;
  char *aline; /* current line, old line */
  int i, j, newlabs=0, newcol=0, actcol, colread = 0, lread=0, i_lab=0, rlen;
  size_t blen=64;
  REAL  test;
  char ignore[]= " ,;\"\n\t\0", *ptr, *token = NULL;
  BOOLEAN statist_labels = FALSE, header_OK = FALSE;
#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif

  if(ask_fileformat){
    show_file_head(sourcename);
    set_fileformat();
  }

  FOPEN(sourcename, "rt", source);

  aline = (char*)mymalloc(blen);
  rlen = get_line(&aline, &blen, source);
  lread ++;
  if(noheader){ /* Don't scan the first lines looking for a header */
    while (rlen > -1 && (emptyline(aline) || aline[0] == COMMENT)){
      rlen = get_line(&aline, &blen, source);
      lread ++;
    }
  } else{

    /* Seek column labels, skipping true commentaries and empty lines */
    while (rlen > -1 && (emptyline(aline) || aline[0] == COMMENT)){
      if(aline[0] == COMMENT && (aline[1] == '%' || aline[1] == '!')){
	newlabs = parsecomment(aline, TRUE);
	if(newlabs)
	  statist_labels = TRUE;
      }
      rlen = get_line(&aline, &blen, source);
      lread ++;
    }
    if(rlen == -1){
      if(silent){
	out_err(FAT, ERR_FILE, ERR_LINE,
	    _("Couldn't find data in file \"%s\"!"), sourcename);
      } else{
	out_err(ERR, ERR_FILE, ERR_LINE,
	    _("Couldn't find data in file \"%s\"!"), sourcename);
	finish_readingsource(source, aline, FALSE);
	return;
      }
    }

    /* Read the first line as if it contains the column names because the user
     * passed the command line option --header */
    if(!statist_labels && has_header){
      newlabs = (parsecomment(aline, FALSE));
      if(newlabs > 0)
	header_OK = TRUE;
      for(i = 0; i < newlabs; i++){
	j = i + ncol;
	if(!((alias[j][0] >= 'A' && alias[j][0] <= 'Z') 
	      || (alias[j][0] >= 'a' && alias[j][0] <= 'z'))){
	  out_err(WAR, ERR_FILE, ERR_LINE,
	      _("Name of column %d doesn't begin with an ascii letter: \"%s\"."),
	      j+1, alias[j]);
	  break;
	}
      }
    } else{
      /* Read the first line as if it contains the column names, but drop the
       * names if they don't appear to be valid ones. */
      if(!statist_labels && detect_header){
	newlabs = (parsecomment(aline, FALSE));
	if(newlabs > 0){
	  header_OK = TRUE;
	  for(i = 0; i < newlabs; i++){
	    j = i + ncol;
	    if(!((alias[j][0] >= 'A' && alias[j][0] <= 'Z') 
		  || (alias[j][0] >= 'a' && alias[j][0] <= 'z'))){
	      for(j = ncol; j < (ncol+newlabs); j++){
		myfree(alias[j]);
		alias[j] = NULL;
	      }
	      header_OK = FALSE;
	      newlabs = 0;
	      break;
	    }
	  }
	}
      }
    }

    if(newlabs > 0){
      for(i = ncol; i < (ncol+newlabs); i++)
	out_d(_("Column %i = %s\n"), i+1, alias[i]);
    }

    /* Read another line if we successfully read the column names, although
     * the data file doesn't have the "#%" string. */
    if(header_OK){
      rlen = get_line(&aline, &blen, source);
      lread ++;
      if(verbose && detect_header && !has_header){
	out_d("\n");
	if(newlabs == 1)
	  out_d(_("One valid column name found!"));
	else
	  out_d(_("Valid column names found!"));
	out_d("\n\n");
      }
    }
  } /* End of "if(noheader)" */

  /* Parse the first line of data, but don't put the data in the temp files
     yet. Just count the number of columns. */
  char *linecopy = (char*)m_calloc(2, blen);
  strcpy(linecopy, aline);
  clean_the_line(linecopy);
  
  ptr = linecopy;
  if(sep){  /* A field separator was defined. */
    newcol = 1;
    j = strlen(ptr);
    for(i = 0; i < j; i++)
      if(ptr[i] == sep)
	newcol++;
  } else{
    while ((token = strtok(ptr, ignore))!= NULL) {
      ptr = NULL;
      while ((i_lab<n_lab) && (newcol==labelcol[i_lab])) {
	i_lab++;
	token = strtok(ptr, ignore);
      }
      if (token==NULL) {
	break;
      }
      if ((strcmp(token, NODATA) == 0) || (sscanf(token, "%lf", &test) == 1)) {
	newcol++;
      }
      else {
	if(silent){
	  out_err(FAT, ERR_FILE, ERR_LINE,
	      _("Illegal format of value '%s' in line %i!\n"
		"Couldn't read file %s!"), token, lread, sourcename);
	} else{
	  out_err(ERR, ERR_FILE, ERR_LINE,
	      _("Illegal format of value '%s' in line %i!\n"
		"Couldn't read file %s!"), token, lread, sourcename);
	  finish_readingsource(source, aline, TRUE);
	  return;
	}
      }
    }
  }

  /* Check whether n_columns == n_labels     */
  if ((newlabs != 0) && (newlabs != newcol)) {
    if(silent){
      out_err(FAT, ERR_FILE, ERR_LINE,
	  _("Number of columns (%d) does not equal number of labels (%d)!"),
	  newcol, newlabs);
    } else{
      out_err(ERR, ERR_FILE, ERR_LINE,
	  _("Number of columns (%d) does not equal number of labels (%d)!"),
	  newcol, newlabs);
      finish_readingsource(source, aline, TRUE);
      return;
    }
  }

  create_columns(newcol);

  /* Finally, read data */
  out_d(_("Reading %i columns ...\n"), newcol);
  BOOLEAN endofline;
  if(sep)
    token = (char*)mymalloc(256);
  do {
    if ((!emptyline(aline)) && (aline[0]!=COMMENT)) {
      colread = 0;
      i_lab = 0;
      clean_the_line(aline);
      ptr = aline;
      j = 0;
      endofline = FALSE;

      while (!endofline) {
	if(sep){
	  i = 0;
	  while(!(*ptr == sep || *ptr == '\n')){
	    token[i] = *ptr;
	    i++;
	    ptr++;
	  }
	  token[i] = 0;
	  if(i == 0)
	    strcpy(token, NODATA);
	  if(*ptr == '\n'){
	    endofline = TRUE;
	  }
	  ptr++;
	  actcol = j + (ncol-newcol);
	} else{
	  if((token = strtok(ptr, ignore)) == NULL)
	    break;
	  ptr = NULL;
	  actcol = j + (ncol-newcol);
	  while ((i_lab<n_lab) && (j==labelcol[i_lab])) {
	    i_lab++;
	    token = strtok(ptr, ignore);
	  }
	  if (token==NULL) {
	    break;
	  }
	}
	if (j>=newcol) {
	  if(silent){
	    out_err(FAT, ERR_FILE, ERR_LINE,
		_("Too many columns in row %i. (%d columns)"), lread, j+1);
	  } else{
	    out_err(ERR, ERR_FILE, ERR_LINE,
		_("Too many columns in row %i. (%d columns)"), lread, j+1);
	    delete_last_columns(newcol);
	    finish_readingsource(source, aline, TRUE);
	    return;
	  }
	}

	if (strcmp(token, NODATA) == 0 || strcmp(token, "nan") == 0){
	  FWRITE(&SYSMIS, sizeof(REAL), 1, tmpptr[actcol]);
	  nn[actcol] ++;
	  colread ++;
	}
	else if (sscanf(token, "%lf", &test)==1) {
	  FWRITE(&test, sizeof(REAL), 1, tmpptr[actcol]);
	  nn[actcol] ++;
	  colread ++;
	}
	else {
	  if(silent){
	    out_err(FAT, ERR_FILE, ERR_LINE,
		_("Illegal format of value '%s' in line %i!"), token, lread);
	  } else{
	    out_err(ERR, ERR_FILE, ERR_LINE,
		_("Illegal format of value '%s' in line %i!"), token, lread);
	    delete_last_columns(newcol);
	    finish_readingsource(source, aline, TRUE);
	    return;
	  }
	}
	j++;
      }
    }

    if (colread != newcol) {
      if(silent){
	out_err(FAT, ERR_FILE, ERR_LINE,
	    _("Row %i contains just %i instead of %i columns!"),
	    (lread), colread, newcol);
      } else{
	out_err(ERR, ERR_FILE, ERR_LINE,
	    _("Row %i contains just %i instead of %i columns!"),
	    (lread), colread, newcol);
	delete_last_columns(newcol);
	finish_readingsource(source, aline, TRUE);
	return;
      }
    }
    rlen = get_line(&aline, &blen, source);
    lread ++;
  } while (rlen != -1);

  out_d(_("\nRead data sets: \n") );
  for (j=0; j<newcol; j++) {
    actcol = j + (ncol-newcol);
    out_d(_("Column %s: %i\n"), alias[actcol], nn[actcol]);
  }

  finish_readingsource(source, aline, FALSE);
}

/* ==================================================================== */


void newsourcefile() {
   char answer[3], newsourcename[80];
   FILE *source;

   ls();
   out_i(_("Name of data file: ") );
   GETRLINE;
   sscanf(line, "%s", newsourcename);
   out_d("\n\n");

   while ((source = fopen(newsourcename,"rt")) == NULL)
   {
      out_i(_("File \"%s\" not found!\n"), newsourcename);
      out_i(_("Please enter new file name: ") );
      GETRLINE;
      sscanf(line, "%s", newsourcename);
      out_d("\n");
   }
   FCLOSE(source);
   show_file_head(newsourcename);
   if(ncol > 0){
     out_i(_("Shall the old data be removed? (%s) "), _("y/N") );
     GETNLINE;
     if(!(empty)){
       sscanf(line, "%s", answer);
       if (answer[0] == _("y")[0]  || answer[0] == _("Y")[0]) {
	 erasetempfiles();
       }
     }
   }
   if(sourcename)
     myfree(sourcename);
   sourcename = (char*) mycalloc(strlen(newsourcename) + 1, sizeof(char));
   strcpy(sourcename, newsourcename);
   readsourcefile();
   if (log_set) {
     fprintf(logfile, "-----------------------------------------------------\n");
     fprintf(logfile,
     	_("\nNew source file: %s\n\n") , sourcename);
   }
   attach_labels_to_columns();
 }


/* =================================================================== */

int getcols(int min, int max, BOOLEAN eraserow){
  char salias[80];
  int i, j, w, nc, nr;
  BOOLEAN inputok, found;
  salias[79] = 0;
  if(ncol == 0)
    return 0;
  if(ncol < min){
    if(ncol == 1)
      strncpy(salias, _("but this data file has just 1 column!"), 79);
    else
      snprintf(salias, 80, _("but this data file has only %i columns!"), ncol);
    out_err(ERR, ERR_FILE, ERR_LINE,
	_("This analysis requires at least %i columns,\n  %s"), min, salias);
    return 0;
  }

  out_d("\n");
  out_d(_("Columns: ") );
  if(format_columns_out){
    set_winsize();
    out_d("\n");
    w = 0;
    for(j = 0; j < ncol; j++)
      if(strlen(alias[j]) > w)
	w = strlen(alias[j]);
    w += 2;
    nc = SCRCOLS / w;
    nr = 1 + (ncol / nc);
    if(ncol > nc && (ncol % nc) != 0)
      nr++;
    snprintf(salias, 80, "%%-%is", w);
    for(i = 0; i < nr; i++){
      for(j = 0; j < nc; j++){
	w = nc * i + j;
	if(w < ncol)
	  out_d(salias, alias[w]);
      }
      out_d("\n");
    }
  } else{
    for (j = 0; j < ncol; j++) {
      out_d("%s ", alias[j]);
    }
  }
  out_d("\n");

  i = 0;
  do{
    inputok = FALSE;
    found  = FALSE;
    while (!inputok) {
      if(max > 1)
	out_i(_("Column for variable %i: "), (i+1));
      else
	out_i(_("Column name: "));
      GETBLINE;
      sscanf(line, "%s", salias);

      if (strcmp(line, _(_ALL_)) == 0) {
	if(max < ncol){
	  out_err(ERR, ERR_FILE, ERR_LINE,
	      _("Please, choose at most %i columns!"), max);
	  break;
	}
	for (j = 0; j < ncol; j++) {
	  acol[j] = j;
	}
	alloc_cols(ncol, eraserow);
	return ncol;
      }

      /* check if column name is matched exactly */
      for (j=0; j<ncol; j++) {
	if (strcmp(alias[j], salias)==0) {
	  acol[i] = j;
	  inputok = TRUE;
	  i++;
	  break;
	}
      }
      if (inputok) {		
	break;		/* exact column alias entered -> go on 	*/
      }

      /* try to complete entered column alias	*/
      for (j=0; j<ncol; j++) {
	if (str_in_str(alias[j], salias)) {
	  if (found) {
	    out_err(ERR, ERR_FILE, ERR_LINE,
		_("Column name '%s' is not unique!"), salias);
	    inputok = FALSE;
	    i--;
	    break;
	  }
	  else {
	    found = TRUE;
	    inputok = TRUE;
	    acol[i] = j;
	    i++;
	  }
	}
      }
      if ((!inputok) && (!found)) {
	out_err(ERR, ERR_FILE, ERR_LINE,
	    _("Column %s does not exist!"), salias);
      }
      else if ((!inputok) && (found)) {
	found = FALSE;
      }
    }
  } while(!empty && i < max);

  if(i < min){
    if(i > 0)
      out_err(ERR, ERR_FILE, ERR_LINE,
	  _("At least %i columns have to be selected!"), min);
    return 0;
  }
  if(eraserow == 1 && !(equal_rows(i))){
    out_err(ERR, ERR_FILE, ERR_LINE, _("The columns must have "
	  "the same number of data points for this analysis!"));
    return 0;
  }

  out_d("\n");
  if (log_set) {
    fprintf(logfile, "-----------------------------------------------------------\n\n");
  }

  alloc_cols(i, eraserow);
  return(i);
}

void printcols() {
  int i, j, k, n, r, p, q, w;
  char b[50], b2[50], *header;
  BOOLEAN labelfound = FALSE;
  b[49] = 0;

  /* Choosing columns */
  set_winsize();
  k = (SCRCOLS / 11) - 1;
  if(ncol < k){ /* show all columns if the screen is width enough */
    for (j = 0; j < ncol; j++) {
      acol[j] = j;
    }
    alloc_cols(ncol, 3);
    n = ncol;
  } else{
    n = getcols(1, ncol, 3);
  }
  if(n == 0)
    return;

  /* Determining number of rows */
  k = nn[acol[0]];
  j = strlen(alias[acol[0]]) + 10;
  for(i = 1; i < n; i++){
    j += strlen(alias[acol[i]]);
    if(nn[acol[i]] > k)
      k = nn[acol[i]];
  }

  /* Creating header */
  if(j < ((n + 1) * 16))
    j = (n + 1) * 16;
  header = (char*)m_calloc(j, sizeof(char));
  strcpy(header, "       ");
  for(i = 0; i < n; i++){
    if(names[acol[i]] && names[acol[i]]->n > 0)
      snprintf(b, 50, "%-10s ", alias[acol[i]]);
    else
      snprintf(b, 50, "%10s ", alias[acol[i]]);
    strcat(header, b);
  }
  strcat(header, "\n");
 
  /* Printing data */
  out_r(_("Data from columns:\n"));
  colorize(ClHeader);
  out_r(header);
  colorize(ClDefault);
  j = 0;
  p = 2; /* already printed lines = column_names + wait_message */
  int sz;
  while(j < k){
    colorize(ClLineNum);
    out_r("%5i: ", (j + 1));
    colorize(ClDefault);
    for(i = 0; i < n; i++){
      if(j < nn[acol[i]]){
	if(xx[acol[i]][j] == SYSMIS)
	  out_r("%10c ", '.');
	else{
	  if(names[acol[i]])
	    for(q = 0; q < names[acol[i]]->n; q++)
	      if(names[acol[i]]->v[q] == xx[acol[i]][j]){
		labelfound = TRUE;
		strncpy(b, names[acol[i]]->l[q], 10);
		b[10] = 0;
		sz = 10;
		if(is_utf8){
		  /* If there are non-ascii chars, we truncated the label prematurely*/
		  while(stringLen(b) < 10 && sz < 40){
		    sz++;
		    strncpy(b, names[acol[i]]->l[q], sz);
		  }

		  /* Avoiding truncating multibyte char, at least in Latin 1. */
		  w = strlen(b);
		  if(b[w-1] == (char)0xC3){
		    sz++;
		    strncpy(b, names[acol[i]]->l[q], sz);
		  }
		  sz = 10 + strlen(b) - stringLen(b);
		}

		sprintf(b2, "%%-%ds ", sz);
		out_r(b2, b);
	      }
	  if(labelfound)
	    labelfound = FALSE;
	  else
	    out_r("%10g ", xx[acol[i]][j]);
	}
      } else{
	out_r("%10s ", " ");
      }
    }
    out_r("\n");
    p++;
    if (p == (SCRLINES - 1) && !(silent)){
      p = 2;
      out_i(_("---> Please, choose: <RETURN> to continue,\n"
	    "     <Any letter> to stop, or a row number: ") );
      GETNLINE
	if(!empty){
	  if((line[0] >= 'a' && line[0] <= 'z') 
	      || (line[0] >= 'A' && line[0] <= 'z'))
	    return;
	  r = getint();
	  if(r > 0)
	    j = r - 2;
	}
      colorize(ClHeader);
      out_r(header);
      colorize(ClDefault);
    }
    j++;
  }
}

void printcol(REAL x[], int n) {
   int i, k;
   out_r(_("Data from column \"%s\":\n"), get_label(x));
   for (i=0; i<n; i++) {
     k=i+1;
     if (x[i] == SYSMIS)
       out_r("%5i.)  %s\n", k, NODATA);
     else
       out_r("%5i.)  %g\n", k, x[i]);
     if ((i+1) % (SCRLINES - 1) == 0) {
	mywait();
	if (!empty) {
	  return;
	}
     }
   }
   out_r("-------------------------------------------\n\n");
}

/* =================================================================== */


PREAL readcol(int i) {
  PREAL px;

  if (nn[i] == 0) {
    out_err(FAT, ERR_FILE, ERR_LINE,
    	_("Column %i does not exist!"), i+1);
  }
  px = (REAL*)mycalloc(nn[i], sizeof(REAL));
  rewind(tmpptr[i]);
  FREAD(px, sizeof(REAL), nn[i], tmpptr[i]);
  x_read[i] = TRUE;
  return px;
}

/* ==================================================================== */


void alloc_cols(int n_alloc, BOOLEAN eraserow) {
  int k;
  int cr = 0; /* current row */
  int tr = 0; /* total number of rows already checked */
  BOOLEAN RowHasMis = FALSE;

  /* delete all columns from memory */
  for (k=0; k<MCOL; k++){
    if((x_read[k])){
      free_column(k);
    }
  }

  /* put selected columns in memory */
  for (k=0; k<n_alloc; k++)
    if (!x_read[acol[k]]){
      xx[acol[k]] = readcol(acol[k]);
    }

  /* Delete rows with missing values or simply delete missing values */
  /* if eraserow == 3, do nothing                                    */
  if (eraserow == TRUE){
    while(tr < nn[acol[0]]){
      for (k=0; k<n_alloc; k++)
        if(xx[acol[k]][tr] == SYSMIS) RowHasMis = TRUE;
      if (RowHasMis){
        tr++;
        RowHasMis = FALSE;
      }
      else{
        for (k=0; k<n_alloc; k++)
          xx[acol[k]][cr] = xx[acol[k]][tr];
        cr++;
        tr++;
      }
    }
    for (k=0; k<n_alloc; k++)
      vn[acol[k]] = cr;
    out_r( _("%d rows with missing values were deleted for this analysis\n\n"),
          (nn[acol[0]] - cr));
  }
  else if (eraserow == FALSE) {
    for (k=0; k<n_alloc; k++){
      tr = 0;
      cr = 0;
      while (tr < nn[acol[k]]){
        if(xx[acol[k]][tr] == SYSMIS)
          tr++;
        else {
          xx[acol[k]][cr] = xx[acol[k]][tr];
          cr++;
          tr++;
        }
      }
      vn[acol[k]] = cr;
      out_r( _("Column %s: %d data points\n"),
          alias[acol[k]], cr);
    }
  }

  if (log_set)
    for (k=0; k<n_alloc; k++)
      fprintf(logfile, _("Variable %i = Column %s\n"), (k+1), alias[acol[k]] );

  /* rewinding of pointers to tmpfiles */
  for (k=0; k<n_alloc; k++)
    rewind(tmpptr[acol[k]]);
}


BOOLEAN make_new_col(char *analias, int n) {
  int i;

  for (i=0; i<ncol; i++) {
    if (strcmp(analias, alias[i])==0) {
      out_err(ERR, ERR_FILE, ERR_LINE,
		_("Column %s exists already!"), analias);
      return FALSE;
    }
  }
  create_columns(1);
  if(alias[ncol - 1])
    myfree(alias[ncol - 1]);
  alias[ncol - 1] = (char*)mymalloc((strlen(analias)+1));
  strcpy(alias[ncol - 1], analias);
  out_r(_("New column %s created!\n"), alias[ncol - 1]);
  nn[ncol - 1] = n;
  return TRUE;
}

 /* =================================================================== */

int col_exist(char *analias, BOOLEAN is_error) {
  int i;

  for ( i=0; i<ncol; i++ ) {
    if ( alias[i] && strcmp(analias, alias[i])==0)  {
      if(is_error)
	out_err(ERR, ERR_FILE, ERR_LINE,
	    _("Column %s exists already!"), analias);
      return i;
    }
  }
  return -1;
}

/* =================================================================== */

/* =================================================================== */

#ifndef STATIST_X
char *get_label(PREAL x) {
  int i;

  for (i=0; i<ncol; i++) {
    if (x == xx[i]) {
      if(names[i] && names[i]->ctitle)
	return names[i]->ctitle;
      else
	return alias[i];
    }
  }
  out_err(ERR, ERR_FILE, ERR_LINE,
  	_("No label found for column!") );
  return NULL;
}
#endif

char *get_name(PREAL x) {
  int i;

  for (i=0; i<ncol; i++) {
    if (x == xx[i] && alias[i]){
	return alias[i];
    }
  }
  out_err(ERR, ERR_FILE, ERR_LINE,
  	_("No name found for column!") );
  return NULL;
}

void log_transform() {
  char analias[80];
  PREAL y;
  int i, n = 0;

  out_i(_("Please select column for log-transformation\n") );
  i = getcols(1, 1, 3);
  if(i == 0)
    return;
  strncpy(analias, "log_", 79);
  strncat(analias, alias[acol[0]], 79-strlen(analias));
  if(col_exist(analias, TRUE) != -1)
    return;
  y = (REAL*)m_calloc(nn[acol[0]], sizeof(REAL));
  for (i=0; i<nn[acol[0]]; i++) {
    if (xx[acol[0]][i] > 0.0)
      y[i] = log10(xx[acol[0]][i]);
    else{
      y[i] = SYSMIS;
      if(xx[acol[0]][i] != SYSMIS)
	n++;
    }
  }

  if(n == 1)
    out_err(MWA, ERR_FILE, ERR_LINE, _("One value was less or equal to zero"
	  " and was transformed into missing value!"));
  if(n > 1)
    out_err(MWA, ERR_FILE, ERR_LINE, _("%i values were less or equal to zero"
	  " and were transformed into missing values!"), n);

  if (!(make_new_col(analias, nn[acol[0]]))) {
    return;
  }
  FWRITE(y, sizeof(REAL), nn[acol[0]], tmpptr[ncol - 1]);
}


void ln_transform() {
  char analias[80];
  PREAL y;
  int i, n = 0;

  out_i(_("Please select column for log-transformation\n") );
  i = getcols(1, 1, 3);
  if(i == 0)
    return;
  strncpy(analias, "ln_", 79);
  strncat(analias, alias[acol[0]], 79-strlen(analias));
  if(col_exist(analias, TRUE) != -1)
    return;
  y = (REAL*)m_calloc(nn[acol[0]], sizeof(REAL));
  for (i=0; i<nn[acol[0]]; i++) {
    if (xx[acol[0]][i] > 0.0)
      y[i] = log(xx[acol[0]][i]);
    else{
      y[i] = SYSMIS;
      if(xx[acol[0]][i] != SYSMIS)
	n++;
    }
  }

  if(n == 1)
    out_err(MWA, ERR_FILE, ERR_LINE, _("One value was less or equal to zero"
	  " and was transformed into missing value!"));
  if(n > 1)
    out_err(MWA, ERR_FILE, ERR_LINE, _("%i values were less or equal to zero"
	  " and were transformed into missing values!"), n);

  if (!(make_new_col(analias, nn[acol[0]]))) {
    return;
  }
  FWRITE(y, sizeof(REAL), nn[acol[0]], tmpptr[ncol - 1]);
}
/* =================================================================== */


void power_10_transform() {
  char analias[80];
  PREAL y;
  int i;

  out_i(_("Please select column for exponentiation\n") );
  i = getcols(1, 1, 3);
  if(i == 0)
    return;
  strncpy(analias, "10^_", 79);
  strncat(analias, alias[acol[0]], 79-strlen(analias));
  if(col_exist(analias, TRUE) != -1)
    return;
  y = (REAL*)m_calloc(nn[acol[0]], sizeof(REAL));
  for (i=0; i<nn[acol[0]]; i++) {
    if(xx[acol[0]][i] == SYSMIS)
      y[i] = SYSMIS;
    else
      y[i] = pow(10.0, xx[acol[0]][i]);
  }

  if (!(make_new_col(analias, nn[acol[0]]))){
    return;
  }
  FWRITE(y, sizeof(REAL), nn[acol[0]], tmpptr[ncol - 1]);
}

void power_e_transform() {
  char analias[80];
  PREAL y;
  int i;

  out_i(_("Please select column for exponentiation\n") );
  i = getcols(1, 1, 3);
  if(i == 0)
    return;
  strncpy(analias, "e^_", 79);
  strncat(analias, alias[acol[0]], 79-strlen(analias));
  if(col_exist(analias, TRUE) != -1)
    return;
  y = (REAL*)m_calloc(nn[acol[0]], sizeof(REAL));
  for (i=0; i<nn[acol[0]]; i++) {
    if(xx[acol[0]][i] == SYSMIS)
      y[i] = SYSMIS;
    else
      y[i] = exp(xx[acol[0]][i]);
  }

  if (!(make_new_col(analias, nn[acol[0]]))){
    return;
  }
  FWRITE(y, sizeof(REAL), nn[acol[0]], tmpptr[ncol - 1]);
}


/* =================================================================== */


void inv_transform() {
  char analias[80];
  PREAL y;
  int i;

  out_i(_("Please select column for inversion\n") );
  i = getcols(1, 1, 3);
  if(i == 0)
    return;
  strncpy(analias, "inv_", 79);
  strncat(analias, alias[acol[0]], 79-strlen(analias));
  if(col_exist(analias, TRUE) != -1)
    return;
  y = (REAL*)m_calloc(nn[acol[0]], sizeof(REAL));
  for (i=0; i<nn[acol[0]]; i++) {
    if(xx[acol[0]][i] == SYSMIS)
      y[i] = SYSMIS;
    else
      y[i] = 1./xx[acol[0]][i];
  }

  if (!(make_new_col(analias, nn[acol[0]]))){
    return;
  }
  FWRITE(y, sizeof(REAL), nn[acol[0]], tmpptr[ncol - 1]);
}


/* =================================================================== */


void z_transform() {
  char analias[80];
  PREAL y;
  REAL mean, sdv;
  int i;

  out_i(_("Please select column for z-transformation\n") );
  i = getcols(1, 1, TRUE);
  if(i == 0)
    return;
  strncpy(analias, "z_", 79);
  strncat(analias, alias[acol[0]], 79-strlen(analias));
  if(col_exist(analias, TRUE) != -1)
    return;
  y = (REAL*)m_calloc(nn[acol[0]], sizeof(REAL));
  sdv = get_sdv(xx[acol[0]], nn[acol[0]]);
  mean = get_mean(xx[acol[0]], nn[acol[0]]);
  if(nn[acol[0]] != vn[acol[0]])
    alloc_cols(1, 3);
  for (i=0; i<nn[acol[0]]; i++) {
    if(xx[acol[0]][i] == SYSMIS)
      y[i] = SYSMIS;
    else
      y[i] = (xx[acol[0]][i]-mean)/sdv;
  }

  if (!(make_new_col(analias, nn[acol[0]]))){
    return;
  }
  FWRITE(y, sizeof(REAL), nn[acol[0]], tmpptr[ncol - 1]);
}


/* =================================================================== */

void sort_col() {
  char analias[80];
  PREAL y;
  int i;

  out_i(_("Please select column to be sorted\n") );
  i = getcols(1, 1, 3);
  if(i == 0)
    return;
  strncpy(analias, "sort_", 79);
  strncat(analias, alias[acol[0]], 79-strlen(analias));
  if(col_exist(analias, TRUE) != -1)
    return;
  y = (REAL*)m_calloc(nn[acol[0]], sizeof(REAL));
  for (i=0; i<nn[acol[0]]; i++) {
    y[i] = xx[acol[0]][i];
  }
  qsort(y, nn[acol[0]], sizeof(REAL), real_compar_up);
  if (!(make_new_col(analias, nn[acol[0]]))){
    return;
  }
  FWRITE(y, sizeof(REAL), nn[acol[0]], tmpptr[ncol - 1]);
}


/* =================================================================== */


void readcol_from_term() {
  char aline[80], answer[10];
  int n=0;
  REAL temp;
  BOOLEAN ok, stop=FALSE;

  if (ncol > 0) {
    out_i(_("Shall all data be deleted? (%s) "), _("y/N") );
    GETNLINE;
    if(!(empty)){
      sscanf(line, "%s", answer);
      if (answer[0] == _("y")[0] || answer[0] == _("Y")[0]) {
	erasetempfiles();
      }
    }
  }

  out_i(_("Column %i is being read, stop input with '.'\n"), (ncol+1));
  aline[0] = '1';
  create_columns(1);

  while (!stop) {
    ok = FALSE;
    while (!ok) {
      out_d(_("Value %i: "), (n+1));
      fgets(aline, 79, stdin);
      if ( (aline[0]=='.') && (strlen(aline)==2) ) {
	stop = TRUE;
      }
      if ( (sscanf(aline, "%lf", &temp)==1) || (stop) ) {
	ok = TRUE;
      }
      else {
	out_err(ERR, ERR_FILE, ERR_LINE,
	  _("Illegal input, please repeat: ") );
      }

      if ( (ok) && (!stop) ) {
	n++;
	FWRITE(&temp, sizeof(REAL), 1, tmpptr[ncol - 1]);
      }
    }
  }
  if (n>0) {
    nn[ncol - 1] = n;
  } else{
    delete_column(ncol - 1);
  }
}



BOOLEAN str_in_str(const char *s1, const char *s2) {
  int i, n = strlen(s2);

  for (i=0; i<n; i++) {
    if (s1[i] != s2[i]) {
      return FALSE;
    }
  }
  return TRUE;
}


BOOLEAN emptyline(const char *s) {
  int i, n = strlen(s);

  for (i=0; i<n; i++) {
    if (!isspace((int)s[i])) {
      return FALSE;
    }
  }
  return TRUE;
}


BOOLEAN formatToken(char *token, char *result){
  int i = 0;
  REAL test;

  /* Remove leading blanks */
  while(*token == ' ' && *token != 0)
    token++;
  if(*token == 0 || *token == '\r' || *token == '\n'){
    sprintf(result, "%s", NODATA);
    return FALSE;
  }

  /* Get the token, including blanks */
  while(*token != 0){
    result[i] = *token;
    token++;
    i++;
  }

  /* Remove trailing blanks */
  i--;
  while(result[i] == ' '){
    result[i] = 0;
    i--;
  }
  i++;

  result[i] = 0;
  if(sscanf(result, "%lf", &test) == 1){
    if(test == floor(test) && test >= -9999999999999999.0 && test <= 9999999999999999.0)
      snprintf(result, 32, "%-16.16g", test);
    else
      snprintf(result, 32, "%-16.10g", test);
    token = result;
    while(*token != 0){
      token++;
      if(*token == ' '){
	*token = 0;
	break;
      }
    }
    return FALSE;
  } else{
    strcpy(token, result);
    sprintf(result, "\"%s\"", token);
    return TRUE;
  }
}

void xcols_usage(char * name){
  out_d(_("\nThe option --xcols tells Statist to extract columns from a fixed "
      "width data file.\n\n"));
  out_d(_("\nUsage:\n"
      "%s --xcols config_file data_base dest_file\n\n"), name);
  exit(0);
}

void extract_cols(int argc, char *argv[]){
  FILE *f1, *f2;
  char b[1000];
  int i, j, k, n = 0, nrows = 0, rlen, max = 100, pos = 0, *begin, *end;
  size_t blen = 64;
  char *b1, **lbel, *token, *ftokn;
  unsigned int l;
  BOOLEAN *alpha;

#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif

  if(sep == 0)
    sep = ' ';

  for(i = 1; i < argc; i++)
    if(strcmp(argv[i], "--xcols") == 0){
      pos = i;
      break;
    }

  if((argc - pos) < 4)
    xcols_usage(argv[0]);

  lbel = (char**) mymalloc(max * sizeof(char*));
  begin = (int*) mymalloc(max * sizeof(int));
  end = (int*) mymalloc(max * sizeof(int));
  alpha = (BOOLEAN*)mymalloc(max * sizeof(BOOLEAN));
  b1 = (char*) mymalloc(blen * sizeof(char));

  /* read config_file */
  FOPEN(argv[pos + 1], "r", f1);
  rlen = get_line(&b1, &blen, f1);
  while (rlen != -1){
    if(b1[0] == '#' || strlen(b1) < 3){
      rlen = get_line(&b1, &blen, f1);
      continue;
    }
    i = 0;
    while(i < 997 && b1[i] != ' ' && b1[i] != '\t'){
      b[i] = b1[i];
      i++;
    }
    b[i] = 0;
    lbel[n] = (char*) mymalloc(i + 1);
    strcpy(lbel[n], b);
    while(!(b1[i] >= '0' && b1[i] <= '9'))
      i++;
    j = 0;
    while(b1[i] >= '0' && b1[i] <= '9'){
      b[j] = b1[i];
      i++; j++;
    }
    b[j] = 0;
    begin[n] = atoi(b) - 1;
    while(!((b1[i] >= '0' && b1[i] <= '9') || b1[i] == '\n'))
      i++;
    j = 0;
    while(b1[i] >= '0' && b1[i] <= '9'){
      b[j] = b1[i];
      i++; j++;
    }
    b[j] = 0;
    if(b[0])
      end[n] = atoi(b) - 1;
    else
      end[n] = begin[n];
    n++;
    if(n == max){
      max += 100;
      lbel = (char**) myrealloc(lbel, (max * sizeof(char*)));
      begin = (int*) myrealloc(begin, (max * sizeof(int)));
      end = (int*) myrealloc(end, (max * sizeof(int)));
      alpha = (BOOLEAN*) myrealloc(alpha, (max * sizeof(BOOLEAN)));
    }
    rlen = get_line(&b1, &blen, f1);
  }
  FCLOSE(f1);

  j= 0;
  l = 3;
  for(i = 0; i < n; i++){
    alpha[i] = FALSE;
    j += strlen(lbel[i]) + 30;
    if(j > l)
      l = j;
    if((end[i] - begin[i] + 1) > l)
      l = end[i] - begin[i] + 1;
  }
  if(l < 128)
    l = 128;
  else
    l *= 3;
  token = (char*)mymalloc(l * sizeof(char));
  ftokn = (char*)mymalloc(l * sizeof(char));

  /* read from origin, and write to destination */
  FOPEN(argv[pos + 2], "r", f1);
  FOPEN(argv[pos + 3], "w", f2);
  out_d(_("Extracting columns from \"%s\" to \"%s\"...\n"),
      argv[pos + 2], argv[pos + 3]);
  /* Don't put the "#%" string in the first line if the user doesn't seem
   * to use it. */
  if(!(has_header || detect_header))
    fprintf(f2, "#%%");
  for(i = 0; i < (n - 1); i++)
    fprintf(f2, "%s%c", lbel[i], sep);
  fprintf(f2, "%s\n", lbel[n-1]);
  rlen = get_line(&b1, &blen, f1);
  while (rlen != -1){
    for(i = 0; i < n; i++){
      k = 0;
      for(j = begin[i]; j <= end[i]; j++){
	token[k] = b1[j];
	k++;
      }
      token[k] = 0;
      if(formatToken(token, ftokn))
	alpha[i] = TRUE;
      if(i < (n - 1))
	fprintf(f2, "%s%c", ftokn, sep);
      else
	fprintf(f2, "%s\n", ftokn);
    }
    rlen = get_line(&b1, &blen, f1);
    nrows++;
  }
  FCLOSE(f1);
  FCLOSE(f2);
  myfree(b1);
  myfree(begin);
  myfree(end);
  out_d(_("Done: %d columns, %d rows.\n"), n, nrows);
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  j = 0;
  for(i = 0; i < n; i++)
    if(alpha[i])
      j = 1;
  if(j){
    out_err(WAR, ERR_FILE, ERR_LINE,
	_("Non-numeric values were found."));
    out_r(_("List of columns with non-numeric values:\n"));
    for(i = 0; i < n; i++)
      if(alpha[i])
	out_r(" %s", lbel[i]);
    out_r("\n");
  }
  for(i = 0; i < n; i++)
    myfree(lbel[i]);
  myfree(lbel);
  myfree(alpha);
}

void xsample_usage(char * name){
  out_d(_("\nThe option --xsample tells Statist to extract a random sample of\n"
	"rows from a given data file.\n\n"));
  out_d(_("Usage:\n\n"
      "  %s --xsample percentage data_base dest_file\n\n"
      "where \"percentage\" is an integer between 1 and 99.\n\n"), name);
  exit(1);
}

void extract_sample(int argc, char * argv[]){
  int percent = -1;
  char *s;
  int i, k, n = 0, N = 0, rlen, pos = 0;
  FILE * f1;
  FILE * f2;
  size_t blen = 64;
#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif

  for(i = 1; i < argc; i++)
    if(strcmp(argv[i], "--xsample") == 0){
      pos = i;
      break;
    }

  if((argc - pos) < 4)
    xsample_usage(argv[0]);
  percent = atoi(argv[pos  + 1]);
  if(percent > 99 || percent < 1){
    out_err(ERR, ERR_FILE, ERR_LINE,
	_("\"%s\" is not a valid value for percentage."), argv[pos + 1]);
    xsample_usage(argv[0]);
  }

  s = (char*)mymalloc(blen);

  /* read from source, and write to destine */
  srand(time(NULL));
  k = percent * 10;
  FOPEN(argv[pos + 2], "r", f1);
  FOPEN(argv[pos + 3], "w", f2);

  out_r(_("Creating a new database with a random sample of approximately\n"
      "%i%% of \"%s\" rows...\n"), percent, argv[3]);

  rlen = get_line(&s, &blen, f1);
  while(rlen != -1 && (s[0] == '#' || (s[0] >= 'A' && s[0] <= 'Z') ||
      (s[0] >= 'a' && s[0] <= 'z') || (s[0] == '"' &&
	((s[1] >= 'A' && s[1] <= 'Z') || (s[1] >= 'a' && s[1] <= 'z'))))){
    fputs(s, f2);
    rlen = get_line(&s, &blen, f1);
  }
  while(rlen != -1){
    i = rand() % 1000;
    if(i < k){
      fputs(s, f2);
      n++;
    }
    rlen = get_line(&s, &blen, f1);
    N++;
  }
  FCLOSE(f1);
  FCLOSE(f2);
  myfree(s);
  out_r(_("Done: selected %d out of %d rows.\n"), n, N);
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
}

/* Export current database as fixed width data file */
void exp_fwdf(){
  int i, j, k, *w;
  char *p, s[32], q[32], dfname[MLINE], cfname[MLINE];
  FILE *df, *cf;
  REAL r;
  if(ncol < 2){
    out_err(ERR, ERR_FILE, ERR_LINE,
	_("The current data file has less than 2 columns!"));
    return;
  }
  for(i = 1; i < ncol; i++)
    if(nn[0] != nn[i]){
      out_err(ERR, ERR_FILE, ERR_LINE,
	  _("There are columns with different number of rows!"));
      return;
    }

  /* Calculating the necessary width for each column */
  w = (int*)m_calloc(ncol, sizeof(int));
  for(i = 0; i < ncol; i++){
    acol[0] = i;
    alloc_cols(1, FALSE);
    for(j = 0; j < vn[i]; j++){
      r = xx[i][j];
      if(r == floor(r) && r >= -9999999999999999.0 && r <= 9999999999999999.0)
	snprintf(s, 32, "%16.16g", r);
      else
	snprintf(s, 32, "%16.10g", r);
      p = s;
      while(p[0] == ' ')
	p++;
      k = strlen(p);
      if(k > w[i])
	w[i] = k;
    }
  }
  
  out_i(_("Please enter name of the export file: ") );
  GETRLINE;
  sscanf(line, "%s", dfname);
  out_i(_("Please enter name of the list of columns file: ") );
  GETRLINE;
  sscanf(line, "%s", cfname);
  FOPEN(dfname, "wt", df);
  FOPEN(cfname, "wt", cf);
  j = 1;
  k = 0;

  /* saving the list of columns */
  for(i = 0; i < ncol; i++){
    k += w[i];
    fprintf(cf, "%s %i-%i\n", alias[i], j, k);
    j += w[i];
  }
  FCLOSE(cf);
  out_d(_("File \"%s\" saved!"), cfname);
  out_d("\n");
  
  /* saving the fixed width datafile */
#ifndef NO_GETTEXT
  SET_C_LOCALE;
#endif
  k = sizeof(REAL);
  for(i = 0; i < nn[0]; i++){
    for(j = 0; j < ncol; j++){
      FREAD(&r, k, 1, tmpptr[j]);
      if(r == SYSMIS){
	sprintf(s, "%%%is", w[j]);
        fprintf(df, s, " ");
      } else{
	if(r == floor(r) && r >= -9999999999999999.0 && r <= 9999999999999999.0)
	  snprintf(s, 32, "%16.16g", r);
	else
	  snprintf(s, 32, "%16.10g", r);
	p = s;
	while(p[0] == ' ')
	  p++;
	snprintf(q, 32, "%%%is", w[j]);
	fprintf(df, q, p);
      }
    }
    fprintf(df, "\n");
  }
#ifndef NO_GETTEXT
  RESET_LOCALE;
#endif
  
  /* Finishing */
  FCLOSE(df);
  out_r(_("File \"%s\" saved!"), dfname);
  out_r("\n\n");
}


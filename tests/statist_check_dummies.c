/* Dummy functions which helps us to more easily test some statist's functions
 * seperately 
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <check.h>


void finish() {
	/* no need to close tempfiles */
	/* no need to close pipes */
}


void out_err(int errno, char *modulname, int lno, char *fmt, ...) {
  /* va_list argptr; */

  fprintf(stderr,"Got a statist out_err() call.\n");
  exit(1);
}

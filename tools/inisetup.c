/* This file is part of statist
**
** It is distributed under the GNU General Public License.
** See the file COPYING for details.
**
**  $Id: inisetup.c,v 1.1 2005/10/29 17:39:23 jakson Exp $
***************************************************************/

#include<stdio.h>

/* This program creates a .bat file that popup a DOS window with the correct
 * %PATH% variable to run statist. */

int main(int argc, char **argv){
  char s[300];
  int i;
  strcpy(s, argv[0]);
  i = strlen(s) - 1;
  while(s[i] != '\\')
    i--;
  s[i] = 0;
  FILE *F;
  F = fopen("prompt.bat", "w");
  fprintf(F, "set PATH=\"%s;%%PATH%%\"\ncd \\\n@echo Type EXIT to close this window.\nCOMMAND.COM\n@cls\n@exit\n", s);
  fclose(F);
  return 0;
}

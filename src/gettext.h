/* This file is part of statist
**
** It is distributed under the GNU General Public License v>=2.
** See the file COPYING for details.
**
** (c) 1997 Dirk Melcher
** (c) Bernhard Reiter  http://www.usf.Uni-Osnabrueck.DE/~breiter
**
***************************************************************/

/* gettext.h for STATIST */

#ifndef STATIST_GETTEXT_HEADER
#define STATIST_GETTEXT_HEADER

#ifdef NO_GETTEXT
	#define _(english) (english)
#else
	#include <libintl.h>
	#include <locale.h>
	#define _(string) gettext (string)
#endif /* NO_GETTEXT */

#define N_(english) (english)


#endif /* STATIST_GETTEXT_HEADER */

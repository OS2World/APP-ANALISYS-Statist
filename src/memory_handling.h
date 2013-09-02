/* This file is part of statist
**
** It is Free Software distributed under the GNU General Public License.
** See the file COPYING for details.
** 
** (c) 1997 Dirk Melcher
** (c) 2001 Bernhard Reiter
**
** published by Bernhard Reiter  http://www.usf.Uni-Osnabrueck.DE/~breiter
** $Id: memory_handling.h,v 1.2 2005/08/16 08:07:44 jakson Exp $
***************************************************************/

#ifndef MEMORY_HANDLING_H
#define MEMORY_HANDLING_H

/* wrappers for malloc/calloc/realloc/free combination */
void *myrealloc(void *ptr, int newsize);
void *mycalloc(int nitems, int size);
void *mymalloc(int size);
void myfree(void *ptr);

/* keep track of memory allocated with m_calloc(), m_freeall frees all */
void *m_calloc(int nitems, int size);
int m_freeall();

#endif

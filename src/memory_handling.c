/* This file is part of statist
**
** It is distributed under the GNU General Public License.
** See the file COPYING for details.
**
** (c) 1997 Dirk Melcher old email: Dirk.Melcher@usf.Uni-Osnabrueck.DE
** (c) 2001 Bernhard Reiter
**
** published by Bernhard Reiter  http://www.usf.Uni-Osnabrueck.DE/~breiter
** $Id: memory_handling.c,v 1.3 2005/08/31 12:48:16 jakson Exp $
***************************************************************/


/* Memory handling functions
 *
 * */

#include <stdlib.h>
#include "statist.h"

#include "memory_handling.h"

typedef struct {
	size_t size;
	void * previous_chunk;
	void * memory;
} memory_chunk ;

/* pointer to the last memory chunk */
static memory_chunk *ptr_memory_chunks=NULL;



void *myrealloc(void * pointer, int newsize) {
  /* realloc() with error check */
  pointer = realloc(pointer, newsize);
  if (pointer == NULL) {
     out_err(FAT, ERR_FILE, ERR_LINE,
	 _("Not enough memory. Terminating!") );
     finish();
     exit(1);
  }
  return pointer;
}

void *mycalloc(int nitems, int size) {
  /* calloc() with error check */
  void *pointer;
  pointer = calloc(nitems, size);
  if (pointer == NULL) {
     out_err(FAT, ERR_FILE, ERR_LINE,
	 _("Not enough memory. Terminating!") );
     finish();
     exit(1);
  }
  return pointer;
}


void *mymalloc(int size) {
  /* malloc() with error check */
  void *pointer;
  pointer = malloc(size);
  if (pointer == NULL) {
     out_err(FAT, ERR_FILE, ERR_LINE,
	 _("Not enough memory. Terminating!") );
     finish();
     exit(1);
  }
  return pointer;
}

void myfree(void *ptr) {
  /* pair function to mycalloc() and mymalloc() */
  free(ptr);
}


void *m_calloc(int nitems, int size) {
  /* calloc() should be used for temp variables,
   * keeps a record of allocated blocks
   * you need to call m_freeall to free them all
   */
  memory_chunk *chunk;
  chunk = mycalloc(1,sizeof(memory_chunk));

  /* get mem and fill struct */
  chunk->memory = mycalloc(nitems, size);
  chunk->size = nitems * size;
  chunk->previous_chunk=ptr_memory_chunks;
  /* put entry pointer to last chunk */
  ptr_memory_chunks=chunk;

  return chunk->memory;
}


int m_freeall(void) {
  /* frees all blocks allocated by several calls to m_calloc() */
  memory_chunk *chunk=NULL;

  while(ptr_memory_chunks!=NULL) {
    chunk=ptr_memory_chunks;

    myfree(chunk->memory);
    ptr_memory_chunks=chunk->previous_chunk;

    /* free struct */
    myfree(chunk);
  }

  return 0;
}


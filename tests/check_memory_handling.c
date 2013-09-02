#include <stdlib.h>
#include <memory_handling.h>

#include <check.h>


START_TEST(test_small_malloc)
{
	void * mem=NULL;
	mem=mymalloc(1000);

	fail_unless(mem!=NULL,"Got no memory!");

	myfree(mem);
}
END_TEST


START_TEST(test_small_calloc)
{
	void * mem=NULL;
	int i=0;

	mem=mycalloc(100,100);

	fail_unless(mem!=NULL,"Got no memory!");

	for(i=0; i<100*100;i++) {
		if (((char*) mem)[i] != 0)  
			fail("Memory was not clear!");
	}

	myfree(mem);
}
END_TEST

START_TEST(test_small_m_calloc)
{
	void * mem=NULL;
	int i=0;

	mem=m_calloc(100,100);

	fail_unless(mem!=NULL,"Got no memory!");

	for(i=0; i<100*100;i++) {
		if (((char*) mem)[i] != 0)  
			fail("Memory was not clear!");
	}

	m_freeall();
}
END_TEST


START_TEST(test_mult_m_calloc)
{
	void **memarray=NULL;
#define CHUNKS (10) 
#define CHUNKSIZE (1000)
	int i=0;
	int j=0;

	memarray=m_calloc(CHUNKS,sizeof(void *));

	fail_unless(memarray!=NULL,"Got no memory!");

	for(i=0; i<CHUNKS;i++) {
		memarray[i]=m_calloc(1, CHUNKSIZE);

		for(j=0; j<CHUNKSIZE;j++) {
			if (((char*) memarray[i])[j] != 0)  
				fail("Memory was not clear!");
		}
	}

	m_freeall();
}
END_TEST

/******************************************************************************
 * Test Suites
 */

Suite *memory_suite (void) 
{ 
	 Suite *s = suite_create("Memory Handling"); 
	 TCase *tc_core = tcase_create("Simple");
	 TCase *tc_m_calloc = tcase_create("M_calloc()");
	   
	 suite_add_tcase (s, tc_core);
	 suite_add_tcase (s, tc_m_calloc);
	    
	 tcase_add_test (tc_core, test_small_malloc); 
	 tcase_add_test (tc_core, test_small_calloc); 

	 tcase_add_test (tc_m_calloc, test_small_m_calloc );
	 tcase_add_test (tc_m_calloc, test_mult_m_calloc );
	 return s; 
}

/******************************************************************************
 * Main 
 */

int main (void) 
/* standard code from the check-0.9.2 tutorial			*/
{ 
	int nf; 

	Suite *s = memory_suite(); 

	SRunner *sr = srunner_create(s); 
	srunner_run_all (sr, CK_NORMAL);
	/* srunner_run_all (sr, CK_VERBOSE); */
	nf = srunner_ntests_failed(sr); 
	srunner_free(sr); 
	return (nf == 0) ? EXIT_SUCCESS : EXIT_FAILURE; 
}



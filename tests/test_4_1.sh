#!/bin/sh
# integration test
# testing menu item 4-1
# using the --bernhard option the table is language independent

expected1="
#Result general statistical information in a table
#n	mean	m-conf	m+conf	median	me_c_lo	me_c_up	quar_lo	quar_up	sdv	varc(%)	sdv_err	min	max
9	5	2.89512	7.10488	5	2	8	2.5	7.5	2.73861	54.772256	0.912871	1	9	"


menu="4
1
"

exitprogram="0
0
"


parameter1="a



"
string1="${menu}${parameter1}${exitprogram}"


actual=`echo "${string1}" | ./run_statist.sh --bernhard --silent t1.dat | sed -e 1d`

    if [ x"${expected1}" != x"${actual}" ]; then
        echo "Problem with menu 4-1" ;
        echo "Expected:";
        echo "${expected1}";
        echo "Got:";
        echo "${actual}";
        exit 1;
    fi
    


exit 0

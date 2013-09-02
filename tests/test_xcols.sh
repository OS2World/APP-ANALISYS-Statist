#!/bin/sh

# Based on test_4_1.sh written by Bernhard Reiter.

# With this script we will read a fixed width datafile with information about
# people, identify the couples and calculate both the mean age difference and
# the mean income difference between husbands and wives. The entire process was
# made by hand before and the results must coincide.

# Extracting columns from fixed width datafile
./run_statist.sh --silent --xcols xcols.config xcols.dat people.dat


# Creating key variable for merging the data files later
awk '{if(/state/){print "#%key" "\t" $0} else
       {
         if(NF == 0) {print $0}
         else{
          {key = $1$2$3} {print key "\t" $0}
         }
       }
}' people.dat > people_with_key.dat


# Creating the two data files that will be merged
awk '{if(/state/ || NF == 0){print $0} else
        {if($5 < 3 && $6 == 0) {print $0}}
}' people_with_key.dat > men.dat

awk '{if(/state/ || NF == 0){print $0} else
        {if($5 < 3 && $6 == 1) {print $0}}
}' people_with_key.dat > women.dat


# Merging the two data files
join -e "" women.dat men.dat > couples1.dat

# Making the file more human readable
sed 's/\ /\t/g
s/#%state/state/g' couples1.dat > couples2.dat

# Creating variables "Age difference between husband and wife" and "Income
# difference between husband and wife"
awk '{if(/state/){print $0 "\tage_d\tincome_d"} else
       {
         if(NF == 0) {print $0}
         else{
          {ad = $14 - $7} {id = $15 - $8} {print $0 "\t" ad "\t" id}
         }
       }
}' couples2.dat > couples.dat


expected="
#Result general statistical information in a table
#n	mean	m-conf	m+conf	median	me_c_lo	me_c_up	quar_lo	quar_up	sdv	varc(%)	sdv_err	min	max
7	103.286	-56.6998	263.271	100	-111	450	-56	284.5	173.093	167.586728	65.4231	-111	450	"


menu="4
1
income_d



0
0
"
actual=`echo "${menu}" | ./run_statist.sh --silent --bernhard couples.dat | sed -e 1d`

    if [ x"${expected}" != x"${actual}" ]; then
        echo "Problem with --xcols" ;
        echo "Expected:";
        echo "${expected}";
        echo "Got:";
        echo "${actual}";
        exit 1;
    fi
    
rm -f couples*.dat men.dat people.dat people_with_key.dat women.dat

exit 0


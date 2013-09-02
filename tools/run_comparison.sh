#!/bin/bash

# This file is part of statist
#
# It is distributed under the GNU General Public License.
# See the file COPYING for details.
#
# (C) 2006 Jakson Alves de Aquino <jalvesaq@gmail.com>
#
#  $Id: run_comparison.sh,v 1.3 2006/09/10 02:42:12 jakson Exp $

# Warning: This script is the worst possible example of how to do statistics.
# The analyzes carried here have the only goal of comparing two different
# versions of statist.

# This script will run two different versions of statist to check whether the
# new version still produces the correct results. The old version is the one
# tagged as "bernhard-last-maint-2005" (here called $statist_102) and the new
# one is the 1.4.1.  You have to compile these versions of statist and write
# the path to them in the variables below. The path to statist "examples"
# subdirectory is necessary because the file city.csv will be used.

statist_102="./statist-1.0.2"
statist_141="../src/statist"
examples_dir="../examples"

LANG=C
LANGUAGE=C
LC_ALL=C

# Before running statist-1.0.2 we must create four new data files to run
# analysis that require the same number of rows. Statist-1.0.2 has the option
# -delrow, but, alone, it isn't enough.

# NOTE 1: We can't use statist-1.0.2 to "export a file" with the columns
# because the rows with missing values will be dislocated. Example:

# 34 44                     34 44
# 35 M       will become    35 53
# 36 53                     36 .

# As you can see, the "53" is no longer aligned with "36".

# NOTE 2: statist-1.0.2 would be unable to read the two columns that it had
# just exported due to the "dot" at the third column. The correct (for statist)
# would be to put a "M" in the second line, not a dot at the third one. So,
# let's use awk:

sed -e 's/age/#%age/' $examples_dir/city.csv > city.dat
awk '{if($2 == "sex") {print "#%sex god"} else if($2 == 1 || $2 == 0) print $2, $6}' city.dat > zero_one.dat
awk '{print $1, $7}' city.dat > two_columns.dat
awk '{print $1, $7, $8}' city.dat > three_columns.dat
awk '{print $1, $2, $4, $7, $8}' city.dat > five_columns.dat


# Now we'll run statist-1.0.2 six times.

# FIRST: We'll use two_columns.dat to run all menu items that require two
# columns with exactly the same number of rows. Note: In a real research we
# wouldn't have this privilege. We should use awk to save all combinations of
# two columns that we were interested in.

# SECOND: We'll use three_columns.dat to run all menu items that require
# three columns with exactly the same number of rows.

# THIRD: We'll run partial linear correlation with 5 items.

# FORTH: We'll run Chi^2-fourfold-test (zero_one.dat).

# FIFTH: We'll run all remaining menu items.

# SIXTH: We'll run the probit analysis.

# menu choices for the first run:
echo '2
1
age
deg
2
age
deg
5
age
deg
5
0
3
2
age
deg
8
age
deg
0
0
' > run_it

$statist_102 -delrow -silent -noplot two_columns.dat < run_it > st102a


# menu choices for the second run (three_columns):
echo '2
3
3
inc
deg
age
9
3
inc
deg
age
10
3
inc
deg
age
5
0
3
5
3
inc
deg
age
0
0
' > run_it

$statist_102 -delrow -silent -noplot three_columns.dat < run_it > st102b


# menu choices for the third run (five_columns):
echo '2
4
5
inc
age
deg
hap
sex
6
5
inc
age
deg
hap
sex
7
5
inc
age
deg
hap
sex
0
0
' > run_it

$statist_102 -silent -delrow -noplot five_columns.dat < run_it > st102c


# menu choices for the fourth run (Chi^2):
echo '3
4
sex
god
0
0
' > run_it

$statist_102 -silent -delrow -noplot zero_one.dat < run_it > st102d


# menu choices for the fifth run (city.dat):
echo '2
8
sex
age
0
3
1
deg
inc
6
deg
inc
7
3
deg
inc
age
9
deg
10
sex
god
0
4
1
age



3
inc
n
4
age
0
0' > run_it

$statist_102 -silent -noplot city.dat < run_it > st102e

# menu choices for the sixth run (probit):
echo '4
2
y
dose
ef
2
N
dos
num
ef
0
0
' > run_it

$statist_102 -silent -noplot $examples_dir/probit.dat < run_it > st102f

cat st102* > result_statist-1.0.2
rm st102*


##########################################################################
# Now, finally, we'll write the menu choices for statist-1.4.1 and run it.
# Note that statist-1.4.1 doesn't have (and doesn't need) the option -delrow.

echo '2
1
age
deg
2
age
deg
5
age
deg
5
0
3
2
age
deg
8
age
deg
0
2
3
inc
deg
age

9
inc
deg
age

10
inc
deg
age

5
0
3
5
inc
deg
age

y
0
2
4
inc
age
deg
hap
sex
6
inc
age
deg
hap
sex

7
inc
age
deg
hap
sex

0
3
4
sex
god
0
2
8
sex
age
0
3
1
deg
inc
6
deg
inc
7
deg
inc
age

9
deg
10
sex
god
0
4
1
age



3
inc
n
4
age
0
0' > run_it

$statist_141 --na-string "M" --noplot --silent city.dat < run_it > st140a


# menu choices for probit:
echo '4
2
y
dose
ef
2

dos
num
ef
0
0
' > run_it


$statist_141 --noplot --silent $examples_dir/probit.dat < run_it > st140b

cat st140* > result_statist-1.4.1

rm st140* city.dat run_it five_columns.dat three_columns.dat two_columns.dat zero_one.dat

# Note: It's expected different results for "Partial linear correlation" with 5
# variables because statist-1.0.2 was using an unintialized variable. This bug
# was fixed in statist-1.3.1 (thanks to valgrind).

echo
echo "The results were saved in the files"
echo "result_statist-1.0.2 and result_statist-1.4.1."
echo
exit 0



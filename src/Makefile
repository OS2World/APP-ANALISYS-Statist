# This file is part of statist
# 
# It is distributed under the GNU General Public License.
# See the file COPYING for details.
# 
# Copyright (C) 1997-1999 Dirk Melcher <Dirk.Melcher@usf.Uni-Osnabrueck.DE>
# Copyright (C) 1998-2006 Bernhard Reiter <bernhard@intevation.de>
# Copyright (C) 2005-2006 Jakson Aquino <jalvesaq@gmail.com>
#
# $Id: Makefile,v 1.24 2006/11/18 19:47:23 jakson Exp $
#*************************************************************/
# Makefile fuer UNIX-Version und gnuplot-Graphik
# Makefile adjusted for unix with gnuplot graphics.

include ../Makefile.inc

# Choose your compiler, the example flags are for the gcc
CC = gcc

LIBTOOL = ar crus

COMPILERCFLAGS = -c -Wall -D__EMX__ -DOS2  -D__ST_MT_ERRNO__ 
COMPILERCFLAGS += -g
COMPILERCFLAGS += -O2
#COMPILERCFLAGS += -DDEBUG

LIBS += -lm
LFLAGS = -Zexe -Zomf -Zmap -Zargs-wild -o
#LFLAGS = -g -o


# Only use the next line, if you want to create a lib ("shlib" or "lib")
#STATISTX = -DSTATIST_X

###############################################################################
# ADDITIONAL OPTIONS FOR MS-DOS/WINDOWS:
#
#COMPILERCFLAGS += -DMSDOS
#
# To compile with MinGW (or Dev-C++), uncomment the options -DMSDOS,
# -I"C:/DEV-CPP/include", and  -L"C:/DEV-CPP/lib".
#
# If you will compile with gettext support, using MinGW (or Dev-C++) you have
# to indicate where libintl.h, libintl3.dll and libiconv.dll are.  So,
# uncomment the options -I"C:/GnuWin32/include and -I"C:/GnuWin32/bin; You
# have also to uncomment the option -lintl3. With CygWin the only option you
# have to set is -lintl.
#
# You need a version of gnuplot that comes with pgnuplot.exe. If not, uncomment
# the option -DNOPIPE.
#
#COMPILERCFLAGS += -DNOPIPE
#COMPILERCFLAGS += -I"C:/DEV-CPP/include"
#COMPILERCFLAGS += -I"C:/GnuWin32/include"
#
#LIBS += -L"C:/DEV-CPP/lib"
#LIBS += -L"C:/GnuWin32/bin"
#LIBS += -lintl3
#LIBS += -lintl
#
# If your compiler don't find "values.h" or if it don't find the definition of
# DBL_MAX, uncomment the option -DNO_VALUES_H.
#
# If you don't want to use gettext, you have to define NO_GETTEXT, and Statist
# will use its own macro. If you want real gettext, use one of the following
# instead:
#
#EXTRACFLAGS = -DNO_GETTEXT
EXTRACFLAGS += -DNO_VALUES_H
EXTRACFLAGS += -DNO_IOCTIL_H
#
# END OF OPTIONS FOR MS-DOS/WINDOWS
###############################################################################


# Flags for statist.c:
STATISTCFLAGS = -DLOCALEDIR="\"$(MESSAGE_CATALOGS)\""
STATISTCFLAGS += -DDOCDIR="\"$(DOCDIR)\"" 

CFLAGS= $(COMPILERCFLAGS) $(EXTRACFLAGS) $(STATISTX)

#MAIN = $(HOME)/tools/statist
MAIN=statist.exe
SH_LIB = libstatist.so.1.3.1
STAT_LIB = libstatist.a


OBJECTS = statist.o menue.o procs.o funcs.o data.o plot.o \
memory_handling.o
LIB_OBJECTS = menue.o procs.o funcs.o data.o plot.o memory_handling.o

.POSIX:

# Decide between "main", "shlib" and "lib"
all: main
#all: lib
#all: shlib

main:   $(MAIN)

check: main
	echo Running tests... in the appropriate directory
	(cd ../tests ; \
	$(MAKE) check \
	)

shlib:	$(SH_LIB)

lib:	$(STAT_LIB)


$(MAIN) : $(OBJECTS) Makefile
	$(CC) $(LFLAGS) $(MAIN) $(OBJECTS) $(LIBS)

# change AB:
$(SH_LIB):		$(LIB_OBJECTS) Makefile
	$(CC) -shared $(LIB_OBJECTS) -o $(SH_LIB)

$(STAT_LIB):	$(LIB_OBJECTS) Makefile
	 $(LIBTOOL) $(STAT_LIB) $(LIB_OBJECTS)

statist.o : statist.c statist.h menue.h data.h funcs.h
	$(CC) $(CFLAGS) $(STATISTCFLAGS) statist.c

menue.o : menue.c menue.h procs.h data.h statist.h
	$(CC) $(CFLAGS) menue.c

procs.o : procs.c procs.h funcs.h statist.h
	$(CC) $(CFLAGS) procs.c

funcs.o : funcs.c funcs.h statist.h
	$(CC) $(CFLAGS) funcs.c

data.o : data.c data.h statist.h
	$(CC) $(CFLAGS) data.c

plot.o : plot.c plot.h statist.h funcs.h
	$(CC) $(CFLAGS) plot.c

memory_handling.o : memory_handling.c memory_handling.h statist.h
	$(CC) $(CFLAGS) memory_handling.c

install: $(MAIN)
	install -d $(PREFIX)/bin
	install -s $(MAIN) $(PREFIX)/bin/
	(cd ../po; make install)
	(cd ../doc; make install)

uninstall:
	rm -f $(PREFIX)/bin/statist
	(cd ../po; make uninstall)
	(cd ../doc; make uninstall)

clean:
	-rm -f $(OBJECTS) $(MAIN) *.log *.aux core $(STAT_LIB) $(SH_LIB)
	(cd ../tests ; \
	$(MAKE) clean \
	)

include Makefile.inc

.POSIX:

all:
	(cd src; make)
	(cd po; make)
	(cd tests; make)
	(cd doc; make)

install:
	(cd src; make install)

uninstall:
	(cd src; make uninstall)

check:
	(cd src; make)
	(cd tests; make)

clean:
	(cd src; make clean)
	(cd doc; make clean)
	(cd po; make clean)


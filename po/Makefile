# This file is part of statist
# 
# It is Free Software distributed under the GNU General Public License v>=2.
# See the file COPYING for details.
# 
# Author: 2005 Jakson Alves de Aquino <jakson.alves at ig.com.br>, 2005

include ../Makefile.inc

POSTFIX = /LC_MESSAGES
POSTFIXMO = $(POSTFIX)/statist.mo

SOURCES = ../src/data.c ../src/menue.c ../src/memory_handling.c \
  ../src/funcs.c ../src/plot.c ../src/procs.c \
  ../src/statist.c ../src/data.h  ../src/funcs.h ../src/gettext.h \
  ../src/memory_handling.h ../src/menue.h ../src/plot.h ../src/procs.h \
  ../src/statist.h
POFILES = pt_BR.po es_ES.po de_DE.po fr_FR.po it_IT.po
MOFILES = pt_BR.gmo es_ES.gmo de_DE.gmo fr_FR.gmo it_IT.gmo

.POSIX:

all: gmo

help :
	@echo "Type:"
	@echo "make po        : to update the .po files."
	@echo "make gmo       : to create the .gmo files."
	@echo "make noaccents : to copy the .po files into new .po files where"
	@echo "                 all words are without accents."
	@echo "make install   : to copy  the .gmo files to"
	@echo "                 $(MESSAGE_CATALOGS)/xx_XX/$(POSTFIX)/"

statist.pot: $(SOURCES)
	xgettext -k_ -kN_ -d statist -o statist.pot -s $(SOURCES)

%.po: statist.pot
	msgmerge -F -U $@ statist.pot
	touch $@

%.gmo: %.po
	msgfmt -c -v -o $@ $<

po: $(POFILES)

gmo: $(MOFILES)

noaccents:
	@echo "I'm sorry! This option wasn't implemented yet."

install: $(MOFILES)
	install -d $(MESSAGE_CATALOGS)/pt_BR$(POSTFIX)
	install -d $(MESSAGE_CATALOGS)/de$(POSTFIX)
	install -d $(MESSAGE_CATALOGS)/es$(POSTFIX)
	install -d $(MESSAGE_CATALOGS)/it$(POSTFIX)
	install -d $(MESSAGE_CATALOGS)/fr$(POSTFIX)
	cp pt_BR.gmo $(MESSAGE_CATALOGS)/pt_BR$(POSTFIXMO)
	cp de_DE.gmo $(MESSAGE_CATALOGS)/de$(POSTFIXMO)
	cp es_ES.gmo $(MESSAGE_CATALOGS)/es$(POSTFIXMO)
	cp it_IT.gmo $(MESSAGE_CATALOGS)/it$(POSTFIXMO)
	cp fr_FR.gmo $(MESSAGE_CATALOGS)/fr$(POSTFIXMO)

uninstall:
	rm -f $(MESSAGE_CATALOGS)/pt_BR$(POSTFIXMO)
	rm -f $(MESSAGE_CATALOGS)/de$(POSTFIXMO)
	rm -f $(MESSAGE_CATALOGS)/es$(POSTFIXMO)
	rm -f $(MESSAGE_CATALOGS)/it$(POSTFIXMO)
	rm -f $(MESSAGE_CATALOGS)/fr$(POSTFIXMO)

clean:
	-rm *.po~ *.gmo


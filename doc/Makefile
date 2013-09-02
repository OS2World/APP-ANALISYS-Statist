# This file is part of statist
# 
# It is Free Software distributed under the GNU General Public License v>=2.
# See the file COPYING for details.
# 
# Author: 2005 Jakson Alves de Aquino <jakson.alves at ig.com.br>, 2005

include ../Makefile.inc

MANDIR = $(PREFIX)/share/man/man1

.PHONY : all pdf install halfclean clean

all:  pdf


pdf: programm.pdf stat_man.pdf manual-en.pdf manual-pt.pdf

boxplot-en.pdf: boxplot-en.fig
	fig2dev -L pdf boxplot-en.fig boxplot-en.pdf

boxplot-pt.pdf: boxplot-pt.fig
	fig2dev -L pdf boxplot-pt.fig boxplot-pt.pdf

stat_man.pdf: stat_man.tex
	latex stat_man
	bibtex stat_man
	latex stat_man
	latex stat_man
	dvipdfm stat_man

programm.pdf: programm.tex
	pdflatex programm
	pdflatex programm

manual-en.pdf: manual-en.tex boxplot-en.pdf
	pdflatex manual-en
	pdflatex manual-en

manual-pt.pdf: manual-pt.tex boxplot-pt.pdf
	pdflatex manual-pt
	pdflatex manual-pt

install: stat_man.pdf manual-en.pdf manual-pt.pdf
	install -d $(DOCDIR)/examples
	install -d $(MANDIR)
	cp -p statistrc $(DOCDIR)/sample.statistrc
	cp -p stat_man.pdf $(DOCDIR)
	cp -p manual-en.pdf $(DOCDIR)
	cp -p manual-pt.pdf $(DOCDIR)
	cp -p statist.1 $(MANDIR)
	cp -p ../CHANGES $(DOCDIR)
	cp -p ../COPYING $(DOCDIR)
	cp -p ../KNOWN_BUGS $(DOCDIR)
	cp -p ../README $(DOCDIR)
	cp -p ../CREDITS $(DOCDIR)
	(cd ../examples; cp -p *.dat *.labels README $(DOCDIR)/examples)

uninstall:
	-rm $(DOCDIR)/sample.statistrc
	-rm $(DOCDIR)/stat_man.pdf
	-rm $(DOCDIR)/manual-en.pdf
	-rm $(DOCDIR)/manual-pt.pdf
	-rm $(MANDIR)/statist.1
	-rm $(DOCDIR)/CHANGES
	-rm $(DOCDIR)/COPYING
	-rm $(DOCDIR)/KNOWN_BUGS
	-rm $(DOCDIR)/README
	-rm $(DOCDIR)/CREDITS
	-rm $(DOCDIR)/examples/*.dat
	-rm $(DOCDIR)/examples/*.labels
	-rm $(DOCDIR)/examples/README
	-rmdir $(DOCDIR)/examples
	-rmdir $(DOCDIR)

clean:
	-rm *.log *.toc *.aux *.bbl *.blg *.dvi *.eps

reallyclean:
	-rm *.log *.toc *.aux *.bbl *.blg *.dvi *.eps *.pdf

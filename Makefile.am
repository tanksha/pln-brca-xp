include guile.am
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

SOURCES = \
    pln-bio/bio-utils.scm \
    pln-bio/ure-utils.scm \
    pln-bio/preprocess.scm \
    pln-bio/intensional-similarity.scm


abstopdir=$(abspath $(lastword $(srcdir)))
rulesdir = $(abstopdir)/rules/
kbsdir = $(abstopdir)/kbs/


pln-bio/bio-utils.scm: Makefile pln-bio/bio-utils.scm.in
	mkdir -p pln-bio
	sed -e "s|@kbsdir\@|\"$(kbsdir)\"|; s|@rulesdir\@|\"$(rulesdir)\"|" \
		$(srcdir)/pln-bio/bio-utils.scm.in > pln-bio/bio-utils.scm

clean-go:
	-$(RM) $(GOBJECTS)
.PHONY: clean-go

CLEANFILES += \
	$(GOBJECTS)	\
	*.log

EXTRA_DIST = \
	$(SOURCES)

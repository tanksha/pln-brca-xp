include build-aux/guile.am
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

SOURCES = \
    pln-bio/bio-utils.scm \
    pln-bio/ure-utils.scm \
	pln-bio/rule-utils.scm \
    pln-bio/preprocess.scm \
	pln-bio/expr.scm \
	pln-bio/combo-preprocess.scm \
	pln-bio/bp-deduction.scm \
    pln-bio/similarity-xp.scm \
	pln-bio/concept-deduction.scm


abstopdir=$(abspath $(lastword $(srcdir)))/


pln-bio/bio-utils.scm: Makefile pln-bio/bio-utils.scm.in
	mkdir -p pln-bio
	sed -e "s|@abstopdir\@|\"$(abstopdir)\"|" \
		$(srcdir)/pln-bio/bio-utils.scm.in > pln-bio/bio-utils.scm

clean-go:
	-$(RM) $(GOBJECTS)
.PHONY: clean-go

CLEANFILES += \
	$(GOBJECTS)	\
	*.log

EXTRA_DIST = \
	$(SOURCES)

include guile.am
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

SOURCES = \
    pln-bio/bio-utils.scm \
    pln-bio/ure-utils.scm \
    pln-bio/preprocess.scm \
    pln-bio/main/go-intensional-similarity.scm 

clean-go:
	-$(RM) $(GOBJECTS)
.PHONY: clean-go

CLEANFILES += \
	$(GOBJECTS)	\
	*.log

EXTRA_DIST = \
	$(SOURCES)

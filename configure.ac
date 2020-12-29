AC_PREREQ(2.52)
AC_CONFIG_SRCDIR([env.in])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([-Wall])

#
# Check for Guile
#
GUILE_PKG([3.0 2.2])
GUILE_PROGS
GUILE_SITE_DIR
GUILE_FLAGS

##Check Guile modules

GUILE_MODULE_AVAILABLE(HAVE_OPENCOG, [(opencog) (opencog exec)])
if test "x$HAVE_OPENCOG" = "xno"; then
    AC_MSG_ERROR([The atomspace is required. See https://github.com/singnet/atomspace for details])
fi

GUILE_MODULE_AVAILABLE(HAVE_URE, [(opencog ure)])
if test "x$HAVE_URE" = "xno"; then
    AC_MSG_ERROR([The Unified Rule Engine (URE) is required. See https://github.com/singnet/ure for details])
fi

GUILE_MODULE_AVAILABLE(HAVE_PLN, [(opencog pln)])
if test "x$HAVE_PLN" = "xno"; then
    AC_MSG_ERROR([PLN is required. See https://github.com/singnet/pln for details])
fi

GUILE_MODULE_AVAILABLE(HAVE_BIOSCI, [(opencog bioscience)])
if test "x$HAVE_BIOSCI" = "xno"; then
	AC_MSG_ERROR([The (opencog bioscience) module is needed. See https://github.com/opencog/agi-bio for details.])
fi

AC_CONFIG_FILES([env], [chmod +x env])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT

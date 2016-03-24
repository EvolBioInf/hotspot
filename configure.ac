AC_INIT([hotspot], [0.4])
AM_INIT_AUTOMAKE([-Wall foreign ])

#AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_RANLIB
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_GREP

# Execute all tests using C
AC_LANG(C)

AC_CHECK_LIB(m, cos)

# Check for libdivsufsort
AC_CHECK_HEADERS([divsufsort.h],[have_libdivsufsort=yes],[have_libdivsufsort=no])
AC_CHECK_LIB(divsufsort, divsufsort, [], [have_libdivsufsort=no])

AS_IF([test "x$have_libdivsufsort" = "xno"],
	[AC_MSG_ERROR([Configuration for libdivsufsort failed. Please install libdivsufsort.])
])

# Check for tabix
AC_CHECK_PROG(TABIX_CHECK,tabix,yes)
AS_IF([test x"$TABIX_CHECK" != x"yes"],
    [AC_MSG_ERROR([Please install tabix.])
])

AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])


# Check for various headers including those used by libdivsufsort.
AC_CHECK_HEADERS([limits.h stdlib.h string.h unistd.h stdint.h inttypes.h err.h errno.h fcntl.h])
AC_CHECK_HEADERS([float.h])


# AC_C_INLINE
AC_TYPE_SIZE_T
# AC_TYPE_SSIZE_T
# AC_TYPE_INT32_T
# AC_TYPE_UINT8_T
# AC_HEADER_STDBOOL

AC_FUNC_MALLOC
AC_FUNC_REALLOC


AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([strchr strcspn strerror strstr])

AC_CONFIG_HEADERS([common/config.h:common/config.hin])

AC_CONFIG_FILES([
 Makefile
 common/Makefile
 src_xov/Makefile
 src_aso/Makefile
 src_asp/Makefile
 src_six/Makefile
])
AC_OUTPUT


bin_PROGRAMS= six
six_SOURCES= interface.c interface.h gsl_rng.c gsl_rng.h six.c
six_CPPFLAGS= -W -Wall -I$(top_srcdir)/common
six_LDADD= ../common/libcommon.a  -lm -lgsl -lblas

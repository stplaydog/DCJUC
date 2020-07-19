noinst_HEADERS=insdcj.h insknap.h insmed.h list.h listdis.h listmed.h timer.h insdis.h instance.h listdcj.h listknap.h plist.h utils.h insmc.h listmc.h 
bin_PROGRAMS=optkit
optkit_SOURCES=insdis.cpp insmed.cpp listdcj.cpp listmed.cpp plist.cpp timer.cpp insdcj.cpp insdis.h list.cpp listknap.cpp insknap.cpp instance.cpp listdis.cpp main.cpp utils.cpp insmc.cpp listmc.cpp 
optkit_LDADD=$(OPENMP_LIBi) ../lib/libds.la
#optkit_LIBADD = ../lib/libds.la
if DEBUG
AM_CXXFLAGS = -g -fopenmp -I$(top_srcdir)/lib -std=c++11
else
AM_CXXFLAGS = -g3 -Wall -fopenmp -static -I$(top_srcdir)/lib -std=c++11 -DUSE_DEBUG
endif

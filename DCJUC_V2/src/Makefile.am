noinst_HEADERS=insknap.h insmed.h list.h timer.h insdis.h instance.h plist.h utils.h insmc.h listint.h 
bin_PROGRAMS=optkit
optkit_SOURCES=insdis.cpp insmed.cpp plist.cpp timer.cpp insdis.h list.cpp insknap.cpp instance.cpp main.cpp utils.cpp insmc.cpp listint.cpp 
optkit_LDADD=$(OPENMP_LIBi) ../lib/libds.la
#optkit_LIBADD = ../lib/libds.la
if DEBUG
AM_CXXFLAGS = -g  -fopenmp -I$(top_srcdir)/lib -std=c++11
else
AM_CXXFLAGS = -g -fopenmp -I$(top_srcdir)/lib -std=c++11 -DUSE_DEBUG -gdwarf-2
endif

include ../Makefile.common

median.a: local.o
	$(AR) rcs ../lib/median.a local.o

local.o: local.c 
	$(CC) $(FLAGS) $(LIB) $(INC) -c local.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/median.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/median.a

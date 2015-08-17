include ../Makefile.common

optkit.a: plist.o
	$(AR) rcs ../lib/optkit.a plist.o

plist.o: plist.c 
	$(CC) $(FLAGS) $(LIB) $(INC) -c plist.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/optkit.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/optkit.a

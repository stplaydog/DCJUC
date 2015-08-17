include ../Makefile.common

dist.a: dist.o
	$(AR) rcs ../lib/dist.a dist.o

dist.o: dist.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c dist.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/dist.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/dist.a

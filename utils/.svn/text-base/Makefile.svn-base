include ../Makefile.common

utils.a: graph.o list.o bitarr.o spectral.o  newick.o
	$(AR) rcs ../lib/utils.a graph.o list.o bitarr.o spectral.o newick.o

graph.o: graph.c 
	$(CC) $(FLAGS) $(LIB) $(INC) -c graph.c 
list.o: list.c 
	$(CC) $(FLAGS) $(LIB) $(INC) -c list.c 
bitarr.o: bitarr.c 
	$(CC) $(FLAGS) $(LIB) $(INC) -c bitarr.c 
spectral.o: spectral.c 
	$(CC) $(FLAGS) $(LIB) $(INC) -c spectral.c 
newick.o: newick.c 
	$(CC) $(FLAGS) $(LIB) $(INC) -c newick.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/utils.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/utils.a

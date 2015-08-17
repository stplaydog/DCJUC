include ../Makefile.common

phylogeny.a: nj.o tree.o super.o model.o
	$(AR) rcs ../lib/phylogeny.a nj.o tree.o super.o model.o

nj.o: nj.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c nj.c 
tree.o: tree.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c tree.c 
super.o: super.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c super.c 
model.o: model.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c model.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/phylogeny.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/phylogeny.a

include ../Makefile.common

DCJUC: main.o ../lib/dist.a ../lib/median.a ../lib/phylogeny.a ../lib/generator.a ../lib/utils.a
	$(CC) $(FLAGS) $(LIB) $(INC) -o DCJUC main.o ../lib/dist.a ../lib/median.a ../lib/phylogeny.a ../lib/generator.a ../lib/utils.a $(LINK)

main.o: main.c
	$(CC) $(FLAGS) $(LIB) $(INC) -c main.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/main.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/main.a

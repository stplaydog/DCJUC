include ../Makefile.common

generator.a: sim.o
	$(AR) rcs ../lib/generator.a sim.o

sim.o: sim.c 
	$(CC) $(FLAGS) $(LIB) $(INC) -c sim.c 

all:
	make

clean:
	$(RM) *.o
	$(RM) ../lib/generator.a

clean_all:
	$(RM) *.o
	$(RM) ../lib/generator.a

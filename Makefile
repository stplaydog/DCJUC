include Makefile.common

default:
	cd dist && make
	cd generator && make
	cd median && make
	cd phylogeny && make
	cd utils && make
	cd optkit && make
	cd main && make

all:
	cd dist && make all
	cd generator && make all
	cd median && make all
	cd phylogeny && make all
	cd utils && make all
	cd optkit && make all
	cd main && make all

clean:
	cd dist && make clean
	cd generator && make clean
	cd median && make clean
	cd phylogeny && make clean
	cd utils && make clean
	cd optkit && make clean
	cd main && make clean

clean_all:
	cd dist && make clean_all
	cd generator && make clean_all
	cd median && make clean_all
	cd phylogeny && make clean_all
	cd utils && make clean_all
	cd optkit && make clean_all
	cd main && make clean_all

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

static __inline__ unsigned long long rdtsc(void);

int 
main(int argc, char **argv){
	printf("%lld ", rdtsc());
	return 0;
}

static __inline__ unsigned long long rdtsc(void)
{
	unsigned hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

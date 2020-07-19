#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
//#include <cimmintrin>

using namespace std;


class time_elem
{
public:
	unsigned long long start;
	unsigned long long end;
	unsigned long long cycle;
	unsigned long long operations;
	unsigned long long bytes;
	double time;
};

class time_list
{
public:
	vector<time_elem> list;	
	int num;
public:
	void init_timer(int id);
	void toc_ser (int id);
	void tic_ser (int id);
	unsigned long long timer_ser ();
#if defined(__i386__)

	static __inline__ unsigned long long rdtsc(void)
	{
		unsigned long long int x;
		__asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
		return x;
	}

#elif defined(__x86_64__)

	static __inline__ unsigned long long rdtsc(void)
	{
		unsigned hi, lo;
		__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
		return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
	}

#endif
};





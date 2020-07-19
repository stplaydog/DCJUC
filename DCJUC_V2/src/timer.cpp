#include "timer.h"



	unsigned long long
time_list::timer_ser ()
{
	return rdtsc();
}

void 
time_list::init_timer(int size){
	num = size;
	int i;
	for(i=0; i<size; i++){
		list[i].start = 0;
		list[i].end = 0;
		list[i].cycle = 0;
		list[i].time = 0;
		list[i].operations = 0;
		list[i].bytes = 0;
	}	
}


	void
time_list::tic_ser (int id)
{
	list[id].start = timer_ser ();
}

	void
time_list::toc_ser (int id)
{
	list[id].end = timer_ser ();
	list[id].cycle += (list[id].end - list[id].start);
	list[id].time += (double)(list[id].end - list[id].start)/(double)2.7e9; 
}


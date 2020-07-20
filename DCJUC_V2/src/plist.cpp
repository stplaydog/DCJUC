#include "plist.h"
#include <iostream>
#include <unistd.h>

template <class T>
inline T my_fetch_add(T *ptr, T val);
template <class T>
inline T my_fetch_set(T *ptr, T val) ;

extern int parallel_mode;

PList::PList(int lb, int ub, int num_t){
	global_lower_bound = lb;	
	global_upper_bound = ub;	
	num_threads = num_t;
	for(int i=0; i<10; i++){
                time_elem e;
                t.list.push_back(e);
                t.init_timer(i);
        }
	thresh = 10000;
}

PList::~PList(){
}

void
PList::bnb_parallel_threads_ub(Instance **ins){
	//this is the sequential part
	p_lists[0]->bnb(ins, 0);
	global_upper_bound = p_lists[0]->upper_bound;
	global_lower_bound = p_lists[0]->lower_bound;
	distribute_works();
	omp_lock_t *writelock = (omp_lock_t*)malloc(sizeof(omp_lock_t)*(num_threads+1));
	for(int i=0; i<=num_threads; i++)
		omp_init_lock(&writelock[i]);
	base = p_lists[0]->base;
	int is_active[num_threads];
	int total_active = num_threads;
	for(int i=0; i< num_threads; i++)
		is_active[i] = true;
	bool giving =false;

	t.tic_ser(0);
	//start the parallel process 
#pragma omp parallel num_threads(num_threads) shared(is_active, total_active)
	{
		int tid = omp_get_thread_num();
		while(total_active>0){
			//check to decrease upper bound
			int buck_id = p_lists[tid]->upper_bound - base;
			if(p_lists[tid]->lower_bound>global_lower_bound){
#pragma omp atomic
				global_lower_bound += (p_lists[tid]->lower_bound-global_lower_bound);
			}
			if(p_lists[tid]->lower_bound<global_lower_bound)
				p_lists[tid]->lower_bound = global_lower_bound;
			if(p_lists[tid]->num[buck_id] == 0 ){ 
				//exhausted
				if(p_lists[tid]->upper_bound<p_lists[tid]->lower_bound){
					if(omp_test_lock(&writelock[tid])){
						if(is_active[tid]==true){
							is_active[tid] = false;
#pragma omp atomic 
							total_active--;
						}
						omp_unset_lock(&writelock[tid]);
					}
					usleep(10);
					buck_id = p_lists[tid]->upper_bound-base;
				}
				else{
					p_lists[tid]->upper_bound--;
				}
			}
			//check others
			if(total_active < num_threads && p_lists[tid]->lower_bound<=p_lists[tid]->upper_bound){
				if(p_lists[tid]->num[buck_id]>thresh){
					for(int i=0; i<num_threads; i++){
						if(i==tid)
							continue;
						if(is_active[i]==false){
							if(omp_test_lock(&writelock[i])){
								for(int j=p_lists[tid]->lower_bound; j<=p_lists[tid]->upper_bound; j++)
									steal_jobs(tid, i, j-base);
								p_lists[i]->upper_bound = p_lists[tid]->upper_bound;
								p_lists[i]->lower_bound = p_lists[tid]->lower_bound;
								//if(p_lists[i]->num[p_lists[i]->upper_bound-base] <=0)
								//	printf("here is wrong!\n");
								struct timespec tim, tim2;
								tim.tv_sec = 0;
								tim.tv_nsec = 500;
								nanosleep(&tim , &tim2);
								is_active[i] = true;
								//printf("upper_bound %d lower_bound %d\n", p_lists[tid]->upper_bound, p_lists[tid]->lower_bound);
								//printf("thresh %d from %d to %d total active %d from has %d to has %d\n", thresh, tid, i, total_active, 
								//		p_lists[tid]->num[p_lists[tid]->upper_bound-base], p_lists[i]->num[p_lists[i]->upper_bound-base]);
#pragma omp atomic 
								total_active ++;
								giving = true;
								omp_unset_lock(&writelock[i]);
							}
							break;
						}
					}
				}
			}
			if(is_active[tid] == false || p_lists[tid]->num[buck_id]<=0)
				continue;
			//this is the real step
			p_lists[tid]->expand_ub(ins, tid, buck_id, writelock[num_threads]);
		}
	}	
	t.toc_ser(0);
	printf("time taken %f\n", t.list[0].time);
	ins[num_threads]->print_encode();
}

void
PList::bnb_parallel_threads_lb(Instance **ins){
	//this is the sequential part
	p_lists[0]->bnb(ins, 0);
	global_upper_bound = p_lists[0]->upper_bound;
	global_lower_bound = p_lists[0]->lower_bound;
	distribute_works(); //it doesn't matter if ub based or lb based
	omp_lock_t *writelock = (omp_lock_t*)malloc(sizeof(omp_lock_t)*(num_threads+1));
	for(int i=0; i<=num_threads; i++)
		omp_init_lock(&writelock[i]);
	base = p_lists[0]->base;
	int is_active[num_threads];
	int total_active = num_threads;
	for(int i=0; i< num_threads; i++)
		is_active[i] = true;
	bool giving =false;

	t.tic_ser(0);
	//start the parallel process 
#pragma omp parallel num_threads(num_threads) shared(is_active, total_active)
	{
		int tid = omp_get_thread_num();
		while(total_active>0){
			//check to decrease upper bound
			int buck_id = p_lists[tid]->lower_bound - base;
			if(p_lists[tid]->upper_bound<global_upper_bound){
#pragma omp atomic
				global_upper_bound += (p_lists[tid]->upper_bound-global_upper_bound);
			}
			if(p_lists[tid]->num[buck_id] == 0 ){ 
				//exhausted
				if(p_lists[tid]->upper_bound<p_lists[tid]->lower_bound){
					if(omp_test_lock(&writelock[tid])){
						if(is_active[tid]==true){
							is_active[tid] = false;
#pragma omp atomic 
							total_active--;
						}
						omp_unset_lock(&writelock[tid]);
					}
					usleep(10);
					buck_id = p_lists[tid]->lower_bound-base;
				}
				else{
					p_lists[tid]->lower_bound++;
				}
			}
			//check others
			if(total_active < num_threads && p_lists[tid]->lower_bound<=p_lists[tid]->upper_bound){
				if(p_lists[tid]->num[buck_id]>thresh){
					for(int i=0; i<num_threads; i++){
						if(i==tid)
							continue;
						if(is_active[i]==false){
							if(omp_test_lock(&writelock[i])){
								//printf("from has %d upper_bound %d lower_bound %d\n", p_lists[tid]->num[p_lists[tid]->lower_bound-base],
								//	p_lists[tid]->upper_bound, p_lists[tid]->lower_bound);
								for(int j=p_lists[tid]->lower_bound; j<=p_lists[tid]->upper_bound; j++){
									steal_jobs(tid, i, j-base);
								}
								p_lists[i]->upper_bound = p_lists[tid]->upper_bound;
								p_lists[i]->lower_bound = p_lists[tid]->lower_bound;
								is_active[i] = true;
								//printf("thresh %d from %d to %d total active %d from has %d to has %d\n", thresh, tid, i, total_active, 
								//		p_lists[tid]->num[p_lists[tid]->lower_bound-base], p_lists[i]->num[p_lists[i]->lower_bound-base]);
#pragma omp atomic 
								total_active ++;
								giving = true;
								omp_unset_lock(&writelock[i]);
							}
							break;
						}
					}
				}
			}
			if(is_active[tid] == false || p_lists[tid]->num[buck_id]<=0)
				continue;
			//this is the real step
			p_lists[tid]->expand_lb(ins, tid, buck_id, writelock[num_threads]);
		}
	}	
	t.toc_ser(0);
	printf("time taken %f\n", t.list[0].time);
	ins[num_threads]->print_encode();
}

void
PList::distribute_works(){
	//copy the main class stuff
	int start = p_lists[0]->lower_bound;
	int end = p_lists[0]->upper_bound;
	for(int i = start; i<=end; i++){
		int buck_id = i-base;
		int avg = p_lists[0]->num[buck_id]/num_threads;
		if(avg==0)
			continue;
		for(int j=1; j<num_threads; j++){
			p_lists[0]->copy(p_lists[j], buck_id, avg);
			p_lists[j]->upper_bound = p_lists[0]->upper_bound;
			p_lists[j]->lower_bound = p_lists[0]->lower_bound;
		}
	}
}

bool
PList::steal_jobs(int from, int to, int buck_id){
	int size = p_lists[from]->num[buck_id]/2;
	if(size==0){
		p_lists[to]->num[buck_id]=0;
		return false;
	}
	p_lists[from]->copy(p_lists[to], buck_id, size);
	return true;
}

template <class T>
inline T my_fetch_add(T *ptr, T val) {
#ifdef GCC_EXTENSION
	return __sync_fetch_and_add(ptr, val);
#endif
#ifdef USE_TPAR
	T t;
#pragma omp atomic capture
	{ t = *ptr; *ptr += val; }
	return t;
#endif
}

template <class T>
inline T my_fetch_set(T *ptr, T val) {
#ifdef GCC_EXTENSION
	return __sync_fetch_and_add(ptr, val);
#endif
	//if(parallel_mode == MODE_TPAR){
	T t;
	//#pragma omp atomic capture
	//	{ t = *ptr; *ptr = val; }
	//	return t;
	//		T t;
	//#pragma omp atomic capture
	//{ t = *ptr; *ptr = val; }
	//return t; 
	//}
	//if(parallel_mode == MODE_SEQ)
	//	return -1;
	//if(parallel_mode == MODE_PAR)
	//	return -1;
}


#pragma once

#ifndef _PLIST_H
#define _PLIST_H


#include "list.h"
#include "instance.h"

using namespace std;

class PList {
public:
	int global_upper_bound;
	int global_lower_bound;
	List **p_lists;
	int num_threads;
	int thresh;
	int base;
	time_list t;
public:
	PList(int lb, int ub, int num_t);
	~PList();
	void bnb_parallel_threads_ub(Instance **ins);	
	void bnb_parallel_threads_lb(Instance **ins);	
	bool steal_jobs(int from, int to, int buck_id);
	void distribute_works();
};

#endif

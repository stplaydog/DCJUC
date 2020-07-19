#pragma once


#ifndef _LIST_H
#define _LIST_H

#include <cstdio>
#include <cstdlib>
#include "omp.h"
#include "instance.h"
#include <vector>
#include <assert.h>
//#include "timer.h"
#define MODE_SEQ 1
#define MODE_PAR 2
#define MODE_TPAR 3
#define LIST_SIZE 1000000
#define CACHE_FILL 64

typedef int myint;
typedef long long lint;

using namespace std;

class List{
public:
	myint upper_bound;
	myint lower_bound;
	myint num_threads;
	myint buck_size; 
    myint base;
	int max_val;
	int min_val;
	bool is_ub;
	bool is_enumerate_all;
	bool is_parallel;
    lint *content_size;
    lint *num;
	int **count;
	int **idx_sum;
	int **end;
	int **start;
	int **l_lb;
	int **l_ub;
	bool **eliminate;
	int **nei_score;
	int *elim_idx;
	time_list t;
	int search_space;
	int parallel_base;
	long long int read_num;
	long long int write_num;
	long long int read_cnt;
	long long int write_cnt;
	bool lk_terminate;
	bool is_lk_max;
	bool is_lk;
    //this for debug purpose only
#ifdef USE_DEBUG
    int **g_count_list;
    int *g_count_idx;
#endif

public:
	//main functions
	List(myint buck_size, myint list_size, 
		myint base, myint num_t, 
		bool is_ub);
	List(myint buck_size, myint list_size, 
		myint base, myint num_t, 
		bool is_ub, bool is_enumerate_all);
	virtual ~List();
	bool compute_partition(int buck_id, Instance** ins);
	void reorder_list(int buck_id, Instance** ins);
	void bnb_parallel_bucket(Instance** ins);
	void bnb(Instance** ins, int ins_id);
	void expand_ub(Instance **ins, int tid, int buck_id, omp_lock_t writelock);
	void expand_lb(Instance **ins, int tid, int buck_id, omp_lock_t writelock);
	void print_buck();
	//lin kernighan algorithm implementation
	void lk(Instance** ins, int heu_level, int term_move, bool is_opt);
	bool add_child_combinations(Instance **ins, int current_level, int use_heu);
	bool add_child_max(Instance **ins, int current_level, int use_heu);
	void reset_list();
    void add_g_count_list(int pos, int g_count);
	//virtual functions
	virtual void add(lint pos, myint buck_id, 
			int num_code, int *encode, int ins_id);
	virtual void get(lint pos, myint buck_id, 
			int *encode, int *num, int ins_id);
	virtual void prepare_parallel_list(int buck_id);
	virtual int copy_zero(int buck_id);
	virtual void reset_num();
	virtual void copy(List *other, int buck_id, int size);
	virtual void print_bucket(int buck_id);
    virtual void free_buck(int buck_id);
    virtual void reallocate_buck(int buck_id);
};

#endif

#pragma once
#ifndef _H_INSKNAP
#define _H_INSKNAP

#include "instance.h"
#include <cstdio>
#include <cstdlib>

class InsKnapsack : public Instance {
public:
	int total_space;
	int value;
	int used_space;
	int num_item;
	int chosen_now;
	int num_branch;
	int min;
	int max;
	int max_now;
	int min_now;
	int *items_space;
	int *items_value;
	bool *in_bag;
	int *branch_id;
	int *sorted;
	int *map;
	InsKnapsack(const char *file, int tid, const char *lf = NULL);
	InsKnapsack(const InsKnapsack &other);
	void to_branch(int which_branch);
	void from_branch();
	void compute_bound();
	int get_num_branches();
	int get_encode(int* encode);
	void to_encode(int* encode, int size);
	void print_encode();
	void copy_code(int *encode, int size);
	int get_value();
private:
	void compute_min_max();
};

#endif

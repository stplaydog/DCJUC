#pragma once
#ifndef _H_INSDCJ
#define _H_INSDCJ

#include "instance.h"
#include <cstdio>
#include <cstdlib>

class InsDCJ : public Instance {
public:
	int c1;
	int c2;
	int size;
	int gene;
	int g1_renamed;
	int rnm_id;
	int old_rnm;
	int last_dup;
	int num_dup[2];
	int num_branch;
	int branch_selected;
	int gf_size;
	int g_size[2];
	int e_size[2];
	//int reverse_idx[1000];
	//order related
	int **order;
	int **renamed_order;
	//duplication related
	int **dup_pos;
	int **dup_idx;
	//branch related
	int *branch_code; //use the first genome as the standard
	int *branch_idx;
	int *branch_pos;
	//some global information
	//for distance computation
	int **adj;
	int *tmp_second;
	int *vet_type;
	//
public:
	InsDCJ(char *file, int tid);
	InsDCJ(const InsDCJ &other);
	InsDCJ(int **order, int one, int two, int size_one, int size_two, int gf_size);
	~InsDCJ();
	void to_branch(int which_branch);
	void from_branch();
	void compute_bound();
	int get_num_branches();
	int get_encode(int* encode);
	void to_encode(int* encode, int size);
	void print_encode();
	void copy_code(int *encode, int size);
	int get_value();
	void from_order_to_adj();
	void vis_graph();
	void print_adj();
	void compute_num_elem();
private:
	int compute_indel_dis();
};

#endif

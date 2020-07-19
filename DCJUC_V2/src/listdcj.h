#pragma once

#include "list.h"

#ifndef _H_LISTDCJ
#define _H_LISTDCJ

using namespace std;

class DCJList : public List{
public:
	int **content;
	int **c_idx;
	int *c_pos;
public:	
	DCJList(myint b_size, myint list_size, 
			myint b, myint num_t, 
			bool is_ub, int num_elem);
	~DCJList();
	void add(lint pos, myint buck_id, 
			int num_code, int *encode, int ins_id);
	void get(lint pos, myint buck_id, 
			int *encode, int *size, int ins_id);
	void prepare_parallel_list(int buck_id);
	int copy_zero(int buck_id);
	void reset_num();
	void copy(List *other, int buck_id, int size);
	void print_bucket(int buck_id);
};

#endif

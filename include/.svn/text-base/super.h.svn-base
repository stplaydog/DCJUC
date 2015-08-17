#pragma once

#include "bitarr.h"
#include "tree.h"
#include "nj.h"


#define TYPE_EDGE_M 0 
#define TYPE_EDGE 1
#define TYPE_DIS 0 
#define TYPE_BOUND 1

#define SEP_USEFUL 0
#define SEP_USELESS 1
#define LEAF_COMMON 2
#define LEAF_INDP 3

struct separation{
	pbarr bit_code;
	pbarr bit_exist;
	int edge_id;
};
typedef struct separation tsep;
typedef struct separation *psep;

struct group{
	int **belong;
	int group_size;
};
typedef struct group tgrp;
typedef struct group *pgrp;

void
run_super_tree(char *file);
void
partition(pdism d_mat, int *p_group);
void init_sub_order(po from, po to, 
		int *p_group, int group_num);

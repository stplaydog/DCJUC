#pragma once

#include "graph.h"

#define LIST_SIZE 1000

struct tree_edge{
	int from;
	int to;
	int num_indel;
	int num_dup;
	int num_dcj;
	int score;
};
typedef struct tree_edge tte;
typedef struct tree_edge *pte;

struct orders{
	int num_genome;
	int num_med_genome;
	int num_input_genome;
	int med_genome_idx;
	int v_size;
	int *e_size;
	int *valid;
	int **v_idx;
	int **e_idx;
};
typedef struct orders to;
typedef struct orders *po;

struct tree{
	int num_edge;
	int edge_idx;
	int score;
	int upper_bound;
	int lower_bound;
	tte *edges;
	po orders;
        int num_genome;
        int edge_num;
};
typedef struct tree tt;
typedef struct tree *pt;


struct tree_list{
	pt **list;
	int *idx_list;
	int *size_list;
	int base;
};
typedef struct tree_list ttl;
typedef struct tree_list *ptl;

pt 
bnb_mp_tree(po o, 
		int upper_bound, int lower_bound);
void 
bnb_mp_sub_tree(po o, 
		int upper_bound, int lower_bound);
void 
free_orders(po o);
void
init_null_tree(po o, pt tree);

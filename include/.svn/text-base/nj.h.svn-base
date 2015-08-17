#pragma once

#include "graph.h"
#include "tree.h"

#define MAX_INT 10000000

struct nj_edge{
	int from;
	int to;
};
typedef struct nj_edge tnje;
typedef struct nj_edge *pnje;

struct distance_mat{
	int num_spc;
	int num_remain_spc;
	double score;
	int num_new_spc;
	double **dis_mat;
	double **q_mat;
	int *map;
	pnje tree;
};
typedef struct distance_mat tdism;
typedef struct distance_mat *pdism;


int 
run_nj(char *file);
int 
run_sub_nj(po o);
int
run_upper_bound(po o);
int 
run_lower_bound(po o);

void run_example();

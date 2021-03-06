#pragma once
#include "tree.h"
#include "graph.h"

#define VALID1 1
#define VALID2 2
#define VALID3 3


struct adj{
	int **list;
	int *idx_list;
	int list_size;
};
typedef struct adj ta;
typedef struct adj *pa;

void 
run_simulate(char *outfile, int gene_len, int gnm_num, 
		double theta, double gamma, double phi);
void
init_adj(pa a, int size1, int size2);

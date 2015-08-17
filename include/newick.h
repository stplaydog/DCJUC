#pragma once
#include "tree.h"

#define E_TYPE_INIT 1
#define E_TYPE_UNINIT 0
#define V_TYPE_INIT 1
#define V_TYPE_UNINIT 0
#define MAX_NW_LEN 10000

#define ERROR_PRINT() {printf("Error on (%d) line in (%s) file\n", __LINE__, __FILE__); exit(123);}

struct newtree
{
	int num_v;
	int *v_type;
	int **edges;
	int **e_type;
	int *e_idx;
	char **v_string;
};
typedef struct newtree tnt;
typedef struct newtree *pnt;

void 
from_edge_to_newick(FILE *writer, pt trees, 
                        int *map);
void
from_newick_to_edge(char *file, pt tree);

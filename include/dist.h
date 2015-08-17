#include "graph.h"
#include "tree.h"

#define MAX_VET_SIZE 32768
#define PI_OPEN 1
#define GAMMA_OPEN 2
#define CYCLE 3
#define END 4
#define NONE 5
#define NUL -1


struct bnb_elem{
	int edge_one;
	int edge_two;
	int edge_three;
	int edge_four;
};
typedef struct bnb_elem t_elem;
typedef struct bnb_elem *p_elem;

struct exemplar_csr{
	int v_size;
	int e_size;
	int *v_idx;
	int *e_idx;
	int *adj_mat;
	int *e_check;
};
typedef struct exemplar_csr t_exem;
typedef struct exemplar_csr *p_exem;

struct search_pos{
	int pos;
	p_elem elems;
	int size;
	int idx;
	struct search_pos *next;
	struct search_pos *prev;
};
typedef struct search_pos t_pos;
typedef struct search_pos *p_pos;

void compute_indel_exemplar_dis(pg g, int c1, int c2, int *dis, int gene, int pos);
int compute_indel_dis(int **adj, int *vet_type, int c1, int c2, int size, int gene);
void to_adj_arr(pg g, int **adj, int *vet_type, int c1, int c2);
int evaluate(pg g);
int comupte_all_possible_exemplar(pg g, int c1, 
		int c2, int gene);
int 
compute_all_possible_exemplar_orders(po o, int c1, 
		int c2, int gene);
void init_for_dis();
int 
cpec_mem_save(pg g1, int c1, int c2, int gene);

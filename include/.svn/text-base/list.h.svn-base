#pragma once

#include "graph.h"

#define CONN_TYPE_ONE 0
#define CONN_TYPE_TWO 1
#define AVG_CHILD_PER_PARENT 10

struct list_elem{
	int from_one;
	int to_one;
	int from_two;
	int to_two;
	int new_connect_type; 
	//type 1 from_one -> to_two
	//type 2 from_one -> from_two 
	int idx_parent;
};
typedef struct list_elem tle;
typedef struct list_elem *ple;

struct list_bnb_elem{
	int from;
	int to;
	int idx_parent;
};
typedef struct list_bnb_elem tbe;
typedef struct list_bnb_elem *pbe;

struct search_list{
	int buck_size;
	int *list_size;
	int p_size;
	int *idx_child;
	int *child_size;
	tle **child;
	int ***parent;
	int **parent_size;
	int *idx_parent;
	int *f_check;
};
typedef struct search_list tsl;
typedef struct search_list *psl;

struct bnb_list{
	int buck_size;
	int list_size;
	int p_size;
	int *idx_parent;
	int *idx_child;
	int ***parent;
	int **parent_size;
	int *child_size;
	tbe **child;
};
typedef struct bnb_list tbl;
typedef struct bnb_list *pbl;

void
create_search_list(int buck_size, int list_size, 
		int p_size, psl l);
int
add_parent(pg g, psl l,
		int buck_id);
void
add_child(int from_one, int to_one, 
		int from_two, int to_two, 
		int conn_type, int idx_parent,
		int buck_id, psl l);
int retrive(pg g, psl l, int current_level);
void
double_list(psl l, int buck_id);
void reset_list(psl l);

//search list functions for bnb algorithm
void
create_search_list_bnb(int buck_size, int list_size, 
		int p_size, pbl l);
void
add_nodes_bnb(pg g, pbl l, int base,
		int *upper_bound, int* lower_bound);
int retrive_bnb(pg g, pbl l, int current_level);
void
double_list_bnb(pbl l, int buck_id);
int 
get_list_size(psl l);

void free_list(psl l);

#include "bnb.h"


void get_bounds(pg g, int* upper_bound, int* lower_bound);

void bnb(pg g){
	int upper_bound;
	int lower_bound;
	get_bounds(g, &upper_bound, &lower_bound);
	//create bnb search list
	pbl l = (pbl)malloc(sizeof(tbl));	
	int buck_size = upper_bound - lower_bound + 1;
	int list_size = 1000; 
	int p_size = g->e_size[3];
	create_search_list_bnb(buck_size, list_size, p_size, l);
	int base = lower_bound;
	//put list
	add_nodes_bnb(g, l, base, &upper_bound, &lower_bound);
	while(upper_bound > lower_bound){
		retrive_bnb(g, l, upper_bound);
		add_nodes_bnb(g, l, base, &upper_bound, &lower_bound);
	}
	printf("the result is %d", upper_bound);
}

void get_bounds(pg g, int* upper_bound, int* lower_bound){
	int sum =0; 
	int dis1, dis2, dis3;
	int effective_size=0, i, j, color;
	for(i=0;i<g->v_size;i++)
		if(g->v_check[i] != NUL)
			effective_size++;
	prepare_dis(g);
	dis1 = compute_all_possible_exemplar(g, 0, 3, effective_size);
	dis2 = compute_all_possible_exemplar(g, 1, 3, effective_size);
	dis3 = compute_all_possible_exemplar(g, 2, 3, effective_size);
	sum = dis1+dis2+dis3;	
	int min=0;
	if(dis1<dis2 && dis1<dis3)
		min=dis1;
	else if(dis2 < dis3)
		min = dis2;
	else
		min = dis3;
	*upper_bound = sum/2;
	*lower_bound = sum - min;
}

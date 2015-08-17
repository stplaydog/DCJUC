#include "nj.h"
#include "tree.h"

int debug_nj = FALSE;

void
init_dis_mat(pdism d_mat, int size);
void
compute_distance_matrix(po o, pdism d_mat);
void
compute_q_matrix(pdism d_mat);
void 
join_two_car(pdism d_mat, po o);
void 
join_two_simp(pdism d_mat);
void 
re_order_matrix_car(pdism d_mat, po o,
		int from, int to);
void 
re_order_matrix_simp(pdism d_mat, int from, int to);
void
output_tree(pdism d_mat);
void 
init_dis_mat_by_interm_orders(pdism d_mat, po o);
int
compute_min_disjoint_distance(po o);
int 
compute_order_score(po o, int *map, 
		int g_size);
int 
perform_switch(po o, int *map, 
		int pos, int g_size);
void
print_dis_mat(pdism d_mat);
void
print_q_mat(pdism d_mat);
void init_example(pdism d_mat);
void
free_d_mat(pdism d_mat);
void 
init_dup_car(po o, int c1, int c2);

int edge_count;

/********************************************************
 * this is the main entry for the neighbor joining method
 * ******************************************************/
int 
run_nj(char *file){
	int i, j;
	po o = (po)malloc(sizeof(to));	
	init_orders(file, o);
	pdism d_mat = (pdism)malloc(sizeof(tdism));
	init_dis_mat(d_mat, o->num_input_genome);
	compute_distance_matrix(o, d_mat); // we only need to compute distance once
	compute_q_matrix(d_mat);
	edge_count = 0;
	while(d_mat->num_remain_spc > 2){
#ifdef SIMP_NJ
		join_two_simp(d_mat);
#else
		join_two_car(d_mat, o);
#endif
		compute_q_matrix(d_mat);
	}
	//join the last two
	d_mat->score += d_mat->dis_mat[0][1];
	pnje tmp = d_mat->tree;
	int result = d_mat->score;
	free_d_mat(d_mat);
	free_orders(o);
	return result;
}

/****************************************************
 *this is to compute the upper bound of the tree
 * after performing some maximum parsimony operations
 * *************************************************/
int 
run_upper_bound(po o){
	int i, j;
	po oo = (po)malloc(sizeof(to));
	int num_remain = 0;
	for(i=0; i<o->num_genome; i++){
		if(o->valid[i] == FALSE)
			continue;
		num_remain++;
	}
	if(num_remain==0)
		return;
	oo->num_genome = o->num_genome;
	oo->num_med_genome = o->num_med_genome;
	oo->num_input_genome = o->num_input_genome;
	oo->med_genome_idx = o->med_genome_idx;
	oo->v_size = o->v_size;
	oo->v_idx = (int**)malloc(sizeof(int*)*o->num_genome);
	oo->e_idx = (int**)malloc(sizeof(int*)*o->num_genome);
	oo->e_size = (int*)malloc(sizeof(int)*o->num_genome);
	oo->valid = (int*)malloc(sizeof(int)*o->num_genome);
	for(i=0; i<o->num_genome; i++){
		oo->valid[i] = o->valid[i];
		oo->e_size[i] = o->e_size[i];
		oo->v_idx[i] = (int*)malloc(sizeof(int)*o->v_size);
		if(o->e_size[i] <= 0)
			continue;
		oo->e_idx[i] = (int*)malloc(sizeof(int)*o->e_size[i]);
		for(j=0;j<o->v_size; j++)
			oo->v_idx[i][j] = o->v_idx[i][j];
		for(j=0; j<o->e_size[i]; j++)
			oo->e_idx[i][j] = o->e_idx[i][j];
	}
	//
	pdism d_mat = (pdism)malloc(sizeof(tdism));
	init_dis_mat_by_interm_orders(d_mat, oo);
	compute_q_matrix(d_mat);
	while(d_mat->num_remain_spc > 2){
#ifdef SIMP_NJ
		join_two_simp(d_mat);
#else
		join_two_car(d_mat, oo);
		//printf("%d %d %p %p\n", oo->v_idx[i][0], oo->e_size[i], oo->v_idx[i], oo->e_idx[i]);
#endif
		compute_q_matrix(d_mat);
	}
	//join the last two
	d_mat->score += d_mat->dis_mat[0][1];
	pnje tmp = d_mat->tree;
	free_orders(oo);
	int result = d_mat->score;
	free_d_mat(d_mat);
	return result;
}

/* ***********************************************
 * a local search algorithm to get the lower bound 
 * ***********************************************/
int 
run_lower_bound(po o){
	int i, j;
	int best = MAX_INT;
	//init an order
	//perform local search strategy
	int *map = (int*)malloc(sizeof(int)*o->num_genome);
	for(i=0; i<o->num_genome; i++)
		map[i] = i;
	int improved = TRUE;
	int switch_num =0;
	//while(improved == TRUE && switch_num++<10){
	//	for(i=0;i<o->num_genome; i++){
	//		//printf("%d ", i);
	//		if(o->valid[map[i]] == FALSE)
	//			continue;
	//		improved = perform_switch(o, map, i, o->v_size/2);
	//		if(improved == TRUE)
	//			break;
	//	}
	//	//printf("\n");
	//}
	//printf("come here\n");
	//for(i=0; i<o->num_genome; i++)
	//	printf("%d ", map[i]);
	//printf("\n");
	//compute the final score
	int score;	
	score = compute_order_score(o, map, o->v_size/2);
	free(map);
	return score;
}

/*
 *
 * **/
int 
perform_switch(po o, int *map, 
		int pos, int g_size){
	int i, j;
	int order = map[pos];
	int order_pre = map[pos==0?o->num_genome-1:pos-1];
	int pre = pos==0?o->num_genome-1:pos-1;
	i=2;
	while(o->valid[order_pre] == FALSE && order_pre!= order){
		order_pre = map[pos-i<0?o->num_genome+(pos-i):pos-i];
		pre = pos-i<0?o->num_genome+(pos-i):pos-i;
		i++;
	}
	int order_neigh = map[(pos+1)%(o->num_genome)];
	int neigh = (pos+1)%(o->num_genome);
	i=2;
	while(o->valid[order_neigh] == FALSE && order_neigh != order){
		order_neigh = map[(pos+i)%(o->num_genome)];
		neigh = (pos+i)%(o->num_genome);
		i++;
	}
	int order_neigh_after = map[(pos+2)%(o->num_genome)];
	int neigh_after = (pos+2)%(o->num_genome);
	i=3;
	while(o->valid[order_neigh_after] == FALSE && order_neigh_after != order_neigh){
		order_neigh_after = map[(pos+i)%(o->num_genome)];
		neigh_after = (pos+i)%(o->num_genome);
		i++;
	}
	//for(i=0; i<o->num_genome; i++)
	//	printf("%d ", o->valid[i]);
	//printf("\n");
	//printf("order %d pre %d neigh %d neigh_after%d \n", order, order_pre, order_neigh, order_neigh_after);
	int o_score = 
		compute_all_possible_exemplar_orders(o, order, order_pre, g_size)+
		compute_all_possible_exemplar_orders(o, order_neigh, order_neigh_after, g_size);
	int n_score = 
		compute_all_possible_exemplar_orders(o, order_pre, order_neigh, g_size)+
		compute_all_possible_exemplar_orders(o, order, order_neigh_after, g_size);
		//printf("come here\n");
	//perform the switch
	if(n_score > o_score){
		map[pos] = order_neigh;
		map[neigh] = order;
		return TRUE;
	}
	return FALSE;
}

/*
 *
 * **/
int 
compute_order_score(po o, int *map, 
			int g_size){
	int i, j;
	int score =0;
	for(i=0;i<o->num_genome;i++){
		if(o->valid[map[i]] == FALSE)
			continue;
		int pos = i;
		int order = map[pos];
		int order_neigh = map[(pos+1)%(o->num_genome)];
		j=2;
		while(o->valid[order_neigh] == FALSE){
			order_neigh = map[(pos+j)%(o->num_genome)];
			j++;
		}
		//printf("order %d order_neigh %d\n", order, order_neigh);
		score += compute_all_possible_exemplar_orders(o, order, order_neigh, g_size);
	}
	return score/2;
}

/*****************************
 * init data structures
 * ***************************/
void
init_dis_mat(pdism d_mat, int size){
	int i;
	d_mat->num_spc = size;
	d_mat->num_remain_spc = size;
	d_mat->score = 0;
	d_mat->num_new_spc = 0;

	d_mat->dis_mat = (double**)malloc(sizeof(double*)*size);
	d_mat->q_mat = (double**)malloc(sizeof(double*)*size);
	for(i=0; i<size; i++){
		d_mat->dis_mat[i] = (double*)malloc(sizeof(double)*size);
		d_mat->q_mat[i] = (double*)malloc(sizeof(double)*size);
	}

	d_mat->tree = (pnje)malloc(sizeof(tnje)*d_mat->num_spc*2); 

	d_mat->map = (int*)malloc(sizeof(int)*size);
	for(i=0; i<size; i++)
		d_mat->map[i] = i;
}

/************************************************************
 * this is to initilizing the distance matrix for input 
 * genomes that some of the genomes have already been joined 
 * to the mp tree
 * *********************************************************/
void 
init_dis_mat_by_interm_orders(pdism d_mat, po o){
	int i, j;	
	int num_remain = 0;
	for(i=0; i<o->num_genome; i++)
		if(o->valid[i] == TRUE)
			num_remain++;
	int *rev_map = (int*)malloc(sizeof(int)*o->num_genome);
	int idx=0;
	init_dis_mat(d_mat, num_remain);	
	for(i=0; i<o->num_genome; i++)
		if(o->valid[i]==TRUE){
			d_mat->map[idx] = i;
			rev_map[i] = idx;
			idx++;
		}
	//for(i=0; i<num_remain; i++)
	//	printf("%d:%d ", i, d_mat->map[i]);
	//printf("\n");
	for(i=0; i<o->num_genome; i++){
		for(j=i+1;j<o->num_genome; j++){
			if(o->valid[i]==TRUE && o->valid[j]==TRUE)
				d_mat->dis_mat[rev_map[j]][rev_map[i]] = d_mat->dis_mat[rev_map[i]][rev_map[j]] = compute_all_possible_exemplar_orders(o, i, j, o->v_size/2);
		}
		if(o->valid[i]==TRUE)
			d_mat->dis_mat[rev_map[i]][rev_map[i]] = 0;
	}
	free(rev_map);
}

/****************************************************************
 *compute the min distance between the set of used order and the set
 * of unused order.
 * ****************************************************************/
int
compute_min_disjoint_distance(po o){
	int i, j;	
	int min = MAX_INT;
	for(i=0;i<o->num_input_genome; i++){
		if(o->valid[i] == TRUE)
			continue;
		for(j=i+1; j<o->num_input_genome; j++){
			if(o->valid[j] == FALSE)
				continue;
			int dis = compute_all_possible_exemplar_orders(o, i, j, o->v_size/2);
			if(dis < min)
				min = dis;
		}
	}
	return min;
}

/* **********************************************************************
 * to begin with, we need to compute a distance matrix for the future use
 * **********************************************************************/
void
compute_distance_matrix(po o, pdism d_mat){
	int i, j;	
	int size = o->num_input_genome;
	for(i=0; i<size; i++){
		for(j=i+1;j<size; j++){
			d_mat->dis_mat[j][i] = d_mat->dis_mat[i][j] = compute_all_possible_exemplar_orders(o, i, j, o->v_size/2);
		}
		d_mat->dis_mat[i][i] = 0;
	}
	//for(i=0; i<size; i++){
	//	for(j=0;j<size; j++){
	//		printf("%d ", d_mat->dis_mat[j][i]);
	//	}
	//	printf("\n");
	//}
}

/* ****************************************
 * q matrix to find which two edges to join
 * ****************************************/
void
compute_q_matrix(pdism d_mat){
	int i, j, k;	
	int size = d_mat->num_remain_spc;
	for(i=0;i<size; i++){
		for(j=i+1;j<size; j++){
			double d_ij = d_mat->dis_mat[i][j];	
			double d_ik=0;
			double d_jk=0;
			for(k=0; k<size; k++){
				d_ik+=d_mat->dis_mat[i][k];
				d_jk+=d_mat->dis_mat[j][k];
				//printf("%3.1f|%3.1f ", d_mat->dis_mat[i][k], d_mat->dis_mat[j][k]);
			}
			d_mat->q_mat[i][j] = d_mat->q_mat[j][i] = (size-2)*d_ij-d_ik-d_jk;
			
			//printf("d_ij %3.1f d_ik %3.1f d_jk %3.1f qmat %3.1f\n", (size-2)*d_ij, d_ik, d_jk, d_mat->q_mat[i][j]);
		}
		d_mat->q_mat[i][i] = 0;
	}
}

/* ***************************************************************
 * this is the key part of the neighbor join algorithm
 * 1) find two species
 * 2) join these two
 * 3) update the matrix including a) reorder b) add one column/row
 * ***************************************************************/
void 
join_two_car(pdism d_mat, po o){
	int i, j, k;
	int min_i, min_j;
	double min_q = MAX_INT;
	int size = d_mat->num_remain_spc;
	//find min
	for(i=0;i<size; i++){
		for(j=i+1;j<size; j++){
			if(d_mat->q_mat[i][j]<min_q){
				min_q = d_mat->q_mat[i][j];
				min_i = i;
				min_j = j;
			}
		}
	}
	int c1 = d_mat->map[min_i];
	int c2 = d_mat->map[min_j];
	//compute the distance of new_v to min_i and min_j
	//double dis_ij = d_mat->dis_mat[min_i][min_j];
	//double dis_ik = 0;
	//double dis_jk = 0;
	//for(k=0; k<size; k++){
	//	dis_ik += d_mat->dis_mat[min_i][k];
	//	dis_jk += d_mat->dis_mat[min_j][k];
	//}	
	//double dis_new_i = dis_ij / 2 + (dis_ik - dis_jk) / (2 * (size - 2));
	//double dis_new_j = dis_ij - dis_new_i;
	//d_mat->score += dis_ij;	
	//printf("c1 %d c2 %d %p\n", c1, c2, o->v_idx[0]);
	init_dup_car(o, c1, c2);
	//if(c1==3 && c2==4)
	//{
	//	FILE *writer = fopen("nj.dot", "w");
	//	int pos = o->med_genome_idx-1;
	//	for(i=0; i<o->v_size; i++){
	//		int start = i==0?0:o->v_idx[pos][i-1];
	//		int end = o->v_idx[pos][i];
	//		int from = i;
	//		for(j=start; j<end; j++){
	//			int to = o->e_idx[pos][j];
	//			fprintf(writer, "%d -- %d [color=red];\n", from, to);
	//		}
	//		start = i==0?0:o->v_idx[c1][i-1];
	//		end = o->v_idx[c1][i];
	//		for(j=start; j<end; j++){
	//			int to = o->e_idx[c1][j];
	//			fprintf(writer, "%d -- %d [color=red];\n", from, to);
	//		}
	//	}
	//	fclose(writer);
	//}	
	//for(i=0; i<d_mat->num_remain_spc; i++)
	//	printf("%d ", d_mat->map[i]);
	//printf("\n");
	double dis_new_i = compute_all_possible_exemplar_orders(o, o->med_genome_idx-1, c1, o->v_size/2);
	double dis_new_j = compute_all_possible_exemplar_orders(o, o->med_genome_idx-1, c2, o->v_size/2);
	edge_count ++;
	d_mat->score += (dis_new_i + dis_new_j);
	//reorder the matrix
	re_order_matrix_car(d_mat, o, min_i, min_j);
	//printf("finished computing c2\n");
	//join the new vertex with min_i and min_j
	int new_node_id = d_mat->num_spc + d_mat->num_new_spc++;
	pnje tmp = d_mat->tree;

}


void 
join_two_simp(pdism d_mat){
	int i, j, k;
	int min_i, min_j;
	double min_q = MAX_INT;
	int size = d_mat->num_remain_spc;
	//find min
	for(i=0;i<size; i++){
		for(j=i+1;j<size; j++){
			if(d_mat->q_mat[i][j]<min_q){
				min_q = d_mat->q_mat[i][j];
				min_i = i;
				min_j = j;
			}
		}
	}
	//compute the distance of new_v to min_i and min_j
	double dis_ij = d_mat->dis_mat[min_i][min_j];
	double dis_ik = 0;
	double dis_jk = 0;
	for(k=0; k<size; k++){
		dis_ik += d_mat->dis_mat[min_i][k];
		dis_jk += d_mat->dis_mat[min_j][k];
	}	
	double dis_new_i = dis_ij / 2 + (dis_ik - dis_jk) / (2 * (size - 2));
	double dis_new_j = dis_ij - dis_new_i;
	d_mat->score += dis_ij;	
	//init_dup_car(o, c1, c2);
	//double dis_new_i = compute_all_possible_exemplar_orders(o, o->med_genome_idx-1, c1, o->v_size/2);
	//double dis_new_j = compute_all_possible_exemplar_orders(o, o->med_genome_idx-1, c2, o->v_size/2);
	//edge_count ++;
	//printf("after join %dth edge, the score is %3.1f\n", edge_count, d_mat->score);
	//d_mat->score += (dis_new_i + dis_new_j);
	//reorder the matrix
	re_order_matrix_simp(d_mat, min_i, min_j);
	//join the new vertex with min_i and min_j
	int new_node_id = d_mat->num_spc + d_mat->num_new_spc++;
	pnje tmp = d_mat->tree;


}

/* ******************************************************
 * after find two species we need to delete them from the 
 * data structure, and reorder the dis_mat and q_mat
 * ******************************************************/
void 
re_order_matrix_car(pdism d_mat, po o,
		int from, int to){
	int i, j, k;
	int new_i=0, new_j=0;
	double **tmp_mat = (double**)malloc(sizeof(double*)*(d_mat->num_remain_spc-1)); 
	for(i=0; i<d_mat->num_remain_spc; i++){
		if(i == from || i == to)
			continue;
		tmp_mat[new_i] = (double*)malloc(sizeof(double)*(d_mat->num_remain_spc-1));
		d_mat->map[new_i] = d_mat->map[i];
		new_j = 0;
		for(j=0; j<d_mat->num_remain_spc; j++){
			if(j == from || j == to)
				continue;
			tmp_mat[new_i][new_j] = d_mat->dis_mat[i][j];
			new_j++;
		}
		new_i++;
	}
	tmp_mat[new_i] = (double*)malloc(sizeof(double)*(d_mat->num_remain_spc-1));
	//add one new row and new column
	int new_pos = d_mat->num_remain_spc-2;
	int d_from_to = d_mat->dis_mat[from][to];
	new_i = 0;
	for(k=0; k<d_mat->num_remain_spc; k++){
		if(k==from || k==to)
			continue;
	//	double d_from_k = d_mat->dis_mat[from][k];
	//	double d_to_k = d_mat->dis_mat[to][k];
		//printf("---c1 %d ---c2 %d\n", o->med_genome_idx-1, d_mat->map[k]);
		double dis_new_k = compute_all_possible_exemplar_orders(o, o->med_genome_idx-1, d_mat->map[k], o->v_size/2);
		//tmp_mat[new_i][new_pos] = tmp_mat[new_pos][new_i] = (d_from_k + d_to_k - d_from_to) / 2;
		tmp_mat[new_i][new_pos] = tmp_mat[new_pos][new_i] = dis_new_k;
		new_i++;
	}
	tmp_mat[new_pos][new_pos]=0;
	d_mat->map[new_pos] = d_mat->num_spc + d_mat->num_new_spc;
	//copy back
	for(i=0; i<=new_pos; i++)
		for(j=0; j<=new_pos; j++)
			d_mat->dis_mat[i][j] = tmp_mat[i][j];
	d_mat->num_remain_spc-=1;
	for(i=0; i<d_mat->num_remain_spc; i++)
		free(tmp_mat[i]);
	free(tmp_mat);
}

void 
re_order_matrix_simp(pdism d_mat, int from, int to){
	int i, j, k;
	int new_i=0, new_j=0;
	double **tmp_mat = (double**)malloc(sizeof(double*)*(d_mat->num_remain_spc-1)); 
	for(i=0; i<d_mat->num_remain_spc; i++){
		if(i == from || i == to)
			continue;
		tmp_mat[new_i] = (double*)malloc(sizeof(double)*(d_mat->num_remain_spc-1));
		d_mat->map[new_i] = d_mat->map[i];
		new_j = 0;
		for(j=0; j<d_mat->num_remain_spc; j++){
			if(j == from || j == to)
				continue;
			tmp_mat[new_i][new_j] = d_mat->dis_mat[i][j];
			new_j++;
		}
		new_i++;
	}
	tmp_mat[new_i] = (double*)malloc(sizeof(double)*(d_mat->num_remain_spc-1));
	//add one new row and new column
	int new_pos = d_mat->num_remain_spc-2;
	int d_from_to = d_mat->dis_mat[from][to];
	new_i = 0;
	for(k=0; k<d_mat->num_remain_spc; k++){
		if(k==from || k==to)
			continue;
		double d_from_k = d_mat->dis_mat[from][k];
		double d_to_k = d_mat->dis_mat[to][k];
		//double dis_new_k = compute_all_possible_exemplar_orders(o, o->med_genome_idx-1, d_mat->map[k], o->v_size/2);
		tmp_mat[new_i][new_pos] = tmp_mat[new_pos][new_i] = (d_from_k + d_to_k - d_from_to) / 2;
		//tmp_mat[new_i][new_pos] = tmp_mat[new_pos][new_i] = dis_new_k;
		new_i++;
	}
	tmp_mat[new_pos][new_pos]=0;
	d_mat->map[new_pos] = d_mat->num_spc + d_mat->num_new_spc;
	//copy back
	for(i=0; i<=new_pos; i++)
		for(j=0; j<=new_pos; j++)
			d_mat->dis_mat[i][j] = tmp_mat[i][j];
	d_mat->num_remain_spc-=1;
	for(i=0; i<d_mat->num_remain_spc-1; i++)
		free(tmp_mat[i]);
	free(tmp_mat);
}

/* ***********************************************
 * print the phylogenetic tree of the given(nexus) format
 * ***********************************************/
void
output_tree(pdism d_mat){
}

void
print_dis_mat(pdism d_mat){
	int i, j, k;
	int size = d_mat->num_remain_spc;
	for(i=0; i<size; i++){
		for(j=0; j<size; j++){
			printf("%3.1f\t", d_mat->dis_mat[i][j]);
		}
		printf("\n");
	}
	printf("\n\n\n");
}

void
print_q_mat(pdism d_mat){
	int i, j, k;
	int size = d_mat->num_remain_spc;
	for(i=0; i<size; i++){
		for(j=0; j<size; j++){
			printf("%3.1f\t", d_mat->q_mat[i][j]);
		}
		printf("\n");
	}
	printf("\n\n\n");
}


void run_example(){
	int i;
	pdism d_mat = (pdism)malloc(sizeof(tdism));
	d_mat->num_spc = 4;
	d_mat->num_remain_spc = 4;
	d_mat->score = 0;
	d_mat->num_new_spc = 0;
	d_mat->dis_mat = (double**)malloc(sizeof(double*)*4);
	d_mat->q_mat = (double**)malloc(sizeof(double*)*4);
	d_mat->map = (int*)malloc(sizeof(int)*4);
	for(i=0; i<4; i++){
		d_mat->dis_mat[i] = (double*)malloc(sizeof(double)*4);
		d_mat->q_mat[i] = (double*)malloc(sizeof(double)*4);
	}
	d_mat->dis_mat[0][0] = 0;
	d_mat->dis_mat[0][1] = 7;
	d_mat->dis_mat[0][2] = 11;
	d_mat->dis_mat[0][3] = 14;
	d_mat->dis_mat[1][0] = 7;
	d_mat->dis_mat[1][1] = 0;
	d_mat->dis_mat[1][2] = 6;
	d_mat->dis_mat[1][3] = 9;
	d_mat->dis_mat[2][0] = 11;
	d_mat->dis_mat[2][1] = 6;
	d_mat->dis_mat[2][2] = 0;
	d_mat->dis_mat[2][3] = 7;
	d_mat->dis_mat[3][0] = 14;
	d_mat->dis_mat[3][1] = 9;
	d_mat->dis_mat[3][2] = 7;
	d_mat->dis_mat[3][3] = 0;

	d_mat->tree = (pnje)malloc(sizeof(tnje)); 
	d_mat->tree->from = -1;
	d_mat->tree->to = -1;

	for(i=0; i<4; i++)
		d_mat->map[i] = i;

	compute_q_matrix(d_mat);
	print_dis_mat(d_mat);
	print_q_mat(d_mat);
	edge_count = 0;
	while(d_mat->num_remain_spc > 2){
#ifdef SIMP_NJ
		join_two_simp(d_mat);
#else
		//join_two_car(d_mat, o);
#endif
		compute_q_matrix(d_mat);
		print_dis_mat(d_mat);
		print_q_mat(d_mat);
	}
}

void
free_d_mat(pdism d_mat){
	int i;
	for(i=0; i<d_mat->num_spc; i++){
		free(d_mat->dis_mat[i]);
		free(d_mat->q_mat[i]);
	}
	free(d_mat->dis_mat);
	free(d_mat->q_mat);
	free(d_mat->map);
	free(d_mat->tree);
	free(d_mat);
}

void 
init_dup_car(po o, int c1, int c2){
	int i, j, k;
	int *med_deg = (int*)malloc(sizeof(int)*o->v_size);
	int sum =0;
	int pos = o->med_genome_idx++;
	o->valid[pos] = TRUE;
	for(i=0; i<o->v_size; i++){
		int start_one = i==0?0:o->v_idx[c1][i-1];
		int end_one = o->v_idx[c1][i];
		int deg_one = end_one -start_one;
		int start_two = i==0?0:o->v_idx[c2][i-1];
		int end_two = o->v_idx[c2][i];
		int deg_two = end_two -start_two;
		med_deg[i] = deg_one>deg_two?deg_two:deg_one;
		sum += med_deg[i];
		o->v_idx[pos][i] = sum;
	}
	if(sum%2 == 1){
		int last = o->v_size-1;
		for(i=0; i<o->v_size; i++)
			if(med_deg[i] != 0)
				last = i;
		o->v_idx[pos][last] = o->v_idx[pos][last] - 1;
	}
	o->e_idx[pos] = (int*)malloc(sizeof(int)*sum);
	o->e_size[pos] = sum;
	for(i=0;i<o->e_size[pos]; i++)
		o->e_idx[pos][i] = -1;
	//init car
	for(i=0; i<o->v_size; i++){
		int start_one = i==0?0:o->v_idx[c1][i-1];
		int end_one = o->v_idx[c1][i];
		int deg_one = end_one -start_one;
		int start_two = i==0?0:o->v_idx[c2][i-1];
		int end_two = o->v_idx[c2][i];
		int deg_two = end_two -start_two;
		int start_med = i==0?0:o->v_idx[pos][i-1];
		if(deg_one == 1 && deg_two == 1){
			if(o->e_idx[c1][start_one] == o->e_idx[c2][start_two]
				&& o->e_idx[c1][start_one] != CAP){
				o->e_idx[pos][start_med] = o->e_idx[c2][start_two];
				int start_med_to = o->e_idx[c2][start_two]==0?0:o->v_idx[pos][o->e_idx[c2][start_two]-1];
				int end_med_to = o->v_idx[pos][o->e_idx[c2][start_two]];
				for(j= start_med_to; j<end_med_to; j++){
					if(o->e_idx[pos][j] == -1){
						o->e_idx[pos][j] = i;
						break;
					}
				}
			}
		}
	}
	//for others randomly choose
	for(i=0; i<o->v_size; i++){
		int start = i==0?0:o->v_idx[pos][i-1];
		int end = o->v_idx[pos][i];
		for(j=start; j<end; j++){
			if(o->e_idx[pos][j] == -1){
				int to;
				int to_e;
				int count =0;
				while(1){
					count++;
					if(count==100000){
						//int m,n;
						//FILE *writer = fopen("conf.dot", "w");
						//fprintf(writer, "Graph{\n");
						//for(m=0; m<o->v_size; m++){
						//	int start = m==0?0:o->v_idx[pos][m-1];
						//	int end = o->v_idx[pos][m];
						//	for(n=start; n<end; n++){
						//		int from = m;
						//		int to = o->e_idx[pos][n];
						//		fprintf(writer, "%d -- %d [color=red];\n", from, to);
						//	}
						//}
						//fprintf(writer, "}\n");
						//fclose(writer);
						//exit(1);
						to = CAP;
						break;
					}
					to = rand() % o->v_size;
					if(to != i){
						int start = to==0?0:o->v_idx[pos][to-1];
						int end = o->v_idx[pos][to];
						int k;
						int found = FALSE;
						for(k=start;k<end;k++)
							if(o->e_idx[pos][k] == -1){
								found = TRUE;
								to_e = k;
							}
						if(found == TRUE)
							break;
					}
				}
				o->e_idx[pos][j] = to;
				o->e_idx[pos][to_e] = i;
			}
		}
	}
	//printf("come here\n");
	free(med_deg);
}

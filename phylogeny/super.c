#include "super.h"
#include "spectral.h"

int 
check_over_lap(pt tree_to, pt tree_from, 
		int *map_to, int *map_from, 
		int size_from, int size_to);
int
merge_tree_separate(pt tree_to, pt tree_from, 
		int *map_to, int *map_from, 
		int size_from, int size_to, 
		pdism d_mat);
int
init_orders_by_group(po to, po from, 
		int group_num, int *map, 
		pgrp group_id);
void 
init_order_by_edge(pg g, pt tree, 
int e_id, int type);
void 
compute_laplacian_mat(pdism d_mat);
int
merge_tree_overlap(tsep *sep_from, tsep *sep_to,
		pt tree_to, pt tree_from, 
		int *map_to, int *map_from, 
		int size_from, int size_to);
void
add_independent_edge(tte *edges, int edge_idx, 
		int from, int to);
void
separate(tte *edges, int edge_idx, 
		tsep *sep, int sep_idx);
void
merge_sep_by_edges(pt tree_from, pt tree_to, 
		int best_e_from, int best_e_to);
void
calculate_separation_fp(psep sep, pt t, 
		int *map, int e_id);
void 
init_group(pgrp grp, int group_size);
void 
free_group(pgrp grp);

/*
 *
 ****/
void
run_super_tree(char *file){
	int i, j, k;
	//firstly read in orders and init d_mat
	po o = (po)malloc(sizeof(to));	
	init_orders(file, o);
	pdism d_mat = (pdism)malloc(sizeof(tdism));
	init_dis_mat(d_mat, o->num_input_genome);
	compute_distance_matrix(o, d_mat); // we only need to compute distance once
	//print_dis_mat(d_mat);
	//compute the eigen value and do the partition
	pgrp group_id = (pgrp)malloc(sizeof(tgrp));
	init_group(group_id, o->num_input_genome);
	int num_group = 1;
	int group_now = 0;
	int has_large_disk = TRUE;
	double thresh = 0.1; 
	//alpha is used to decide if place overlapped genomes or not
	//beta is used to decide if genome is going to placed in overlapping region or not
	//alpha = 0.08 * sqrt(100/N), beta = 0.005 * sqrt(100/N)
	while((++group_now) <= num_group){
		double alpha = 0.08 * sqrt(100.0/(o->num_input_genome));
		double beta = 0.005 * sqrt(100.0/(o->num_input_genome));
		po oo = (po)malloc(sizeof(to));	
		int *map = (int*)malloc(sizeof(int)*o->num_input_genome);
		int success = init_orders_by_group(oo, o, group_now, map, group_id);
		//printf("group_now %d num in this group %d num_group %d\n",group_now, success, num_group);
		if(success >6){
			pdism dd_mat = (pdism)malloc(sizeof(tdism));
			init_dis_mat(dd_mat, oo->num_input_genome);
			compute_distance_matrix(oo, dd_mat); // we only need to compute distance once
			//print_dis_mat(dd_mat);
			compute_laplacian_mat(dd_mat);
			//print_q_mat(dd_mat);
			double *vector = (double*)malloc(sizeof(double)*dd_mat->num_spc);
			double gap = compute_second_eigen_vector(dd_mat, vector);	
			int num_pos=0, num_neg=0;
			compute_partition(vector, alpha, beta, &num_pos, &num_neg, dd_mat->num_spc, gap);
			//printf("gap %f alpha %f beta %f num_pos %d num_neg %d\n", gap, alpha, beta, num_pos, num_neg);
			int o_group_id = group_now;
			//parition
			for(j=0; j<dd_mat->num_spc; j++)
				group_id->belong[map[j]][o_group_id] = FALSE;
			for(j=0; j<num_neg; j++)
				group_id->belong[map[j]][o_group_id * 2 + 1] = TRUE;
			for(i=0, j=dd_mat->num_spc-1; i<num_pos; i++, j--)
				group_id->belong[map[j]][o_group_id*2] = TRUE;
			num_group += 2;
		}
		free_orders(oo);
		free(map);
	}
	//for(i=0; i<o->num_input_genome; i++){
	//	for(j=0; j<o->num_input_genome; j++){
	//		printf("%d ", group_id->belong[i][j]);	
	//	}
	//	printf("\n");
	//}
	//based on the partition compute tree for each partition
	int **map = (int**)malloc(sizeof(int*)*o->num_input_genome);
	pt *trees = (pt*)malloc(sizeof(pt)*(num_group+1));
	for(i=0; i<num_group+1; i++)
		trees[i] = (pt)malloc(sizeof(tt));
	int *group_size = (int*)malloc(sizeof(int)*(num_group+1));
	for(i=0; i<num_group+1; i++){
		group_size[i] = 0;
                for(j=0; j<o->num_input_genome; j++)
			if(group_id->belong[j][i]==TRUE)
				group_size[i]++;
	}
	for(i=0; i<num_group+1; i++){
		if(group_size[i]>3){
			po oo = (po)malloc(sizeof(to));	
			map[i] = (int*)malloc(sizeof(int)*o->num_input_genome);
			init_orders_by_group(oo, o, i, map[i], group_id);
			int upper_bound = run_upper_bound(oo);
			int lower_bound = run_lower_bound(oo);
			trees[i] = bnb_mp_tree(oo, upper_bound, lower_bound);
			free_orders(oo);
		}
	}
	//output to the file for supertree merge
	char f_buf[100];
	sprintf(f_buf, "%s_newick", file);
	FILE *writer = fopen(f_buf, "w");
	if(writer == NULL){
		printf("file %s does not exists!\n", file);
		ERROR_PRINT();
	}
	for(i=0; i<num_group+1; i++){
		if(group_size[i]<=3){
			if(group_size[i]==2){
				int to[2];
				for(j=0; j<2; j++){
					int c= j;
					for(k=c==0?0:to[0]+1; k<o->num_input_genome; k++)
						if(group_id->belong[k][i]==TRUE){
							to[c] = k;
							break;
						}
				}
				fprintf(writer, "(%d,%d);\n", to[0],to[1]);
			}
			if(group_size[i]==3){
				int to[3];
				for(j=0; j<3; j++){
					int c= j;
					for(k=c==0?0:to[j-1]+1;k<num_group;k++)
						if(group_id->belong[k][i]==TRUE){
							to[c] = k;
							break;
						}
				}
				fprintf(writer, "(%d,%d,%d);\n", to[0],to[1],to[2]);
			}
		}
		else
			from_edge_to_newick(writer, trees[i], map);
	}
	exit(1);

	//free stuff

	//int *merged = (int*)malloc(sizeof(int*)*num_group);
	//for(i=0; i<num_group; i++)
	//	merged[i] = FALSE;
	////merge trees that is non-overlapped
	//for(i=0; i<num_group; i++){
	//	for(j=0; j<num_group; j++){
	//		if(merged[i]==FALSE && merged[j]==FALSE){
	//			int success = merge_tree_separate(trees[i], trees[j], map[i], map[j], 
	//					trees[i]->orders->num_genome, trees[j]->orders->num_genome, d_mat);
	//			if(success == TRUE)
	//				merged[j] = TRUE;
	//		}
	//	}
	//}
	//int **rf_dis = (int**)malloc(sizeof(int*)*num_group); 
	//for(i=0; i<num_group; i++){
	//	rf_dis[i] = (int*)malloc(sizeof(int)*num_group);
	//}
	////merge trees that is overlapped
	////get separation for each tree
	//tsep **separations = (tsep**)malloc(sizeof(tsep*)*num_group);
	//for(i=0; i<num_group; i++){
	//	separations[i] = (tsep*)malloc(sizeof(tsep)*(trees[i]->orders->num_input_genome-3));
	//	int s_pos=0;
	//	for(j=0; j<trees[i]->num_edge; j++){
	//		int from = trees[i]->edges[j].from;
	//		int to = trees[i]->edges[j].to;
	//		//should be inner edges for separation
	//		if(from<trees[i]->orders->num_input_genome && to <trees[i]->orders->num_input_genome){
	//			create_bit_arr(separations[i][s_pos].bit_code, o->num_input_genome);
	//			create_bit_arr(separations[i][s_pos].bit_exist, o->num_input_genome);
	//			separations[i][s_pos].edge_id = j;
	//			s_pos++;
	//		}
	//	}
	//}
	//while(1){
	//	for(i=0; i<num_group; i++){
	//		if(merged[i]==FALSE){
	//			int s_pos=0;
	//			for(j=0; j<trees[i]->num_edge; j++){
	//				int from = trees[i]->edges[j].from;
	//				int to = trees[i]->edges[j].to;
	//				//should be inner edges for separation
	//				if(from<trees[i]->orders->num_input_genome && to <trees[i]->orders->num_input_genome){
	//					calculate_separation_fp(&separations[i][s_pos], trees[i], map[i], j);
	//					separations[i][s_pos].edge_id = j;
	//					s_pos++;
	//				}
	//			}
	//		}
	//		for(j=i+1; j<num_group; j++){
	//			if(merged[i]==FALSE && merged[j]==FALSE)
	//				rf_dis[i][j] = compute_rf_dis(separations[i], separations[j]);
	//			else 
	//				rf_dis[i][j] = 1000000;
	//		}
	//	}
	//	int min=100000, min_i, min_j;
	//	for(i=0; i<num_group; i++){
	//		for(j=i+1; j<num_group; j++){
	//			if(rf_dis[i][j]<min){
	//				min = rf_dis[i][j];
	//				min_i = i;
	//				min_j = j;
	//			}
	//		}
	//	}
	//	merge_tree_overlap(separations[min_i], separations[min_j], 
	//			trees[min_i], trees[min_j], 
	//			map[min_i], map[min_j], 
	//			trees[min_i]->num_edge, trees[min_j]->num_edge);
	//	merged[min_j] = TRUE;
	//	int un_merged_count = 0;
	//	for(i=0; i<num_group; i++)
	//		if(merged[i]==FALSE)
	//			un_merged_count++;
	//	if(un_merged_count == 1)
	//		break;
	//}
	//free_group(group_id);
}

void 
init_group(pgrp grp, int group_size){
	int i, j, k;
	grp->group_size = group_size;
	grp->belong = (int**)malloc(sizeof(int*)*group_size);
	for(i=0; i<group_size; i++){
		grp->belong[i] = (int*)malloc(sizeof(int)*group_size);
		for(j=0; j<group_size; j++){
			grp->belong[i][j] = FALSE;
		}
	}
	for(i=0; i<group_size; i++)
		grp->belong[i][1] = TRUE;
}

void 
free_group(pgrp grp){
	int i;
	for(i=0; i<grp->group_size; i++)
		free(grp->belong[i]);
	free(grp->belong);
	free(grp);
}

/*
 *
 ****/
int 
check_over_lap(pt tree_to, pt tree_from, 
		int *map_to, int *map_from, 
		int size_from, int size_to){
	int i, j, k;
	//get vertex set
	for(i=0; i<size_from; i++)
		for(j=0; j<size_to; j++)
			if(map_to[i] == map_from[j])
				return TRUE;
	return FALSE;
}

/*
 *
 ****/
int
merge_tree_overlap(tsep *sep_from, tsep *sep_to,
		pt tree_to, pt tree_from, 
		int *map_to, int *map_from, 
		int size_from, int size_to){
	int i, j, k;
	int is_overlaped = check_over_lap(tree_to, tree_from, 
			map_to, map_from, 
			size_from, size_to);
	if(is_overlaped == FALSE)
		return FALSE;
	//construct a consensus tree
	int *tmp_map = (int*)malloc(sizeof(size_from));
	po tmp_o = (po)malloc(sizeof(to));
	int num_input_from = tree_from->orders->num_input_genome;
	int num_input_to = tree_from->orders->num_input_genome;
	int **common_edge = (int**)malloc(sizeof(int*)*2);	
	common_edge[0] = (int*)malloc(sizeof(int)*tree_from->num_edge);
	common_edge[1] = (int*)malloc(sizeof(int)*tree_to->num_edge);
	//collapse edges of both genomes to make a consensus tree
	for(i=0; i<tree_from->num_edge; i++){
		if(tree_from->edges[i].from < num_input_from || tree_from->edges[i].to < num_input_from)
			common_edge[0][i] = TRUE;
		common_edge[0][i] = FALSE;
	}
	for(i=0; i<tree_to->num_edge; i++){
		if(tree_to->edges[i].from < num_input_to || tree_to->edges[i].to < num_input_to)
			common_edge[0][i] = TRUE;
		common_edge[1][i] = FALSE;
	}
	for(i=0; i<size_from; i++){
		int num = separation_exist(sep_to, sep_from, size_to, i);
		if(num==1)
			common_edge[0][sep_from[i].edge_id] = TRUE;
	}
	for(i=0; i<size_to; i++){
		int num = separation_exist(sep_from, sep_to, size_from, i);
		if(num==1)
			common_edge[1][sep_to[i].edge_id] = TRUE;
	}
	tte *edges = (tte*)malloc(sizeof(tte)*100);
	int edge_idx = 0;
	for(i=0; i<tree_to->num_edge; i++)
		if(common_edge[1][i]==TRUE)
			separate(edges, edge_idx, sep_from, i);
	//add un-overlapped edges	
	for(i=0; i<tree_from->num_edge; i++)
		if(common_edge[0][i]==LEAF_INDP)
			add_independent_edge(edges, edge_idx, i, i);
	//solve ambiguity
	
}

void
add_independent_edge(tte *edges, int edge_idx, int from, int to){
}

void
separate(tte *edges, int edge_idx, tsep *sep, int sep_idx){
}
/*
 *
 ****/
int
merge_tree_separate(pt tree_to, pt tree_from, 
		int *map_to, int *map_from, 
		int size_from, int size_to,
		pdism d_mat){
	int i, j, k;
	int is_overlaped = check_over_lap(tree_to, tree_from, 
			map_to, map_from, 
			size_from, size_to);
	if(is_overlaped != FALSE)
		return FALSE;
	int best_e_from, best_e_to;
	/*****this is in rec-dcm-eigen paper*/
	//get the min pairwise distance	
	int pair_from, pair_to, min_dis=100000000;		
	for(i=0; i<tree_from->num_edge; i++){
		for(j=0; j<tree_to->num_edge; j++){
			int dis = d_mat->dis_mat[map_from[i]][map_to[j]];
			if(dis < min_dis){
				min_dis = dis;
				pair_from = i;
				pair_to = j;
			}
		}
	}
	//get the initial upper bound
	int upper_bound = 100000000;
	for(i=0; i<tree_from->num_edge; i++){
		for(j=0; j<tree_to->num_edge; j++){
			if((tree_from->edges[i].from==pair_from || 
						tree_from->edges[i].to==pair_from) &&
					(tree_to->edges[i].from==pair_to ||
					 tree_to->edges[i].to==pair_to)){
				int dis = compute_sep_merge_dis(tree_from, tree_to, i, j, TYPE_DIS);
				if(dis < upper_bound){
					upper_bound = dis;
					best_e_from = i;
					best_e_to = j;
				}
			}
		}
	}
	//this part could be a paper in hicomb, any edge is possible and which is the improvement step
	for(i=0; i<tree_from->num_edge; i++){
		for(j=0; j<tree_to->num_edge; j++){
			int lower_bound = compute_sep_merge_dis(tree_from, tree_to, i, j, TYPE_BOUND);
			//there is no need to compute if the lower bound is larger the upper bound
			if(lower_bound<upper_bound){
				int dis = compute_sep_merge_dis(tree_from, tree_to, i, j, TYPE_DIS);                                
				if(dis < upper_bound){
					upper_bound = dis;
					best_e_from = i;
					best_e_to = j;
				}
			}	
		}
	}
	//finally merge the tree with the best edge pair
	//compute_sep_merge_dis(tree_from, tree_to, i, j, TYPE_BOUND);
	//compute_sep_merge_dis(tree_from, tree_to, i, j, TYPE_DIS);
	merge_sep_by_edges(tree_from, tree_to, best_e_from, best_e_to);
	return TRUE;
}

/*
 * given orders and their group id
 * we can initiate an order set by their assigend group id
 * if group id is set to 0 it's definitely selected 
 ****/
int
init_orders_by_group(po to, po from, 
		int group_num, int *map, 
		pgrp group_id){
	int i, j, k;
	int total_num = 0;
	//calculate how many elements in the group	
	for(i=0; i< from->num_input_genome; i++)
		if(group_id->belong[i][group_num] == TRUE)
			total_num++;	

	to->num_genome = 2*total_num -1;
	to->num_med_genome = total_num -1;
	to->num_input_genome = total_num;
	to->med_genome_idx = total_num;
	to->v_size = from->v_size;
	to->e_size = (int*)malloc(sizeof(int)*to->num_genome);
	to->valid = (int*)malloc(sizeof(int)*to->num_genome);
	to->v_idx = (int**)malloc(sizeof(int*)*to->num_genome);
	to->e_idx = (int**)malloc(sizeof(int*)*to->num_genome);
	for(i=0; i<to->num_genome; i++){
		to->v_idx[i] = (int*)malloc(sizeof(int)*to->v_size);
		to->e_size[i] = 0;
		to->valid[i] = FALSE;
	}
	int pos = 0;
	for(i=0; i<from->num_input_genome; i++){
		if(group_id->belong[i][group_num] != TRUE)
			continue;
		map[pos] = i;
		to->valid[pos] = TRUE;
		to->e_size[pos] = from->e_size[i];
		to->e_idx[pos] = (int*)malloc(sizeof(int)*to->e_size[pos]);
		for(j=0; j<to->e_size[pos]; j++)
			to->e_idx[pos][j] = from->e_idx[i][j]; 
		for(j=0; j<to->v_size; j++)
			to->v_idx[pos][j] = from->v_idx[i][j]; 
		pos++;
	}
	return total_num;
}

/*
 *
 ****/
void 
compute_laplacian_mat(pdism d_mat){
	int i, j, k;
	for(i=0; i<d_mat->num_spc; i++)
		d_mat->q_mat[i][i] = 0;
	for(i=0; i<d_mat->num_spc; i++)
		for(j=0; j<d_mat->num_spc; j++)
			d_mat->q_mat[i][i] += d_mat->dis_mat[i][j];
	for(i=0; i<d_mat->num_spc; i++)
		for(j=i+1; j<d_mat->num_spc; j++)
			d_mat->q_mat[i][j] = d_mat->q_mat[j][i] = d_mat->q_mat[i][j]-d_mat->dis_mat[i][j];
}

int 
compute_sep_merge_dis(pt tree_from, pt tree_to, int e_from, int e_to, int type){
	int i, j, k;
	int old_dis = tree_from->edges[e_from].score + tree_to->edges[e_to].score;
	//init graph
	pg g = (pg)malloc(sizeof(tg));
	g->v_size = tree_from->orders->v_size;
	g->v_idx = (int**)malloc(sizeof(int*)*4);
	g->v_deg = (int**)malloc(sizeof(int*)*4);
	g->v_check = (int*)malloc(sizeof(int)*4);
	g->e_idx = (int**)malloc(sizeof(int*)*4);
	g->e_size = (int*)malloc(sizeof(int)*4);
	for(i=0; i<4; i++){
		g->v_idx[i] = (int*)malloc(sizeof(int)*g->v_size);
		g->v_deg[i] = (int*)malloc(sizeof(int)*g->v_size);
	}
	init_order_by_edge(g, tree_from, e_from, TYPE_EDGE_M);
	init_order_by_edge(g, tree_to, e_to, TYPE_EDGE);
	rename_graph_by_car(g);
	init_median(g);
	if(type == TYPE_DIS)
		lk_opt(g, 2, 3, FALSE);
	int dis1 = compute_all_possible_exemplar_orders(tree_from->orders, tree_from->edges[e_from].from, 
			tree_from->orders->med_genome_idx-1, tree_from->orders->v_size/2);
	int dis2 = compute_all_possible_exemplar_orders(tree_from->orders, tree_from->edges[e_from].to, 
			tree_from->orders->med_genome_idx-1, tree_from->orders->v_size/2);
	int dis3 = evaluate(g);
	//resume the order
	tree_from->orders->med_genome_idx--;		
	free(tree_from->orders->e_idx[tree_from->orders->med_genome_idx]);
	int result = dis1+dis2+dis3-old_dis;
	return result;
}

void 
init_order_by_edge(pg g, pt tree, int e_id, int type){
	int i, j, k;
	if(type == TYPE_EDGE_M){
		//get the 	
		int from = tree->edges[e_id].from;
		int to = tree->edges[e_id].to;
		init_dup_car(tree->orders, from, to);
		g->e_idx[0] = (int*)malloc(sizeof(int)*tree->orders->e_size[tree->orders->med_genome_idx-1]);
		g->e_size[0] = tree->orders->e_size[tree->orders->med_genome_idx-1];
		//copy the v_idx and e_idx
		for(i=0; i<tree->orders->v_size; i++)
			g->v_idx[0][i] = tree->orders->v_idx[tree->orders->med_genome_idx-1][i];	
		for(i=0; i<tree->orders->e_size[tree->orders->med_genome_idx-1]; i++)
			g->e_idx[0][i] = tree->orders->e_idx[tree->orders->med_genome_idx-1][i];	
	}
	else if(type == TYPE_EDGE){
		int color;
		int from = tree->edges[e_id].from;
		int to = tree->edges[e_id].to;
		for(color=0; color<2; color++){
			int c=color==0?from:to;
			int cc=color==0?1:2;
			g->e_idx[cc] = (int*)malloc(sizeof(int)*tree->orders->e_size[c]);
			g->e_size[cc] = tree->orders->e_size[c];
			//copy the v_idx and e_idx
			for(i=0; i<tree->orders->v_size; i++)
				g->v_idx[cc][i] = tree->orders->v_idx[c][i]; 
			for(i=0; i<tree->orders->e_size[c]; i++)
				g->e_idx[cc][i] = tree->orders->e_idx[c][i];
		}
	}
}


void
merge_sep_by_edges(pt tree_from, pt tree_to, 
		int best_e_from, int best_e_to){
	int i, j, k, color;
	//get a tmp order
	pg g = (pg)malloc(sizeof(tg));
	po o_from = tree_from->orders; 
	po o_to = tree_to->orders;
	po tmp = (po)malloc(sizeof(to));	
	tmp->num_genome = o_from->num_genome+o_to->num_genome + 2;
        tmp->num_med_genome = o_from->num_med_genome+o_to->num_med_genome + 2;
        tmp->num_input_genome = o_from->num_input_genome+o_to->num_input_genome;
        tmp->med_genome_idx = o_from->num_genome+o_to->num_genome;
        tmp->v_size = o_from->v_size;
        tmp->e_size = (int*)malloc(sizeof(int)*tmp->num_genome);
        tmp->valid = (int*)malloc(sizeof(int)*tmp->num_genome);
        tmp->v_idx = (int**)malloc(sizeof(int*)*tmp->num_genome);
        tmp->e_idx = (int**)malloc(sizeof(int*)*tmp->num_genome);
	for(i=0; i<o_from->num_genome; i++){
		tmp->e_size[i] = o_from->e_size[i];
		tmp->valid[i] = o_from->valid[i];
		tmp->v_idx[i] = (int*)malloc(sizeof(int)*o_from->v_size);
		tmp->e_idx[i] = (int*)malloc(sizeof(int)*o_from->e_size[i]);
		for(j=0; j<o_from->v_size; j++)
			tmp->v_idx[i][j] = o_from->v_idx[i][j];
		for(j=0; j<o_from->e_size[i]; j++)
			tmp->e_idx[i][j] = o_from->e_idx[i][j];
	}
	for(i=0, k=o_from->num_genome; i<o_from->num_genome; i++, k++){
		tmp->e_size[k] = o_to->e_size[i];
		tmp->valid[k] = o_to->valid[i];
		tmp->v_idx[k] = (int*)malloc(sizeof(int)*o_to->v_size);
		tmp->e_idx[k] = (int*)malloc(sizeof(int)*o_to->e_size[i]);
		for(j=0; j<o_to->v_size; j++)
			tmp->v_idx[k][j] = o_to->v_idx[i][j];
		for(j=0; j<o_to->e_size[i]; j++)
			tmp->e_idx[k][j] = o_to->e_idx[i][j];
	}
	for(color=0, k=tmp->num_genome-2; color<2; color++, k++){
		int c= color==0?0:3;
		tmp->e_size[k] = g->e_size[c];
		tmp->valid[k] = TRUE;
		tmp->v_idx[k] = (int*)malloc(sizeof(int)*g->v_size);
		tmp->e_idx[k] = (int*)malloc(sizeof(int)*g->e_size[c]);
		for(j=0; j<g->v_size; j++)
			tmp->v_idx[k][j] = g->v_idx[c][j];
		for(j=0; j<g->e_size[c]; j++)
			tmp->e_idx[k][j] = g->e_idx[c][j];
		
	}
	//destruct the order of the tree_from and tree_to
	free_orders(o_from);
	//copy back
	o_from = tmp;
}

void
calculate_separation_fp(psep sep, pt t, 
		int *map, int e_id){
	int i, j, k;
	empty_bit_arr(sep->bit_code);
	empty_bit_arr(sep->bit_exist);
	//first compute the existance	
	for(i=0; i<t->orders->num_input_genome; i++)
		add_bit_arr_pos(sep->bit_exist, map[i]);
	
	//then compute the separation code, if something does not exist it's separation is 0
	int *queue = (int*)malloc(sizeof(int)*t->orders->num_input_genome);
	int *visited = (int*)malloc(sizeof(int)*t->orders->num_input_genome); 
	for(i=0; i<t->orders->num_input_genome; i++)
		visited[i] = FALSE;
	int queue_idx = 0;
	queue[queue_idx++] = t->edges[e_id].from;
	visited[queue_idx-1] = TRUE;
	while(queue_idx > 0){
		int frontier = queue[--queue_idx];
		add_bit_arr_pos(sep->bit_code, map[frontier]);
		for(i=0; i<t->num_edge; i++){
			int from = t->edges[i].from;
			int to = t->edges[i].to;
			if(from==frontier && visited[to]==FALSE){
				queue[queue_idx++] = to;
			}
			if(to==frontier && visited[from]==FALSE){
				queue[queue_idx++] = from;
			}
		}
	}
}

int
compute_rf_dis(tsep *from, tsep *to, int size_from, int size_to){
	int i, j, k;
	//from - to
	int num_same;
	for(i=0; i<size_to; i++){
		num_same += separation_exist(from, to, size_from, i);
	}
	int result = size_from + size_to - num_same;
}

int
separation_exist(tsep *from, tsep *to, int size_from, int idx_to){
	int i, j, k;
	int size = from[i].bit_code->size; 
	//for(i=0; i<size_from; i++){
	//	for(j=0; j<size; j++){
	//		if(CHECK_BIT(from[i].bit_exist->arr, j) && CHECK_BIT(to[idx_to].bit_exist->arr, j)){
	//			if((CHECK_BIT(from[i].bit_code, j) && !CHECK_BIT(to[idx_to].bit_code, j)) || 
	//					(!CHECK_BIT(from[i].bit_code, j) && CHECK_BIT(to[idx_to].bit_code, j)))
	//				return 0;
	//		}
	//	}
	//}
	return 1;
}

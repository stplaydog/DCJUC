#include "local.h"
#include "list.h"
#include <stdlib.h>

int min_deg(pg g, int pos);
int num_zero(pg g, int pos);
int median_deg(pg g, int pos);
int 
get_deg(pg g, int pos);
int 
find_next_edge(pg g, int *rev_v,
		 int from);
void 
add_child_max(pg g, psl l, int current_level, int term_level);
void 
add_child_combinations(pg g, psl l, int current_level, int use_heu);

int turn_on_debug = FALSE;

int DEBUG_RETRIEVE = FALSE;
int HAS_PROBLEM = FALSE;

extern int round_count;

/* ***************************************************************
 * this function is used to initialize the median adj in the graph
 * 1) choose the content of the median adj
 * 2) randomly assign the adj
 * ***************************************************************/
void init_median(pg g){
	int i,j,k;
	//first assigning spots for the csr
	//for(i=0; i<g->v_size; i++)
	//	if(g->zero[i] !=NUL)
	//		printf("%d %d\n", i, g->zero[i]);
	for(i=0;i<g->v_size;i++)
		g->v_deg[3][i] = NUL;
	for(i=0;i<g->v_size;i++){
		if(g->v_deg[3][i] == NUL){
			int deg = get_deg(g, i);
			g->v_deg[3][i] = deg;
		}
	}
	for(i=0;i<g->v_size;i++)
		if(g->v_deg[3][i] == NUL)
			g->v_deg[3][i] = 0;
	//for(i=0; i<g->v_size; i++)
	//	printf("%d:%d ", i, g->v_deg[3][i]);
	//printf("\n");
	//for(i=0;i<g->v_size; i++)
	//	printf("%d|%d|%d|%d ", g->v_deg[0][i], g->v_deg[1][i], g->v_deg[2][i], g->v_deg[3][i]);
	//printf("\n");
	//compute v_idx
	int sum = 0;
	for(i=0; i<g->v_size; i++){
		//printf("%d:%d:%d ", g->v_deg[0][i], g->v_deg[1][i], g->v_deg[2][i]);
		sum += g->v_deg[3][i];
		g->v_idx[3][i] = sum;
	}
	g->e_size[3]=sum;
	//printf("sum is %d\n", sum);
	//randomly init edges
	//srand(time(NULL));
	for(i=0;i<g->e_size[3]; i++)
		g->e_idx[3][i] = -1;
	//if(round_count==2){
	//	printf("actual %d\n", g->e_size[3]);
	//}
	//add zero matching
	//for(i=0;i<g->v_size;i++)
	//	if(g->v_check[i]==NUL){
	//		int start_zero = i==0?0:g->v_idx[3][i-1];
	//		int start_one = i==0?0:g->v_idx[0][i-1];
	//		int start_two = i==0?0:g->v_idx[1][i-1];
	//		int start_three = i==0?0:g->v_idx[2][i-1];
	//		int to_one = g->e_idx[0][start_one];
	//		int to_two = g->e_idx[0][start_two];
	//		int to_three = g->e_idx[0][start_three];
	//		if(to_one == to_two)
	//			g->e_idx[0][start_zero] = to_one;
	//		else if(to_three == to_two)
	//			g->e_idx[0][start_zero] = to_two;
	//		else if(to_one == to_three)
	//			g->e_idx[0][start_zero] = to_one;
	//		
	//	}
	int *rev_v = (int*)malloc(sizeof(int)*g->e_size[3]);
	for(i=0; i<g->v_size; i++){
		int start = i==0?0:g->v_idx[3][i-1];
		int end = g->v_idx[3][i];
		for(j=start; j<end; j++)
			rev_v[j] = i;
	}
	for(i=0; i<g->e_size[3]; i++){
		if(g->e_idx[3][i] == -1){
			int to = find_next_edge(g, rev_v, i);
			//printf("%d %d %d %d\n", rev_v[i], rev_v[to], i, to);
			if(to != CAP){
				g->e_idx[3][i] = rev_v[to];
				g->e_idx[3][to] = rev_v[i];
			}
			else
				g->e_idx[3][i] = CAP;
		}
	}
	//graph_vis_two(g, "after.dot", 0, 3);
	//graph_vis_two(g, "init1.dot", 1, 3);
	//graph_vis_two(g, "init2.dot", 2, 3);
	//exit(1);

	//for(i=0; i<g->v_size; i++){
	//	int start = i==0?0:g->v_idx[3][i-1];
	//	int end = g->v_idx[3][i];
	//	for(j=start; j<end; j++){
	//		if(g->e_idx[3][j] == -1){
	//			int to;
	//			int to_e;
	//			while(1){
	//				to = rand() % g->v_size;
	//				if(g->v_check[to]!=NUL && to != i){
	//					int start = to==0?0:g->v_idx[3][to-1];
	//					int end = g->v_idx[3][to];
	//					int k;
	//					int found = FALSE;
	//					for(k=start;k<end;k++)
	//						if(g->e_idx[3][k] == -1){
	//							found = TRUE;
	//							to_e = k;
	//						}
	//					if(found == TRUE)
	//						break;
	//				}
	//			}
	//			g->e_idx[3][j] = to;
	//			g->e_idx[3][to_e] = i;
	//		}
	//	}
	//}
	free(rev_v);
}

int 
find_next_edge(pg g, int *rev_v,
		int from){
	int i;
	int conflict = 0;
	if(rev_v[from] % 2 == 0)
		conflict = rev_v[from]+1;
	else
		conflict = rev_v[from]-1;
	int result = -1;
	for(i=0; i<g->e_size[3]; i++)
		if(g->e_idx[3][i] == -1 && rev_v[i] != conflict 
				&& rev_v[i] != rev_v[from]){
			result = i;
			break;
		}
	if(result != -1)
		return result;
	return CAP;
}

int 
get_deg(pg g, int pos){
	int d1 = g->v_deg[0][pos];	
	int d2 = g->v_deg[1][pos];	
	int d3 = g->v_deg[2][pos];	
	if(d1==d2 && d1==d3)
		return d1;
	else{
		if(max_deg(g, pos)==1)
			return 1;
		if(max_deg(g, pos)>1)
			return max_deg(g, pos)/2;
		if(max_deg(g, pos)==0)
			return 0;
	}
	return 0;
}


int max_deg(pg g, int pos){
	int one = g->v_deg[0][pos];
	int two = g->v_deg[1][pos];
	int three = g->v_deg[2][pos];
	if(one>two && one>three)
		return one;
	else if(two > three)
		return two;
	else
		return three;
}

/* *********************************************************
 * given a vertex, find the vertex which has the min degree
 * *********************************************************/
int min_deg(pg g, int pos){
	int one = g->v_deg[0][pos];
	int two = g->v_deg[1][pos];
	int three = g->v_deg[2][pos];
	if(one<=two && one<=three)
		return one;
	else if(two <= three)
		return two;
	else
		return three;
}

/* ****************************************************
 * find the number of vertices that have degree of zero
 * ****************************************************/
int num_zero(pg g, int pos){
	int sum = g->v_deg[0][pos]==0?0:1 + g->v_deg[1][pos]==0?0:1 + g->v_deg[2][pos]==0?0:1;
	return sum;
}

/* ****************************************
 * find the median degree of three vertices
 * ****************************************/
int median_deg(pg g, int pos){
	int one = g->v_deg[0][pos];
	int two = g->v_deg[1][pos];
	int three = g->v_deg[2][pos];
	if(one>=two && one <= three || one>=three && one <= two)
		return one;
	else if(two>=three && two <= one || two>=one && two <= three)
		return two;
	else if(three>=two && three <= one || three>=one && three <= two)
		return three;
	return -1;
}

/******************************************************************************
 *heau level is the level before which all of search frontier will be expanded*
 *is_opt is to indicate if perform k-opt algorithm or not k=heu_level**********
 *******************************************************************************/
void lk_opt(pg g, int heu_level, 
		int term_move, int is_opt){
	int i,j,k,m,n;
	//init the list
	psl l = (psl)malloc(sizeof(tsl));	
	int improved=TRUE;
	int accept=FALSE;	
	int best_score = evaluate(g);
	printf("the initial score is %d\n", best_score);
	//upperbound(best_score) is the current score
	//lowerbound is (d_12 + d_23 + d_13)/2
	int buck_size = term_move;
	int list_size = 1000;
	int parent_size = g->e_size[3];
	create_search_list(buck_size, list_size, parent_size, l);
	int count =0;
	//FILE *dbg = fopen("dbg", "w"); 
	int use_heu = TRUE;
	graph_vis_two(g, "before.dot", 0, 3);
	while(improved==TRUE){
		int current_level=0;
		int processed=0;
		accept=FALSE;
		//printf("start another round of optimization\n");
		count = 0;
		improved = FALSE;
		while(accept==FALSE && current_level<term_move){
			if(current_level==0){
				add_child_combinations(g, l, current_level, use_heu);
				if(get_list_size(l) == 0)
					accept=TRUE;
				current_level = 1;
			}
			else{
				//if(get_list_size(l) == 267524+1){
				//	DEBUG_RETRIEVE = TRUE;
				//}	
				int has_remain_in_last_level = retrive(g, l, current_level-1);
				int start_0 = g->v_idx[3][0];
				int start_4 = g->v_idx[3][3];
				if(has_remain_in_last_level == FALSE){
					current_level += 1;
					printf("searching on level %d there are %d elements\n", current_level, l->idx_child[current_level-1]);
				}
				else{
					int score = evaluate(g);
					//fprintf(dbg, "size0 %d size1 %d\n", l->list_size[0], l->list_size[1]);
					if(score < best_score){
						accept=TRUE;
						best_score = score;
						improved = TRUE;
						reset_list(l);
						//printf("distance is %d best %d\n", score, best_score);
					}else{ //still need to be improved
						if(is_opt==FALSE && current_level < heu_level){
							add_child_combinations(g, l, current_level, use_heu);
						}
						else if(is_opt==FALSE && current_level >= heu_level){
							add_child_max(g, l, current_level, term_move);
						}
						//the following part is for k-opt algorithm
						else if(is_opt==TRUE && current_level < heu_level){
							add_child_combinations(g, l, current_level, use_heu);
						}
						else if(is_opt==TRUE && current_level >= heu_level){
							;
						}
					}
				}
			}
		}
	}
	//fclose(dbg);
	printf("the best score is %d\n", best_score);
	free_list(l);
}
/* **********************************************
 * add all possible switch operations to the list
 * **********************************************/
void 
add_child_combinations(pg g, psl l, int current_level, int use_heu){
	int i,j,m,n;
	int color;
	int idx_parent = add_parent(g, l, current_level);
	int **comp_id;
	int **comp_type;
	int thresh = 2;
	if(use_heu == TRUE){
		comp_id = (int**)malloc(sizeof(int*)*3);
		comp_type = (int**)malloc(sizeof(int*)*3);
		for(i=0; i<3; i++){
			comp_id[i] = (int*)malloc(sizeof(int)*g->v_size);
			comp_type[i] = (int*)malloc(sizeof(int)*g->v_size);
		}
		indentify_comps(comp_id, comp_type, g);
	}
	for(i=0;i<g->v_size;i++){
		if(g->v_check[i]==NUL)
			continue;
		int start_one = i==0?0:g->v_idx[3][i-1];
		int end_one = g->v_idx[3][i];
		for(j=start_one;j<end_one;j++){
			int from_one = i;
			int to_one = g->e_idx[3][j];
			for(m=0;m<g->v_size;m++){
				if(g->v_check[m]==NUL)
					continue;
				int start_two = m==0?0:g->v_idx[3][m-1];
				int end_two = g->v_idx[3][m];
				for(n= start_two; n<end_two; n++){
					int from_two = m;
					int to_two = g->e_idx[3][n];	
					if(from_one != to_one && from_one != from_two && from_one != to_two &&
							to_one != from_two && to_one != to_two && from_two != to_two
							&& from_one != CAP && to_one != CAP && from_two != CAP && to_two != CAP){ //this is the rule to chose whom to exchange
						if(use_heu==TRUE){
							int score = 0;
							for(color=0; color<3; color++){
								if(comp_id[color][i] == comp_id[color][m]){
									if(comp_type[color][comp_id[color][i]]==TYPE_SING)
										score+=2;
									if(comp_type[color][comp_id[color][i]]==TYPE_MUL)
										score+=1;
								}
							}
							if(score > thresh){
								add_child(from_one, to_one, from_two, to_two, CONN_TYPE_ONE, idx_parent, current_level, l);
								add_child(from_one, to_one, from_two, to_two, CONN_TYPE_TWO, idx_parent, current_level, l);
							}
						}
						else{
							add_child(from_one, to_one, from_two, to_two, CONN_TYPE_ONE, idx_parent, current_level, l);
							add_child(from_one, to_one, from_two, to_two, CONN_TYPE_TWO, idx_parent, current_level, l);
						}
					}
				}
			}
		}
	}
	if(use_heu == TRUE){
		for(i=0; i<3; i++){
			free(comp_id[i]);
			free(comp_type[i]);
		}
		free(comp_id);
		free(comp_type);
	}
}

/* *****************************
 * add only the best to the list
 * *****************************/

void 
add_child_max(pg g, psl l, int current_level, int term_level){
	int i, j, m, n;
#ifdef USE_APPR
	int best_score = 0;
#else
	int best_score = 100000000;
#endif
	int res_from_one=-1;
	int res_to_one=-1;
	int res_from_two=-1;
	int res_to_two=-1;
	int res_conn_type=-1;
	int color;

	int **comp_id = (int**)malloc(sizeof(int*)*3);
	int **comp_type = (int**)malloc(sizeof(int*)*3);
	for(i=0; i<3; i++){
		comp_id[i] = (int*)malloc(sizeof(int)*g->v_size);
		comp_type[i] = (int*)malloc(sizeof(int)*g->v_size);
	}
	indentify_comps(comp_id, comp_type, g);
	//for(i=0;i<3;i++){
	//	for(j=0;j<g->v_size; j++)
	//		printf("%d:%d:%d ", j, comp_id[i][j], comp_type[i][j]);
	//	printf("\n");
	//}
	//FILE *writer = fopen("comp.dot", "w");
	//fprintf(writer, "Graph{\n");
	//for(i=0; i<g->v_size; i++){
	//	if(g->v_check[i] == NUL)
	//		continue;
	//	int start = i==0?0:g->v_idx[0][i-1];
	//	int end = g->v_idx[0][i];
	//	for(j=start; j<end; j++)
	//		fprintf(writer, "%d:%d -- %d:%d [color=%s];\n", i,comp_id[0][i], g->e_idx[0][j], comp_id[0][g->e_idx[0][j]], "red");
	//	start = i==0?0:g->v_idx[3][i-1];
	//	end = g->v_idx[3][i];
	//	for(j=start; j<end; j++)
	//		fprintf(writer, "%d:%d -- %d:%d [color=%s];\n", i, comp_id[0][i], g->e_idx[3][j], comp_id[0][g->e_idx[3][j]], "blue");
	//}
	//fprintf(writer, "}\n");
	//fclose(writer);



	for(i=0;i<g->v_size;i++){
		if(g->v_check[i]==NUL)
			continue;
		int start_one = i==0?0:g->v_idx[3][i-1];
		int end_one = g->v_idx[3][i];
		for(j=start_one;j<end_one;j++){
			int from_one = i;
			int to_one = g->e_idx[3][j];
			for(m=0;m<g->v_size;m++){
				if(g->v_check[m]==NUL)
					continue;
				int start_two = m==0?0:g->v_idx[3][m-1];
				int end_two = g->v_idx[3][m];
				for(n= start_two; n<end_two; n++){
					int from_two = m;
					int to_two = g->e_idx[3][n];	
					if(from_one != to_one && from_one != from_two && from_one != to_two &&
							to_one != from_two && to_one != to_two && from_two != to_two
							&& from_one != CAP && to_one != CAP && from_two != CAP && to_two != CAP //this is the rule to chose whom to exchange
							&& from_one != NUL && to_one != NUL && from_two != NUL && to_two != NUL){ //this is the rule to chose whom to exchange
#ifdef USE_APPR
						int score = 0;
						for(color=0; color<3; color++){
							if(comp_id[color][i] == comp_id[color][m]){
								if(comp_type[color][comp_id[color][i]]==TYPE_SING)
									score+=2;
								if(comp_type[color][comp_id[color][i]]==TYPE_MUL)
									score+=1;
							}
						}
						if(score > best_score){ //select conn type
							res_from_one = from_one;
							res_to_one = to_one;
							res_from_two = from_two;
							res_to_two = to_two;
							//compute two dis
							shrink(g, from_one, to_one, from_two, to_two, CONN_TYPE_ONE);
							int score1 = evaluate(g);
							resume(g, from_one, to_one, from_two, to_two, CONN_TYPE_ONE);
							shrink(g, from_one, to_one, from_two, to_two, CONN_TYPE_TWO);
							int score2 = evaluate(g);
							resume(g, from_one, to_one, from_two, to_two, CONN_TYPE_TWO);
							//select conn type
							if(score1 < score2)
								res_conn_type = CONN_TYPE_ONE;
							else
								res_conn_type = CONN_TYPE_TWO;
							//printf("%d %d %d %d %d %d\n", res_from_one, res_to_one, res_from_two, res_to_two, res_conn_type, current_level);
							best_score = score;
						} 
#else
						shrink(g, from_one, to_one, from_two, to_two, CONN_TYPE_ONE);
						int score = evaluate(g);
						if(score < best_score){
							res_from_one = from_one;
							res_to_one = to_one;
							res_from_two = from_two;
							res_to_two = to_two;
							res_conn_type = CONN_TYPE_ONE;
						}
						resume(g, from_one, to_one, from_two, to_two, CONN_TYPE_ONE);
						shrink(g, from_one, to_one, from_two, to_two, CONN_TYPE_TWO);
						score = evaluate(g);
						if(score < best_score){
							res_from_one = from_one;
							res_to_one = to_one;
							res_from_two = from_two;
							res_to_two = to_two;
							res_conn_type = CONN_TYPE_TWO;
						}
						resume(g, from_one, to_one, from_two, to_two, CONN_TYPE_TWO);
#endif
					}
				}
			}
		}
	}
	//add
	if(res_conn_type != -1){
		int idx_parent = add_parent(g, l, current_level);
		add_child(res_from_one, res_to_one, res_from_two, res_to_two, res_conn_type, idx_parent, current_level, l);
	}
	for(i=0; i<3; i++){
		free(comp_id[i]);
		free(comp_type[i]);
	}
	free(comp_id);
	free(comp_type);
}

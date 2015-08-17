#include "tree.h"


void 
init_orders(char *file, po o);
void 
read_orders(char *file, po o);
void 
scan_orders(char *file, po o);
pt 
bnb_mp_tree(po o, 
		int upper_bound, int lower_bound);
void
init_graph_by_three_orders(po o, pg g, 
		int genome_from, int genome_to, int genome_avail);
void 
add_edge(pt tree, int which_edge);
void
copy_median_adj_to_order(po o, pg g);
void 
create_tree_list(ptl l, int upper_bound, int lower_bound);
void 
copy_tree(pt tree_from, pt tree_to);
void
copy_edge(tte *edges_to, tte *edges_from, int size);
void 
double_tree_list(ptl l, int pos);
void 
add_tree_list(ptl l, pt tree, int pos);
int 
pop_tree_list(ptl l, pt tree, int pos);
void
free_tree(pt tree);
void 
init_from_orders_to_tree(po o, pt tree,
		int one, int two, int three);
void 
cal_num_edge_operations(pt tree, int edge_idx, 
		int c1, int c2);
int 
max(int one, int two, int three);
void 
graph_vis_orders_two(po  o, char *file, int c1, int c2);
void 
init_best_tree(pt best, int num_genome, int num_edge);
void 
free_best_tree(pt best);
void 
copy_best_tree(pt best, pt tree_from);
void 
vis_best_tree(pt best);

int round_count;

/******************************************************************************************
 * this is the main function to compute the branch and bound based maximum parsimonious tree,
 * which is based on the median computation. the upperbound and lowerbound are coming from
 * the results of the neighbor joining.
 * ****************************************************************************************/
pt 
bnb_mp_tree(po o, 
		int upper_bound, int lower_bound){
	int i, j, k;
	//read into the order
	pt tree = (pt)malloc(sizeof(tt));	
	init_from_orders_to_tree(o, tree, 0, 1, 2); 
	pt best = (pt)malloc(sizeof(tt));	
	init_best_tree(best, tree->orders->num_genome, tree->num_edge);
	upper_bound = tree->upper_bound;
	lower_bound = tree->lower_bound;
	ptl l = (ptl)malloc(sizeof(ttl));
	create_tree_list(l, upper_bound, lower_bound);
	add_tree_list(l, tree, tree->lower_bound);
	int found = FALSE;
	round_count = 0;
	
	while(upper_bound > lower_bound && found==FALSE){
		pt tmp_tree = (pt)malloc(sizeof(tt));
		if(pop_tree_list(l, tmp_tree, lower_bound)==FALSE){
			lower_bound++;
			continue;
		}
		int i;
		round_count++;
		printf("searching ub %d lb %d\n", upper_bound, lower_bound);
		//printf("num_edge %d\n", tmp_tree->num_edge);
		//for(i=0; i<tmp_tree->num_edge; i++)
		//	if(tmp_tree->edges[i].from != NUL)
		//		printf("from %d to %d num_indel %d num_dup %d num_dcj %d score %d\n", 
		//			tree->edges[i].from, tree->edges[i].to, 
		//			tree->edges[i].num_indel, tree->edges[i].num_dup,
		//			tree->edges[i].num_dcj, tree->edges[i].score);
		//could be parallelized, but need some change in the data structure
		for(i=0; i<tmp_tree->num_edge; i++){
			if(tmp_tree->edges[i].from == NUL)
				continue;
			pt neighbor_tree = (pt)malloc(sizeof(tt));
			copy_tree(tmp_tree, neighbor_tree);
			printf("======start adding edge %d===========\n", i);
			add_edge(neighbor_tree, i);
			printf("======end adding edge %d ub %d:%d lb %d:%d===========\n", i, neighbor_tree->upper_bound, upper_bound, neighbor_tree->lower_bound, lower_bound);
			if(neighbor_tree->lower_bound > upper_bound)
				free_tree(neighbor_tree);
			else if(neighbor_tree->upper_bound < upper_bound){
				upper_bound = neighbor_tree->upper_bound;
				copy_best_tree(best, neighbor_tree);
				add_tree_list(l, neighbor_tree, neighbor_tree->lower_bound);	
			}
			else if(neighbor_tree->upper_bound <= lower_bound){
				free_tree(neighbor_tree);
				found = TRUE;
				break;
			}
				
			else
				add_tree_list(l, neighbor_tree, neighbor_tree->lower_bound);
		}
		free_tree(tmp_tree);
	}
	vis_best_tree(best);
	printf("best tree score is %d\n", best->score);
	//free_best_tree(best);
	printf("the best tree has been found, which has a score of %d %d\n", upper_bound, lower_bound);
	return best;
}

/*
 *
 * **/
void 
init_orders(char *file, po o){
	scan_orders(file, o);
	read_orders(file, o);
}


void 
free_orders(po o){
	int i;
	for(i=0; i<o->num_genome; i++){
		free(o->v_idx[i]);
		if(o->e_size[i] > 0){
			free(o->e_idx[i]);
		}
	}
	free(o->e_size);
	free(o->valid);
	free(o->v_idx);
	free(o->e_idx);
	free(o);
}

/* ********************************************************************************
 * this is to scan the graph file to get the basic information about the graph
 * and help to initialize the data structures.
 * ********************************************************************************/
void 
scan_orders(char *file, po o){
	FILE *stream;
	char str[100];
	int v_num;
	int e_num;
	int g_num;
	int count=0;

	int v_id=0;
	int v_to=0;
	int color=0;
	int i,j,sum=0;
	int e_num_read = 0;

	if((stream=fopen(file,"r"))==NULL){
		printf("the file %s you input does not exist!\n", file);
		exit(1);
	}
	if(fscanf(stream, "%d %d %d\n", &v_num, &g_num, &e_num)==EOF)
		printf("error here\n");
	o->v_size = v_num;
	o->num_genome = 2*g_num - 1;
        o->num_med_genome = g_num - 1;
        o->num_input_genome = g_num;
        o->med_genome_idx = g_num;

	o->v_idx = (int**)malloc(sizeof(int*)*o->num_genome);
	o->e_idx = (int**)malloc(sizeof(int*)*o->num_genome);
	o->e_size = (int*)malloc(sizeof(int)*o->num_genome);
	o->valid = (int*)malloc(sizeof(int)*o->num_genome);
	
	for(i=0; i<o->num_genome; i++){
		o->v_idx[i] = (int*)malloc(sizeof(int)*o->v_size);
		for(j=0; j<o->v_size; j++)
			o->v_idx[i][j] = 0;
		o->e_size[i] = 0;
		if(i < o->num_input_genome)
			o->valid[i] = TRUE;
		else
			o->valid[i] = FALSE;
	}
	for(i=0;i<e_num;i++){
		if(fscanf(stream, "%d %d %d\n",&v_id , &v_to, &color)==EOF)
			printf("error here\n");
		o->v_idx[color-1][v_id] +=1;
	}
	//for(i=0; i<g_num; i++){
	//	for(j=0.j<o->v_size; j++)
	//		printf("%d ", o->v_idx[i][j])
	//}
	for(i=0;i<g_num;i++){
		sum=0;
		//printf("sum %d\n", sum);
		for(j=0;j<v_num;j++){
			sum += o->v_idx[i][j];
			o->v_idx[i][j]=sum;
		}
		o->e_size[i] = o->v_idx[i][v_num-1];
		o->e_idx[i] = (int*)malloc(sizeof(int)*o->v_idx[i][v_num-1]);
	}
	fclose(stream);
}

/* *************************************************************************
 * based on the scanning result, this is the real function to read the real
 * graph in to orders.
 * *************************************************************************/
void 
read_orders(char *file, po o){
	FILE *stream;
	char str[100];
	int v_num;
	int e_num;
	int g_num;
	int count=0;

	int v_id=0;
	int v_to=0;
	int e_weight=0;
	int i,j,sum=0;
	int e_num_read = 0;

	if((stream=fopen(file,"r"))==NULL){
		printf("the file %s you input does not exist!\n", file);
		exit(1);
	}
	if(fscanf(stream, "%d %d %d\n", &v_num, &g_num, &e_num)==EOF)
		printf("error here\n");
	//this is for the start of the different postitions
	int **idx = (int**)malloc(sizeof(int*)*g_num);	
	for(i=0;i<g_num;i++){
		idx[i]=(int*)malloc(sizeof(int)*o->v_size);
		for(j=0;j<o->v_size;j++){
			idx[i][j] = j==0?0:o->v_idx[i][j-1];
		}
	}

	for(i=0;i<e_num;i++)
	{
		if(fscanf(stream, "%d %d %d\n",&v_id , &v_to, &e_weight)==EOF)
			printf("error here\n");
		int pos = idx[e_weight-1][v_id];
		o->e_idx[e_weight-1][pos] = v_to;
		idx[e_weight-1][v_id]++;
	}
	//finalize
	for(i=0; i<g_num; i++)
		free(idx[i]);
	free(idx);
	fclose(stream);
}

/* ********************************************************
 * initialize the bnb search list
 * ********************************************************/
void 
create_tree_list(ptl l, int upper_bound, int lower_bound){
	int i, j;
	l->base = lower_bound;
	int size = upper_bound-lower_bound+1;
	l->size_list = (int*)malloc(sizeof(int)*(size));
	l->idx_list = (int*)malloc(sizeof(int)*(size));
	l->list = (pt**)malloc(sizeof(pt*)*(size));
	for(i=0;i<size; i++){
		l->size_list[i] = LIST_SIZE;
		l->idx_list[i] = 0;
		l->list[i] = (pt*)malloc(sizeof(pt)*LIST_SIZE);
		for(j=0; j<LIST_SIZE; j++)
			l->list[i][j] = (pt)malloc(sizeof(tt));
	}
}

void
init_null_tree(po o, pt tree){
	int i,j,k;
	//number of edges |E|=(num_input_genome-3)*2 + 3
	//initialize the tree
	tree->num_edge = o->num_input_genome*2-1;
	tree->score = 0;
	tree->edge_idx = 0;
	tree->upper_bound = 0;
	tree->lower_bound = 0;
	
	tree->edges = (tte*)malloc(sizeof(tte)*tree->num_edge);
	for(i=0; i<tree->num_edge; i++){
		tree->edges[i].from = NUL;
		tree->edges[i].to = NUL;
		tree->edges[i].num_indel = 0;
		tree->edges[i].num_dup = 0;
		tree->edges[i].num_dcj = 0;
		tree->edges[i].score = 0;
	}
	tree->orders = (po)malloc(sizeof(to));
	//initiaize orders	
	tree->orders->num_genome = o->num_genome;
        tree->orders->num_med_genome = o->num_med_genome;
        tree->orders->num_input_genome = o->num_input_genome;
        tree->orders->med_genome_idx = o->med_genome_idx;
        tree->orders->v_size = o->v_size;

        tree->orders->e_size = (int*)malloc(sizeof(int)*tree->orders->num_genome);
        tree->orders->valid = (int*)malloc(sizeof(int)*tree->orders->num_genome);
        tree->orders->v_idx = (int**)malloc(sizeof(int*)*tree->orders->num_genome);
        tree->orders->e_idx = (int**)malloc(sizeof(int*)*tree->orders->num_genome);
	for(i=0;i<tree->orders->num_genome; i++){
		tree->orders->e_size[i] = o->e_size[i];
		tree->orders->valid[i] = o->valid[i];
		tree->orders->v_idx[i] = (int*)malloc(sizeof(int)*o->v_size);
		if(o->e_size[i] <= 0)
			continue;
		tree->orders->e_idx[i] = (int*)malloc(sizeof(int)*o->e_size[i]);
		for(j=0;j<tree->orders->v_size; j++)
			tree->orders->v_idx[i][j] = o->v_idx[i][j];
		for(j=0;j<tree->orders->e_size[i]; j++)
			tree->orders->e_idx[i][j] = o->e_idx[i][j];
	}
}

/* *******************************************************************************
 * init root node tree of the bnb algorithm
 * *******************************************************************************/
void 
init_from_orders_to_tree(po o, pt tree,
		int one, int two, int three)
{
	int i,j,k;
	//number of edges |E|=(num_input_genome-3)*2 + 3
	//initialize the tree
	tree->num_edge = o->num_input_genome*2-1;
	tree->score = 0;
	tree->edge_idx = 0;
	tree->upper_bound = 0;
	tree->lower_bound = 0;
	
	tree->edges = (tte*)malloc(sizeof(tte)*tree->num_edge);
	for(i=0; i<tree->num_edge; i++){
		tree->edges[i].from = NUL;
		tree->edges[i].to = NUL;
		tree->edges[i].num_indel = 0;
		tree->edges[i].num_dup = 0;
		tree->edges[i].num_dcj = 0;
		tree->edges[i].score = 0;
	}
	tree->orders = (po)malloc(sizeof(to));
	//initiaize orders	
	tree->orders->num_genome = o->num_genome;
        tree->orders->num_med_genome = o->num_med_genome;
        tree->orders->num_input_genome = o->num_input_genome;
        tree->orders->med_genome_idx = o->med_genome_idx;
        tree->orders->v_size = o->v_size;

        tree->orders->e_size = (int*)malloc(sizeof(int)*tree->orders->num_genome);
        tree->orders->valid = (int*)malloc(sizeof(int)*tree->orders->num_genome);
        tree->orders->v_idx = (int**)malloc(sizeof(int*)*tree->orders->num_genome);
        tree->orders->e_idx = (int**)malloc(sizeof(int*)*tree->orders->num_genome);
	for(i=0;i<tree->orders->num_genome; i++){
		tree->orders->e_size[i] = o->e_size[i];
		tree->orders->valid[i] = o->valid[i];
		tree->orders->v_idx[i] = (int*)malloc(sizeof(int)*o->v_size);
		if(o->e_size[i] <= 0)
			continue;
		tree->orders->e_idx[i] = (int*)malloc(sizeof(int)*o->e_size[i]);
		for(j=0;j<tree->orders->v_size; j++)
			tree->orders->v_idx[i][j] = o->v_idx[i][j];
		for(j=0;j<tree->orders->e_size[i]; j++)
			tree->orders->e_idx[i][j] = o->e_idx[i][j];
	}
	//init graph and compute median
	int genome_from = 0;
	int genome_to = 1;
	int genome_avail = 2;
	pg g = (pg)malloc(sizeof(tg));
	po oo = tree->orders;
	init_graph_by_three_orders(oo, g, genome_from, genome_to, genome_avail);
	//*****************get the median genome***************
	rename_graph_by_car(g);
	init_median(g);
	//graph_vis_two(g, "tree.dot", 0, 3);
	lk_opt(g, 2, 3, FALSE);
	//copy median graph back
	copy_median_adj_to_order(oo, g);
	//set up edges
	int edge_from_med = tree->edge_idx++;
	cal_num_edge_operations(tree, edge_from_med, 0, oo->med_genome_idx-1);

	int edge_to_med = tree->edge_idx++; 
	cal_num_edge_operations(tree, edge_to_med, 1, oo->med_genome_idx-1);

	int edge_avail_med = tree->edge_idx++; 
	cal_num_edge_operations(tree, edge_avail_med, 2, oo->med_genome_idx-1);

	//for(i=0; i<tree->edge_idx; i++)
	//	printf("from %d to %d num_indel %d num_dup %d num_dcj %d score %d\n", 
	//		tree->edges[i].from, tree->edges[i].to, tree->edges[i].num_indel,
	//		tree->edges[i].num_dup, tree->edges[i].num_dcj, tree->edges[i].score);
	//get the upper and lower_bound
	oo->valid[one] = FALSE;
	oo->valid[two] = FALSE;
	oo->valid[three] = FALSE;
	tree->score = tree->edges[edge_from_med].score + tree->edges[edge_to_med].score + tree->edges[edge_avail_med].score;
	tree->lower_bound = tree->score + run_lower_bound(oo);
	tree->upper_bound = tree->score + run_upper_bound(oo);
	free_graph(g);
}


/* *******************************************************************************
 * find an edge, then find a available input genome and compute the median for 
 * the two genomes on the two side of the edge and the selected genome.
 * *******************************************************************************/
void 
add_edge(pt tree, int which_edge){
	int i;
	//initiate the graph and copy the two end of the adj of the edge to the graph
	po o = tree->orders;
	pg g = (pg)malloc(sizeof(tg)); 
	int genome_from = tree->edges[which_edge].from;
	int genome_to = tree->edges[which_edge].to;
	int genome_avail=-1;
	for(i=0; i<o->num_input_genome; i++)
		if(o->valid[i]==TRUE){
			genome_avail = i;
			break;
		}
	if(genome_avail==-1){
		free(g);
		return;
	}
	init_graph_by_three_orders(o, g, genome_from, genome_to, genome_avail);
	printf("from %d to %d avail %d\n", genome_from, genome_to, genome_avail);
	o->valid[i] = FALSE;	
	//*****************get the median genome***************
	rename_graph_by_car(g);
	init_median(g);
	//graph_vis_two(g, "leak1.dot", 0, 3);
	lk_opt(g, 2, 3, FALSE);
	//copy median adj back to gene orders add new edge etc
	copy_median_adj_to_order(o, g);
	//change the edges
	int edge_from_med = which_edge;
	cal_num_edge_operations(tree, edge_from_med, genome_from, o->med_genome_idx-1);
	int edge_to_med = tree->edge_idx++; 
	cal_num_edge_operations(tree, edge_to_med, genome_to, o->med_genome_idx-1);
	int edge_avail_med = tree->edge_idx++; 
	cal_num_edge_operations(tree, edge_avail_med, genome_avail, o->med_genome_idx-1);
	//get the upper and lower_bound
	o->valid[genome_from] = FALSE;
	o->valid[genome_to] = FALSE;
	o->valid[genome_avail] = FALSE;
	tree->score += tree->edges[edge_from_med].score + tree->edges[edge_to_med].score + tree->edges[edge_avail_med].score;
	tree->upper_bound = tree->score + run_upper_bound(o);
	tree->lower_bound = tree->score + run_lower_bound(o);
	free_graph(g);
}

/* ********************************************************************
 * after selecting an edge and an available genome, we have three genomes
 * we can use these three genomes to initilize a graph and compute a 
 * median which will be added to the tree.
 * ********************************************************************/
void
init_graph_by_three_orders(po o, pg g, 
		int genome_from, int genome_to, int genome_avail){
	int i, j, k;
	//
	g->v_size = o->v_size;
	g->v_idx = (int**)malloc(sizeof(int*)*4);
	for(i=0; i<4; i++){
		g->v_idx[i] = (int*)malloc(sizeof(int)*o->v_size);
	}
	g->v_deg = (int**)malloc(sizeof(int*)*4);
	g->v_deg[0] = (int*)malloc(sizeof(int)*o->v_size);
	g->v_deg[1] = (int*)malloc(sizeof(int)*o->v_size);
	g->v_deg[2] = (int*)malloc(sizeof(int)*o->v_size);
	g->v_deg[3] = (int*)malloc(sizeof(int)*o->v_size);
	g->v_check = (int*)malloc(sizeof(int)*o->v_size);
	//
	g->e_size = (int*)malloc(sizeof(int)*4);
	g->e_size[0] = o->e_size[genome_from];
	g->e_size[1] = o->e_size[genome_to];
	g->e_size[2] = o->e_size[genome_avail];
	int max_e_size = max(g->e_size[0], g->e_size[1], g->e_size[2])*2;
	g->e_idx = (int**)malloc(sizeof(int*)*4);
	g->e_idx[0] = (int*)malloc(sizeof(int)*o->e_size[genome_from]);
	g->e_idx[1] = (int*)malloc(sizeof(int)*o->e_size[genome_to]);
	g->e_idx[2] = (int*)malloc(sizeof(int)*o->e_size[genome_avail]);
	g->e_idx[3] = (int*)malloc(sizeof(int)*max_e_size);
	g->e_check = (int**)malloc(sizeof(int*)*4);
	g->e_check[0] = (int*)malloc(sizeof(int)*o->e_size[genome_from]);	
	g->e_check[1] = (int*)malloc(sizeof(int)*o->e_size[genome_to]);	
	g->e_check[2] = (int*)malloc(sizeof(int)*o->e_size[genome_avail]);	
	g->e_check[3] = (int*)malloc(sizeof(int)*max_e_size);
	//
	g->distance[0] = g->distance[1] = g->distance[2] = 0;

	//v_idx
	for(i=0;i<g->v_size; i++){
		g->v_idx[0][i] = o->v_idx[genome_from][i];
		g->v_idx[1][i] = o->v_idx[genome_to][i];
		g->v_idx[2][i] = o->v_idx[genome_avail][i];
	}
	//v_check
	for(i=0;i<g->v_size; i++)
		g->v_check[i] = TRUE;
	//v_deg
	for(i=0;i<g->v_size; i++){
		int start = i==0?0:g->v_idx[0][i-1];
		int end = g->v_idx[0][i];
		int deg = end - start;
		g->v_deg[0][i] = deg;

		start = i==0?0:g->v_idx[1][i-1];
		end = g->v_idx[1][i];
		deg = end - start;
		g->v_deg[1][i] = deg;

		start = i==0?0:g->v_idx[2][i-1];
		end = g->v_idx[2][i];
		deg = end - start;
		g->v_deg[2][i] = deg;
	}
	//e_idx
	for(i=0;i<g->e_size[0]; i++)
		g->e_idx[0][i] = o->e_idx[genome_from][i];
	for(i=0;i<g->e_size[1]; i++)
		g->e_idx[1][i] = o->e_idx[genome_to][i];
	for(i=0;i<g->e_size[2]; i++)
		g->e_idx[2][i] = o->e_idx[genome_avail][i];
}

/* **************************************************************************
 * get the number of dcj operations, and insertion/deletions, and duplicaions
 * and put these information in an edge of the phylogenetic tree.
 * **************************************************************************/
void 
cal_num_edge_operations(pt tree, int edge_idx, int c1, int c2){
	int i, j, k;
	int color;
	int indel = 0;
	int dup = 0;
	int dcj = 0; //the real dcj is equal to dcj-indel
	po o = tree->orders;
	
	//convert to adj
	int ***adj = (int***)malloc(sizeof(int**)*2);	
	int **idx = (int**)malloc(sizeof(int*)*2);
	for(color =0; color<2; color++){
		int c = color==0?0:1;
		adj[c] = (int**)malloc(sizeof(int*)*o->v_size);
		idx[c] = (int*)malloc(sizeof(int)*o->v_size);
		for(i=0; i<o->v_size; i++){
			adj[c][i] = (int*)malloc(sizeof(int)*o->v_size);
			idx[c][i] = 0;
		}
	}
	for(color=0; color<2; color++){
		int c = color==0?c1:c2;
		for(i=0; i<o->v_size; i++){
			int start = i==0?0:o->v_idx[c][i-1];
			int end = o->v_idx[c][i];
			for(j=start; j<end; j++){
				adj[color][i][idx[color][i]++] = o->e_idx[c][j];
			}
		}
	}
	//compute operations
	int dis = compute_all_possible_exemplar_orders(o, c1, c2, o->v_size/2);
	int num_indel =0;
	int num_dup =0;
	int num_dcj =0;
	for(i=0; i<o->v_size; i++){
		if((idx[0][i] ==0 && idx[1][i]!=0) || (idx[1][i] ==0 && idx[0][i]!=0))
			num_indel ++;
		if((idx[0][i] !=0 && idx[1][i]!=0 && idx[0][i] != idx[1][i]) || 
			(idx[1][i] !=0 && idx[0][i]!=0 && idx[0][i] != idx[1][i]))
			num_dup ++;
	}
	//compute one indel
	int one_indel_c1=0;
	int one_indel_c2=0;
	for(i=0; i<o->v_size; i++){
		if((idx[0][i] ==0 && idx[1][i]!=0)){
			int to = adj[1][i][0];
			if(to==CAP)
				one_indel_c1 ++;
			else if((idx[0][to] ==0 && idx[1][to]!=0))
				one_indel_c1 ++;
		}
		if((idx[0][i] !=0 && idx[1][i]==0)){
			int to = adj[0][i][0];
			if(to==CAP)
				one_indel_c2 ++;
			else if(to != CAP && (idx[0][to] !=0 && idx[1][to]==0))
				one_indel_c2 ++;
		}
	}
	if(one_indel_c1!=0)
		num_indel += 1;
	if(one_indel_c2!=0)
		num_indel += 1;
	num_indel = num_indel - one_indel_c1 - one_indel_c2;
	
	num_dcj = dis - (num_indel)/2;
	tree->edges[edge_idx].from = c1;
	tree->edges[edge_idx].to = c2;
	tree->edges[edge_idx].num_indel = num_indel/2;
	tree->edges[edge_idx].num_dup = num_dup/2;
	tree->edges[edge_idx].num_dcj = num_dcj;
	tree->edges[edge_idx].score = dis;
	//free adj
	for(color=0; color<2; color++){
		for(i=0; i<o->v_size; i++){
			free(adj[color][i]);
		}
		free(idx[color]);
		free(adj[color]);
	}
	free(adj);
	free(idx);
	//graph_vis_orders_two(o, "tree.dot", c1, c2);
	//printf("distance %d \n", tree->edges[edge_idx].score);
	//printf("%d -- %d [label=\"#dcj:%d #ind:%d #dup:%d\"][color=red];\n", 
	//		tree->edges[edge_idx].from, tree->edges[edge_idx].to, tree->edges[edge_idx].num_dcj, tree->edges[edge_idx].num_indel, tree->edges[edge_idx].num_dup);
}

/* *****************************************************************************
 * copy the calculated median genome adj back to the order of the tree
 * ******************************************************************************/
void
copy_median_adj_to_order(po o, pg g){
	int i, j, k;
	int med_pos = 3;
	int med_idx = o->med_genome_idx++;
	//o->v_idx[med_idx] = (int*)malloc(sizeof(int)*o->v_size);
	//from graph to adj
	int **adj = (int**)malloc(sizeof(int*)*g->v_size);
	int *idx = (int*)malloc(sizeof(int)*g->v_size);
	for(i=0; i<g->v_size; i++){
		adj[i] = (int*)malloc(sizeof(int)*g->v_size);
		idx[i] = 0;
	}
	//csr to adj
	for(i=0; i<g->v_size; i++){
		int start = i==0?0:g->v_idx[3][i-1];
		int end = g->v_idx[3][i];
		for(j=start ; j<end; j++)
			adj[i][idx[i]++] = g->e_idx[3][j];		
	}
	//add zero 
	for(i=0; i<g->v_size; i++)
		if(g->zero[i] != NUL)
			adj[i][idx[i]++] = g->zero[i];
	//adj to csr
	//v_idx
	int sum =0;
	for(i=0; i<g->v_size; i++){
		sum += idx[i];
		g->v_idx[3][i] = sum;
	}
	//e_idx
	int pos =0;
	for(i=0; i<g->v_size; i++){
		for(j=0; j<idx[i]; j++){
			g->e_idx[3][pos++] = adj[i][j];
		}
	}
	//to orders
	o->e_size[med_idx] = g->v_idx[3][g->v_size-1];
	o->e_idx[med_idx] = (int*)malloc(sizeof(int)*o->e_size[med_idx]);
	for(i=0; i<o->v_size; i++){
		o->v_idx[med_idx][i] = g->v_idx[3][i];
	}
	for(i=0; i<o->e_size[med_idx]; i++){
		o->e_idx[med_idx][i] = g->e_idx[3][i];
	}
	//graph_vis_two(g, "before.dot", 0, 3);
	//graph_vis_orders_two(o, "after.dot", 0, med_idx);
	//free
	free(idx);
	for(i=0; i<g->v_size; i++)
		free(adj[i]);
} 

void init_best_tree(pt best, int num_genome, int num_edge){
	int i, j;
	best->num_edge = num_edge;
	best->edge_idx = 0;
	best->score = 0;
	best->upper_bound = 0;
	best->lower_bound = 0;

	best->edges = (tte*)malloc(sizeof(tte)*best->num_edge);

	int max_v_size = 10000;
	int max_e_size = 20000;
	best->orders = (po)malloc(sizeof(to));
	best->orders->num_genome = num_genome;
	best->orders->num_med_genome = (num_genome+1)/2-1;
	best->orders->num_input_genome = (num_genome+1)/2;
	best->orders->med_genome_idx = (num_genome+1)/2;
	best->orders->v_size = max_v_size;

	best->orders->e_size = (int*)malloc(sizeof(int)*best->orders->num_genome);

	best->orders->valid = (int*)malloc(sizeof(int)*best->orders->num_genome);
	best->orders->v_idx = (int**)malloc(sizeof(int*)*best->orders->num_genome);
	best->orders->e_idx = (int**)malloc(sizeof(int*)*best->orders->num_genome);

	for(i=0; i<best->orders->num_genome; i++){
		best->orders->valid[i] = FALSE;
		best->orders->e_size[i] = max_e_size;
		best->orders->v_idx[i] = (int*)malloc(sizeof(int)*best->orders->v_size);
		best->orders->e_idx[i] = (int*)malloc(sizeof(int)*best->orders->e_size[i]);
	}
}

void free_best_tree(pt best){
	int i;
	for(i=0; i<best->orders->num_genome; i++){
		free(best->orders->v_idx[i]);
		free(best->orders->e_idx[i]);
	}
	free(best->orders->e_size);
	free(best->orders->valid);
	free(best->orders->v_idx);
	free(best->orders->e_idx);
	free(best->edges);
	free(best->orders);
	free(best);
	
}

void 
vis_best_tree(pt best){
	int i, j;
	//join the rest
	po o = best->orders;
	for(i=0; i<o->num_input_genome; i++){
		if(o->valid[i] == TRUE){
			int best_pos = -1;
			int best_val = 100000;
			for(j=0; j<best->num_edge; j++){	
				if(best->edges[j].from == NUL)
					continue;
				pt neighbor_tree = (pt)malloc(sizeof(tt));
				copy_tree(best, neighbor_tree);
				add_edge(neighbor_tree, j);
				if(neighbor_tree->upper_bound < best_val){
					best_pos = j;
					best_val = neighbor_tree->upper_bound;
				}
				free_tree(neighbor_tree);
			}
			printf("the best is from %d to %d----------\n", best->edges[best_pos].from, best->edges[best_pos].to);
			pt neighbor_tree = (pt)malloc(sizeof(tt));
			copy_tree(best, neighbor_tree);
			add_edge(neighbor_tree, best_pos);
			copy_best_tree(best, neighbor_tree);
		}
	}
	FILE *writer = fopen("tree.dot", "w");
	fprintf(writer, "Graph{\n");
	for(i=0; i<best->edge_idx; i++){
		fprintf(writer, "%d -- %d [label=\"#dcj:%d #ind:%d #dup:%d\"][color=red];\n", 
				best->edges[i].from, best->edges[i].to, best->edges[i].num_dcj, best->edges[i].num_indel, best->edges[i].num_dup);
	}
	fprintf(writer, "}\n");
	fclose(writer);
}

void 
copy_best_tree(pt best, pt tree_from){
	int i, j;
	printf("***********************come here!\n");
	best->num_edge = tree_from->num_edge;
	best->edge_idx = tree_from->edge_idx;
	best->score = tree_from->score;
	best->upper_bound = tree_from->upper_bound;
	best->lower_bound = tree_from->lower_bound;

	copy_edge(best->edges, tree_from->edges, best->num_edge);


	best->orders->num_genome = tree_from->orders->num_genome;
	best->orders->num_med_genome = tree_from->orders->num_med_genome;
	best->orders->num_input_genome = tree_from->orders->num_input_genome;
	best->orders->med_genome_idx = tree_from->orders->med_genome_idx;
	best->orders->v_size = tree_from->orders->v_size;


	for(i=0; i<best->orders->num_genome; i++){
		best->orders->valid[i] = tree_from->orders->valid[i];
		best->orders->e_size[i] = tree_from->orders->e_size[i];
		if(tree_from->orders->e_size[i] != 0){
			for(j=0;j<best->orders->v_size; j++){
				best->orders->v_idx[i][j] = tree_from->orders->v_idx[i][j];
			}
			for(j=0;j<best->orders->e_size[i]; j++){
				best->orders->e_idx[i][j] = tree_from->orders->e_idx[i][j];
			}
		}
	}
}

/* ***********************************************
 * copy from the content of one tree to another
 * ***********************************************/
void 
copy_tree(pt tree_from, pt tree_to){
	int i, j;
	tree_to->num_edge = tree_from->num_edge;
	tree_to->edge_idx = tree_from->edge_idx;
	tree_to->score = tree_from->score;
	tree_to->upper_bound = tree_from->upper_bound;
	tree_to->lower_bound = tree_from->lower_bound;

	tree_to->edges = (tte*)malloc(sizeof(tte)*tree_to->num_edge);
	copy_edge(tree_to->edges, tree_from->edges, tree_to->num_edge);


	tree_to->orders = (po)malloc(sizeof(to));
	tree_to->orders->num_genome = tree_from->orders->num_genome;
	tree_to->orders->num_med_genome = tree_from->orders->num_med_genome;
	tree_to->orders->num_input_genome = tree_from->orders->num_input_genome;
	tree_to->orders->med_genome_idx = tree_from->orders->med_genome_idx;
	tree_to->orders->v_size = tree_from->orders->v_size;

	tree_to->orders->e_size = (int*)malloc(sizeof(int)*tree_to->orders->num_genome);

	tree_to->orders->valid = (int*)malloc(sizeof(int)*tree_to->orders->num_genome);
	tree_to->orders->v_idx = (int**)malloc(sizeof(int*)*tree_to->orders->num_genome);
	tree_to->orders->e_idx = (int**)malloc(sizeof(int*)*tree_to->orders->num_genome);

	for(i=0; i<tree_to->orders->num_genome; i++){
		tree_to->orders->valid[i] = tree_from->orders->valid[i];
		tree_to->orders->e_size[i] = tree_from->orders->e_size[i];
		tree_to->orders->v_idx[i] = (int*)malloc(sizeof(int)*tree_to->orders->v_size);
		if(tree_from->orders->e_size[i] != 0){
			tree_to->orders->e_idx[i] = (int*)malloc(sizeof(int)*tree_to->orders->e_size[i]);
			for(j=0;j<tree_to->orders->v_size; j++){
				tree_to->orders->v_idx[i][j] = tree_from->orders->v_idx[i][j];
			}
			for(j=0;j<tree_to->orders->e_size[i]; j++){
				tree_to->orders->e_idx[i][j] = tree_from->orders->e_idx[i][j];
			}
		}
	}
}

/* **************************
 * copy the edges in a tree
 * ***************************/
void
copy_edge(tte *edges_to, tte *edges_from, int size){
	int i;
	for(i=0; i<size; i++){
		edges_to[i].from = edges_from[i].from;
		edges_to[i].to = edges_from[i].to;
		edges_to[i].num_indel = edges_from[i].num_indel;
		edges_to[i].num_dup = edges_from[i].num_dup;
		edges_to[i].num_dcj = edges_from[i].num_dcj;
	}
}

/****************************************************************************************
 * double the search list, since the size is not large enough to keep the contents 
 * actually this step could be changed, we can use a linked list to store the search list
 * **************************************************************************************/
void 
double_tree_list(ptl l, int pos){
	int i;
	pt *tmp_l = (pt*)malloc(sizeof(pt)*l->size_list[pos]);
	for(i=0; i<l->size_list[pos]; i++){
		l->idx_list[pos]--;
		copy_tree(l->list[pos][l->idx_list[pos]], tmp_l[i]);
		free_tree(l->list[pos][l->idx_list[pos]]);
	}
	l->list[pos] = (pt*)realloc(l->list[pos], sizeof(pt)*l->size_list[pos]*2);		
	for(i=0; i<l->size_list[pos]; i++){
		copy_tree(tmp_l[i], l->list[pos][l->idx_list[pos]]);
		free_tree(tmp_l[i]);
		l->idx_list[pos]++;
	}
	l->size_list[pos] = l->size_list[pos]*2;
}

/* ******************************************
 * add one tree element to the search list
 * ******************************************/
void 
add_tree_list(ptl l, pt tree, int pos){
	int to_pos = pos-l->base;
	if(l->idx_list[to_pos]==l->size_list[to_pos])
		double_tree_list(l, to_pos);
	copy_tree(tree, l->list[to_pos][l->idx_list[to_pos]]);
	l->idx_list[to_pos]++;
	free_tree(tree);
}

/* ******************************************
 * get one tree element from the search list
 * ******************************************/
int 
pop_tree_list(ptl l, pt tree, int pos){
	int from_pos = pos-l->base;
	if(l->idx_list[from_pos]==0)
		return FALSE;
	l->idx_list[from_pos]--;
	copy_tree(l->list[from_pos][l->idx_list[from_pos]], tree);
	free_tree(l->list[from_pos][l->idx_list[from_pos]]);
	return TRUE;
}

/* ***************************************
 * free the memory allocated by a tree
 * ***************************************/
void
free_tree(pt tree){
	int i;
	free(tree->edges);
	free_orders(tree->orders);
	free(tree);
}

int 
max(int one, int two, int three){
	if(one> two && one >three)
		return one;
	else if(two > three)
		return two;
	else
		return three;
}

void graph_vis_orders_two(po o, char *file, int c1, int c2){
	int i=0,j,from;
	FILE *writer = fopen(file, "w");
	fprintf(writer, "graph G{\n");
	for(from=0;from<o->v_size; from++){
		//if(o->v_check[from] != NUL){
		int start = from==0?0:o->v_idx[c1][from-1];
		int end = o->v_idx[c1][from];
		for(j=start;j<end;j++){
			int to = o->e_idx[c1][j];
			fprintf(writer, "%d -- %d [color=red];\n", from, to);
		}   
		//}
	}
	for(from=0;from<o->v_size; from++){
		//if(o->v_check[from] != NUL){
		int start = from==0?0:o->v_idx[c2][from-1];
		int end = o->v_idx[c2][from];
		for(j=start;j<end;j++){
			int to = o->e_idx[c2][j];
			fprintf(writer, "%d -- %d [color=blue];\n", from, to);
		}   
		//}
	}
	fprintf(writer, "}\n");
	fclose(writer);
}

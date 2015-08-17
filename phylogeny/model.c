#include "model.h"
#include "tree.h"
#include "newick.h"

void
graph_to_orders(pg g, pt tree);
void
init_graph_by_two_orders(po o, pg g, 
		int genome_from, int genome_to);
void
init_directive_node(pg g, pt tree,
		pnt nw, int root);
void 
adj_to_graph_csr(pg g, pa adj, int pos);
void 
init_adj_directive_node(po o, pa a, int o_pos);
int 
search(int *dic, int size, int word);

/*
 * this is to implement the GASTS idea for initilize genomes for internal nodes given a model tree
 **/
void
run_initialization(char *file, char *nwfile){
	int i, j, k;
	//initialize orders
	po o = (po)malloc(sizeof(to));	
	init_orders(file, o);
	pt tree = (pt)malloc(sizeof(tt));
	init_null_tree(o, tree); //get a tree that only genomes every internal node is not initialized, there is no edges 
	
	//initialize the pnt tree for reading tree structures 
	//from newick format to edges then from edges to pnt tree
	int *map = (int*)malloc(sizeof(int)*tree->orders->num_input_genome);
	for(i=0; i<tree->orders->num_input_genome; i++)
		map[i] = i;
	from_newick_to_edge(file, tree);
	pnt n_tree = (pnt)malloc(sizeof(tnt));
	from_edges_to_n_tree(n_tree, tree->edges, map, 
			tree->orders->num_input_genome, 
			tree->num_edge);

	//while node is not initialized, select a node that has two
	//neighbors that has been initialized, and one not initialized
	int has_one_with_two_init = TRUE;
	for(i=tree->orders->num_input_genome; 
		i<tree->orders->num_genome; i++){
		pg g = (pg)malloc(sizeof(tg));
		//two neighbors are initiated
		if((n_tree->e_type[i][0]+n_tree->e_type[i][1]+n_tree->e_type[i][2])==2){
			int to[2];
			for(k=0; k<2; k++){
				int c = k;
				for(j=c==0?0:to[1]+1; j<3; j++){
					if(n_tree->e_type[i][j] == E_TYPE_UNINIT &&
							n_tree->v_type[n_tree->edges[i][j]] == V_TYPE_INIT){
						to[c] = n_tree->edges[i][j];
						n_tree->e_type[i][j] = E_TYPE_INIT;
						break;
					}
				}
			}
			int root = -1;
			for(j=0; j<3; j++)
				if(n_tree->v_type[n_tree->edges[i][j]] == V_TYPE_UNINIT)
					root = n_tree->edges[i][j];
			init_graph_by_two_orders(tree->orders, g, 
					to[0], to[1]);
			init_directive_node(g, tree, root);
			graph_to_orders(g, tree);
		}
		//three neighbors are initiated
		else if((n_tree->e_type[i][0]+n_tree->e_type[i][1]+n_tree->e_type[i][2])==3){
			init_graph_by_three_orders(tree->orders, g, 
				n_tree->edges[i][0], n_tree->edges[i][1], 
				n_tree->edges[i][2]);
			graph_to_orders(g, tree);
		}
		n_tree->v_type[i] = V_TYPE_UNINIT;
	}
}

/*
 * After initialization, the iterative refinement is performed,
 * this method is used in GRAPPA, BPAnalysis and MGR. 
 **/
void 
run_refinement(pt tree){
	int i, j, k;
	int *map = (int*)malloc(sizeof(int)*tree->orders->num_input_genome);
	for(i=0; i<tree->orders->num_input_genome; i++)
		map[i] = i;
	pnt n_tree = (pnt)malloc(sizeof(tnt));
	from_edges_to_n_tree(n_tree, tree->edges, map, 
			tree->orders->num_input_genome, 
			tree->num_edge);
	po o = tree->orders;
	int improved = TRUE;
	int *attachment_score = (int*)malloc(sizeof(int)*tree->orders->num_genome);
	for(i=0; i<tree->orders->num_genome; i++){
		attachment_score[i] = 0;
		for(j=0; j<tree->num_edge; j++){
			int from = tree->edges[j].from;
			int to = tree->edges[j].to;
			if(from == i || to == i)
				attachment_score[i] += tree->edges[j].score;
		}
	}
	while(improved == TRUE){
		for(i=o->num_input_genome; i<o->num_genome; i++){
			int old_score = 
			pg g = (pg)malloc(sizeof(tg));
			init_graph_by_three_orders(tree->orders, g, 
					n_tree->edges[i][0], n_tree->edges[i][1], 
					n_tree->edges[i][2]);
			if(g->score < attachment_score[i]){ //update genomes
				//update gene content
				tree->orders->e_size[i] = g->e_size[3];
				free(tree->orders->e_idx[i]);
				tree->orders->e_idx[i] = (int*)malloc(sizeof(int)*g->e_size[3]);
				for(j=0; j<g->v_size; j++)
					tree->orders->v_idx[i][j] = g->v_idx[3][j];
				for(j=0; j<g->e_size[3]; j++)
					tree->orders->e_idx[i][j] = g->e_idx[3][j];
				//update edges
				for(j=0; j<tree->edge_num; j++){
					int from = tree->edges[j].from;
					int to = tree->edges[j].to;
					for(k=0; k<3; k++){
						if((to == n_tree->edges[i][k] && from == i) ||
							(from == n_tree->edges[i][k] && to == i))
							tree->edges[j].score = g->distance[k];
					}
				}
				//update attachment_score
				attachment_score[i] = g->score;
				//update improved
				improved = TRUE;
			}
		}
	}
	free(attachment_score);
}

/*
 * turn the median order of the graph back into the med order.
 * because we iterate the order right the right direction, there
 * is no need to input the position info.
 **/
void
graph_to_orders(pg g, pt tree){
	po oo = tree->orders;
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

	//get the upper and lower_bound
	tree->score = tree->edges[edge_from_med].score + tree->edges[edge_to_med].score + tree->edges[edge_avail_med].score;
	free_graph(g);
}

/*
 * if only two neighbors are initialized, we need to init graph by 
 * using the information from only these two orders.
 **/
void
init_graph_by_two_orders(po o, pg g, 
		int genome_from, int genome_to){
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
	g->e_idx = (int**)malloc(sizeof(int*)*4);
	g->e_idx[0] = (int*)malloc(sizeof(int)*o->e_size[genome_from]);
	g->e_idx[1] = (int*)malloc(sizeof(int)*o->e_size[genome_to]);
	g->e_check = (int**)malloc(sizeof(int*)*4);
	g->e_check[0] = (int*)malloc(sizeof(int)*o->e_size[genome_from]);	
	g->e_check[1] = (int*)malloc(sizeof(int)*o->e_size[genome_to]);	
	//
	g->distance[0] = g->distance[1] = g->distance[2] = 0;

	//v_idx
	for(i=0;i<g->v_size; i++){
		g->v_idx[0][i] = o->v_idx[genome_from][i];
		g->v_idx[1][i] = o->v_idx[genome_to][i];
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
	}
	//e_idx
	for(i=0;i<g->e_size[0]; i++)
		g->e_idx[0][i] = o->e_idx[genome_from][i];
	for(i=0;i<g->e_size[1]; i++)
		g->e_idx[1][i] = o->e_idx[genome_to][i];
}

/*
 * this is actually the hardest part,
 * init orders by it's directive nodes.
 * performing a bfs to get all degrees and try to optimize it.
 **/
void
init_directive_node(pg g, pt tree,
		pnt nw, int root){
	int i, j, k;
	int *queue = (int*)malloc(sizeof(int)*tree->num_genome);
	int q_idx = 0;
	queue[q_idx++] = root;
	po o = tree->orders;
	//bfs to get adj graph
	pa adj = (pa)malloc(sizeof(ta));
	init_adj(adj, g->v_size, o->num_genome);
	while(q_idx>0){
		//pop queue
		int current = queue[--q_idx];
		if(nw->v_type[current] == V_TYPE_INIT)
			init_adj_directive_node(o, a, current);
		else{ // add into queue
			for(i=0; i<3; i++)
				queue[q_idx++] = nw->edges[current][i];
		}
	}
	//bfs to init graph
	adj_to_graph_csr(g, adj, 2);
}

void 
adj_to_graph_csr(pg g, pa adj, int pos){
	int i, j, k;
	//compute v_idx first
	int sum = 0;
	for(i=0; i<g->v_size; i++){
		if(a->idx_list[i] == 0){
			g->v_idx[pos][i] = sum;
			continue;
		}
		for(j=0; j<a->idx_list[i]; j++){
			if(a->list[i][j] != NUL)
				sum++;
		}
		g->v_idx[pos][i] = sum;
	}
	//compute e_idx
	k=0;
	for(i=0;i<g->v_size; i++){
		if(a->idx_list[i] == 0)
			continue;
		for(j=0; j<a->idx_list[i]; j++){
			if(a->list[i][j] != NUL){
				g->e_idx[pos][k] = a->list[i][j];
				k++;
			}
		}
	}
}

void 
init_adj_directive_node(po o, pa a, int o_pos){
	int i, j, k;
	for(i=0; i<o->v_size; i++){
		int start = i==0?0:o->v_idx[pos][i-1];
		int end = o->v_idx[pos][i];
		for(j=start; j<end; j++){
			int to = o->e_idx[pos][j];
			int from = i;
			if(search(a->list[from], idx_list[from], to) == FALSE)
				a->list[from][idx_list[from]++] = to;
		}
	}
}

int 
search(int *dic, int size, int word){
	int i;
	for(i=0; i<size; i++)
		if(dic[i] == word)
			return TRUE;
	return FALSE;
}
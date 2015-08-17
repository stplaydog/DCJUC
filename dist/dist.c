#include "graph.h"
#include "dist.h"

int gene_num[2];
int chr_num[2];
//int vet_type[MAX_VET_SIZE];

void set_visited(int *vet_type, int pos);
void vis_adj(char *file, int **adj, int size);
void 
compress_graph(int **comp_v_idx, int **comp_e_idx, 
		int *comp_v_check, int **v_deg,
		int size, int *new_v_size, 
		int *new_e_size1, int *new_e_size2);
void 
vis_compact_adj(char *file, int **comp_v_idx,
		int **comp_e_idx, int v_size);

int IS_TWO=FALSE;

extern int turn_on_debug;
extern int HAS_PROBLEM;
extern int debug_nj;

//these are global variables for computing distances
//int **comp_v_idx;
//int **comp_e_idx;
//int **comp_v_deg;
//int *comp_v_check;
//int **period_arr;
//int *violate; 
//int ***adj;
//int *vet_type;
//int *tmp_idx;
//int **tmp_v_idx; 
//int **tmp_e_idx; 
//int **tmp_v_deg;

//struct search_list{
//	int list_num;
//	p_elem *list;
//	int *index; //record the index of each element
//};
//typedef struct search_list t_l;
//typedef struct search_list *p_l;




//void init_search_list()
//{
//	p_list search_list = (p_list)malloc(sizeof(t_list));
//	search_list->list_num = 1000;
//	search_list->list = (p_elem*)malloc(sizeof(p_elem)*search_list->list_num);	
//	search_list->index = (int*)malloc(sizeof(int)*search_list->list_num);
//	
//}
//

//use recursive algorithm to solve this problem
void compute_indel_exemplar_dis(pg g, int c1, int c2, int *dis, int gene, int pos){
	int i,j, k, w;
}

void init_for_dis(){
	int i;
	//int max_v_size = 1000;
	//int max_e_size = 2000;
	//int max_period = 100;
	//int possible_comb = 100000;
	//comp_v_idx = (int**)malloc(sizeof(int*)*2);
	//comp_v_idx[0] = (int*)malloc(sizeof(int)*max_v_size);
	//comp_v_idx[1] = (int*)malloc(sizeof(int)*max_v_size);
	//comp_e_idx = (int**)malloc(sizeof(int*)*2);
	//comp_e_idx[0] = (int*)malloc(sizeof(int)*max_e_size);
	//comp_e_idx[1] = (int*)malloc(sizeof(int)*max_e_size);
	//comp_v_deg = (int**)malloc(sizeof(int*)*2);
	//comp_v_deg[0] = (int*)malloc(sizeof(int)*max_v_size);
	//comp_v_deg[1] = (int*)malloc(sizeof(int)*max_v_size);
	//comp_v_check = (int*)malloc(sizeof(int)*max_v_size);
	//period_arr = (int**)malloc(sizeof(int*)*2);
	//period_arr[0] = (int*)malloc(sizeof(int)*max_period);
	//period_arr[1] = (int*)malloc(sizeof(int)*max_period);
	//violate = (int*)malloc(sizeof(int)*possible_comb);
	//adj = (int***)malloc(sizeof(int**)*possible_comb);
	//for(i=0;i<possible_comb;i++){
	//	adj[i] = (int**)malloc(sizeof(int*)*2);
	//	adj[i][0] = (int*)malloc(sizeof(int)*max_v_size);
	//	adj[i][1] = (int*)malloc(sizeof(int)*max_v_size);
	//}
	//vet_type = (int*)malloc(sizeof(int)*max_v_size);
	//tmp_idx = (int*)malloc(sizeof(int)*max_v_size);
	//rev_tmp_idx = (int*)malloc(sizeof(int)*max_v_size);
	//tmp_v_idx = (int**)malloc(sizeof(int*)*2); 
	//tmp_e_idx = (int**)malloc(sizeof(int*)*2); 
	//tmp_v_deg = (int**)malloc(sizeof(int*)*2);
	//tmp_v_idx[0] = (int*)malloc(sizeof(int)*max_v_size); 
	//tmp_v_idx[1] = (int*)malloc(sizeof(int)*max_v_size); 
	//tmp_v_deg[0] = (int*)malloc(sizeof(int)*max_v_size); 
	//tmp_v_deg[1] = (int*)malloc(sizeof(int)*max_v_size); 
	//tmp_e_idx[0] = (int*)malloc(sizeof(int)*max_e_size); 
	//tmp_e_idx[1] = (int*)malloc(sizeof(int)*max_e_size); 
}

int 
cpec_mem_save(pg g1, int c1, int c2, int gene){
	int i, j, k;
	int color;
#ifdef USE_EXPERIMENT
	int new_v_size = 500;
	int num_dup = 20;
	int num_dcj = 100;
	int ***comp_adj = (int***)malloc(sizeof(int**)*2);
	int **comp_idx = (int**)malloc(sizeof(int*)*2);
	for(color=0; color<2; color++){
		comp_adj[color] = (int**)malloc(sizeof(int*)*new_v_size);
		comp_idx[color] = (int*)malloc(sizeof(int)*new_v_size);
		for(j=0; j<new_v_size; j++){
			comp_adj[color][j] = (int*)malloc(sizeof(int)*new_v_size);
			comp_idx[color][j] = 0;
		}
	}
	//init adj
	comp_adj[0][0][comp_idx[0][0]++] = CAP;
	comp_adj[1][0][comp_idx[1][0]++] = CAP;
	comp_adj[0][new_v_size-1][comp_idx[0][new_v_size-1]++] = CAP;
	comp_adj[1][new_v_size-1][comp_idx[1][new_v_size-1]++] = CAP;
	for(color=0; color<2; color++){
		for(i=1; i<new_v_size-1; i+=2){
			comp_adj[color][i][comp_idx[color][i]++] = i+1;
			comp_adj[color][i+1][comp_idx[color][i+1]++] = i;
		}
	}
	//perform DCJ
	srand (time(NULL));
	for(i=0; i<num_dcj; i++){
		int conn_type = rand() % 2;
		int conn_color = rand() % 2;
		int one = rand() % new_v_size;
		int one_to = comp_adj[conn_color][one][0];
		while(one_to == CAP){
			one = rand() % new_v_size;
			one_to = comp_adj[conn_color][one][0];
		}
		int two = rand() % new_v_size;
		int two_to = comp_adj[conn_color][two][0];
		while(two == one || two == one_to || two_to == CAP){
			two = rand() % new_v_size;
			two_to = comp_adj[conn_color][two][0];
		}
		if(conn_type == 0){
			comp_adj[conn_color][one][0] = two;
			comp_adj[conn_color][two][0] = one;
			comp_adj[conn_color][one_to][0] = two_to;
			comp_adj[conn_color][two_to][0] = one_to;
		}
		else if(conn_type == 1){
			comp_adj[conn_color][one][0] = two_to;
			comp_adj[conn_color][two_to][0] = one;
			comp_adj[conn_color][one_to][0] = two;
			comp_adj[conn_color][two][0] = one_to;
		}
	}
	//perform duplication
	for(i=0; i<num_dup; i++){
		int color = rand() % 2;
		int vet = rand() % new_v_size;
		int vet_to;
		while(1){
			vet_to = rand() % new_v_size;
			if(vet_to == vet)
				continue;
			int to_same = FALSE;
			for(j=0; j<comp_idx[color][vet]; j++)
				if(comp_adj[color][vet][j] == vet_to){
					to_same = TRUE;
					break;
				}
			if(to_same == TRUE)
				continue;
			break;
		}
		comp_adj[color][vet][comp_idx[color][vet]++] = vet_to;
		comp_adj[color][vet_to][comp_idx[color][vet_to]++] = vet;
	}
	//get the e_size
	int new_e_size[2];
	for(color=0; color<2; color++){
		new_e_size[color] = 0;
		for(i=0; i< new_v_size; i++)
			new_e_size[color] += comp_idx[color][i];
	}
	//from adj to csr
	int **comp_v_idx = (int**)malloc(sizeof(int*)*2);
	comp_v_idx[0] = (int*)malloc(sizeof(int)*new_v_size);
	comp_v_idx[1] = (int*)malloc(sizeof(int)*new_v_size);
	int **comp_e_idx = (int**)malloc(sizeof(int*)*2);
	comp_e_idx[0] = (int*)malloc(sizeof(int)*new_e_size[0]);
	comp_e_idx[1] = (int*)malloc(sizeof(int)*new_e_size[1]);
	//v_idx
	for(color =0; color<2; color++){
		int sum =0;
		for(i=0; i<new_v_size; i++){
			sum += comp_idx[color][i];
			comp_v_idx[color][i] = sum;
		}
	}
	//e_idx
	int pos =0;
	for(color=0; color<2; color++){
		pos = 0;
		for(i=0; i<new_v_size; i++)
			for(j=0; j<comp_idx[color][i]; j++)
				comp_e_idx[color][pos++] = comp_adj[color][i][j];
	}
	//free
	for(i=0; i<2; i++){
		for(j=0; j<new_v_size; j++){
			free(comp_adj[i][j]);
		}
		free(comp_adj[i]);
		free(comp_idx[i]);
	}
	free(comp_adj);
	free(comp_idx);
#else
	pg g = (pg)malloc(sizeof(tg));
	if(IS_TWO == FALSE){
		prepare_graph_for_distance(g1, g);
	}
	else{
		prepare_graph_for_distance_two(g1, g);
	}
	int **comp_v_idx = (int**)malloc(sizeof(int*)*2);
	comp_v_idx[0] = (int*)malloc(sizeof(int)*g->v_size);
	comp_v_idx[1] = (int*)malloc(sizeof(int)*g->v_size);
	int **comp_e_idx = (int**)malloc(sizeof(int*)*2);
	comp_e_idx[0] = (int*)malloc(sizeof(int)*g->e_size[c1]);
	comp_e_idx[1] = (int*)malloc(sizeof(int)*g->e_size[c2]);
	int **comp_v_deg = (int**)malloc(sizeof(int*)*2);	
	comp_v_deg[0] = (int*)malloc(sizeof(int)*g->v_size);
	comp_v_deg[1] = (int*)malloc(sizeof(int)*g->v_size);
	for(i=0;i<g->v_size; i++){
		comp_v_idx[0][i] = g->v_idx[c1][i];
		comp_v_idx[1][i] = g->v_idx[c2][i];
	}
	for(i=0; i<g->e_size[c1]; i++)
		comp_e_idx[0][i] = g->e_idx[c1][i];
	for(i=0; i<g->e_size[c2]; i++)
		comp_e_idx[1][i] = g->e_idx[c2][i];
	int *comp_v_check = (int*)malloc(sizeof(int)*g->v_size);
	for(i=0; i<g->v_size; i++){
		comp_v_check[i] = g->v_check[i]==NUL?NUL:TRUE;
	}
	int compressed_size;
	int new_v_size;
	int new_e_size1;
	int new_e_size2;
	compress_graph(comp_v_idx, comp_e_idx, comp_v_check, comp_v_deg, g->v_size, &new_v_size, &new_e_size1, &new_e_size2);
	if(IS_TWO==FALSE)
		free_graph(g);
	else
		free_graph_two(g);
	//FILE *dc = fopen("vis_com.dot", "w");
	//fprintf(dc, "Graph{\n");
	//for(i=0;i<new_v_size; i++){
	//	int start = i==0?0:comp_v_idx[0][i-1];
	//	int end = comp_v_idx[0][i];
	//	for(j=start; j<end; j++)
	//		fprintf(dc, "%d -- %d [color=red];\n", i, comp_e_idx[0][j]);
	//	start = i==0?0:comp_v_idx[1][i-1];
	//	end = comp_v_idx[1][i];
	//	for(j=start; j<end; j++)
	//		fprintf(dc, "%d -- %d [color=blue];\n", i, comp_e_idx[1][j]);
	//}
	//fprintf(dc, "}\n");
	//fclose(dc);
#endif

	//vis_compact_adj("exp1.dot", comp_v_idx, comp_e_idx, new_v_size);
	//this is the part to //compute the exemplar distance
	int ***adj_dir = (int***)malloc(sizeof(int**)*2);
	int **idx_dir = (int**)malloc(sizeof(int*)*2);
	for(color =0; color<2; color++){
		adj_dir[color] = (int**)malloc(sizeof(int*)*new_v_size);
		idx_dir[color] = (int*)malloc(sizeof(int)*new_v_size);
		for(i=0; i<new_v_size; i++){
			adj_dir[color][i] = (int*)malloc(sizeof(int)*new_v_size);
			idx_dir[color][i] = 0;
		}
	}
	//init adj
	for(color=0; color<2; color++){
		for(i=0; i<new_v_size; i++){
			int start = i==0?0:comp_v_idx[color][i-1];
			int end = comp_v_idx[color][i];
			for(j=start; j<end; j++){
				int from = i;
				int to = comp_e_idx[color][j];
				if(to>=new_v_size || to<0)
					continue;
				//check if exists we only keep directed edges
				int exist=FALSE;
				if(to != CAP)
					for(k=0; k<idx_dir[color][to]; k++)
						if(adj_dir[color][to][k]==from)
							exist = TRUE;
				if(exist == FALSE)
					adj_dir[color][from][idx_dir[color][from]++] = to;
			}
		}
	}
	//FILE *writer = fopen("exp1.dot", "w");
	//fprintf(writer, "Graph{\n");
	//for(color=0; color<2; color++){
	//	for(i=0; i<new_v_size; i++){
	//		for(j=0; j<idx_dir[color][i]; j++)
	//			fprintf(writer, "%d -- %d [color=%s];\n", i, adj_dir[color][i][j], color==0?"red":"blue");
	//	}
	//}
	//fprintf(writer, "}\n");
	//fclose(writer);
	//compute possible comb
	int possible_comb = 1;
	for(i=0; i<new_v_size; i++){
		if(idx_dir[0][i]!=0 && idx_dir[1][i]==0)
			possible_comb = possible_comb*idx_dir[0][i];
		else if(idx_dir[0][i]==0 && idx_dir[1][i]!=0)
			possible_comb = possible_comb*idx_dir[1][i];
		else if(idx_dir[0][i]!=0 && idx_dir[1][i]!=0)
			possible_comb = possible_comb*idx_dir[0][i]*idx_dir[1][i];
	}
	//iterating
	int iter_num;
	double samp_rate = 0.2;
	int *is_sampled = (int*)malloc(sizeof(int)*possible_comb);
	for(i=0; i<possible_comb*samp_rate; i++)
		is_sampled[i] = FALSE;
	int **dis_adj = (int**)malloc(sizeof(int*)*2);
	dis_adj[0] = (int*)malloc(sizeof(int)*new_v_size);
	dis_adj[1] = (int*)malloc(sizeof(int)*new_v_size);
	int *dis_type = (int*)malloc(sizeof(int)*new_v_size);
	int **selected = (int**)malloc(sizeof(int*)*2);
	selected[0] = (int*)malloc(sizeof(int)*new_v_size);
	selected[1] = (int*)malloc(sizeof(int)*new_v_size);
	int max =0;
	/////////////////////////////
	if(possible_comb > 100)
		iter_num = (double)possible_comb * samp_rate;
	else
		iter_num = (double)possible_comb;
	for(i=0; i<iter_num; i++){
		int sample_pos;
		if(iter_num != possible_comb){
			sample_pos = rand() % possible_comb;
			while(is_sampled[sample_pos] == TRUE)
				sample_pos = rand() % possible_comb;
		}
		else
			sample_pos = i;
		int possible_remain = sample_pos;
		for(color=0; color<2; color++){
			for(j=0; j<new_v_size; j++){
				selected[color][j] = FALSE;
				dis_adj[color][j] = NUL;
			}
		}
		//this is to choose the adjacency
		for(j=0; j<new_v_size; j++){
			int possible_one = idx_dir[0][j]==0?1:idx_dir[0][j];
			int possible_two = idx_dir[1][j]==0?1:idx_dir[1][j];
			int pos_both = possible_remain % (possible_one * possible_two);
			for(color=0; color<2; color++){
				int chose = pos_both % (color==0?possible_one:possible_two);
				int find = FALSE;
				if(idx_dir[color][j]!=0 && selected[color][j] != TRUE){
					for(k=0; k<idx_dir[color][j]; k++){
						if(adj_dir[color][j][chose]!= CAP && 
								selected[color][adj_dir[color][j][chose]] == TRUE){
							chose = (pos_both+k) % (color==0?possible_one:possible_two);
						}
						else{
							find = TRUE;
							break;
						}
					}
				}
				if(find == TRUE){
					if(adj_dir[color][j][chose] != CAP){
						dis_adj[color][j] = adj_dir[color][j][chose];
						dis_adj[color][adj_dir[color][j][chose]] = j;
						selected[color][adj_dir[color][j][chose]] = TRUE;
					}
					else{
						dis_adj[color][j] = adj_dir[color][j][chose];
					}
					selected[color][j] = TRUE;
				}
			}
			possible_remain = possible_remain / (possible_one * possible_two);
		}
		for(j=0;j<new_v_size; j++)
		{
			if(dis_adj[1][j]!=NUL && dis_adj[1][j]!=CAP && dis_adj[0][j]==NUL)
				dis_type[j]=PI_OPEN;
			else if(dis_adj[0][j]!=NUL && dis_adj[0][j]!=CAP && dis_adj[1][j]==NUL)
				dis_type[j]=GAMMA_OPEN;
			else if(dis_adj[0][j]!=NUL && dis_adj[1][j]!=NUL &&
					dis_adj[0][j]!=CAP && dis_adj[1][j]!=CAP)
				dis_type[j]=CYCLE;
			else if((dis_adj[0][j]!=NUL && dis_adj[0][j]==CAP) ||
					(dis_adj[1][j]!=NUL && dis_adj[1][j]==CAP))
				dis_type[j]=END;
			else
				dis_type[j]=NONE;
		}
		int results = compute_indel_dis(dis_adj, dis_type, 0, 1, new_v_size, new_v_size);
		//printf("result %d %d\n", results, sample_pos);
		if(results > max)
			max = results;
	}
	//free comp_***
	free(comp_v_idx[0]);
	free(comp_v_idx[1]);
	free(comp_v_idx);
	free(comp_v_deg[0]);
	free(comp_v_deg[1]);
	free(comp_v_deg);
	free(comp_e_idx[0]);
	free(comp_e_idx[1]);
	free(comp_e_idx);
	free(selected[0]);
	free(selected[1]);
	free(selected);
	free(dis_type);
	free(dis_adj[0]);
	free(dis_adj[1]);
	free(dis_adj);
	free(is_sampled);
	free(comp_v_check);
	//free adj
	for(color=0; color<2; color++){
		free(idx_dir[color]);
		for(i=0; i<new_v_size; i++){
			free(adj_dir[color][i]);
		}
		free(adj_dir[color]);
	}
	free(idx_dir);
	free(adj_dir);
	return max;
}

int 
compute_all_possible_exemplar_compact(pg g1, int c1, int c2, int gene){
	pg g = (pg)malloc(sizeof(tg));
	if(IS_TWO == FALSE){
		prepare_graph_for_distance(g1, g);
	}
	else{
		prepare_graph_for_distance_two(g1, g);
	}
	int i, j, k;
	int **comp_v_idx = (int**)malloc(sizeof(int*)*2);
	comp_v_idx[0] = (int*)malloc(sizeof(int)*g->v_size);
	comp_v_idx[1] = (int*)malloc(sizeof(int)*g->v_size);
	int **comp_e_idx = (int**)malloc(sizeof(int*)*2);
	comp_e_idx[0] = (int*)malloc(sizeof(int)*g->e_size[c1]);
	comp_e_idx[1] = (int*)malloc(sizeof(int)*g->e_size[c2]);
	int **comp_v_deg = (int**)malloc(sizeof(int*)*2);	
	comp_v_deg[0] = (int*)malloc(sizeof(int)*g->v_size);
	comp_v_deg[1] = (int*)malloc(sizeof(int)*g->v_size);
	for(i=0;i<g->v_size; i++){
		comp_v_idx[0][i] = g->v_idx[c1][i];
		comp_v_idx[1][i] = g->v_idx[c2][i];
	}
	for(i=0; i<g->e_size[c1]; i++)
		comp_e_idx[0][i] = g->e_idx[c1][i];
	for(i=0; i<g->e_size[c2]; i++)
		comp_e_idx[1][i] = g->e_idx[c2][i];
	int *comp_v_check = (int*)malloc(sizeof(int)*g->v_size);
	for(i=0; i<g->v_size; i++){
		comp_v_check[i] = g->v_check[i]==NUL?NUL:TRUE;
	}
	int compressed_size;
	int new_v_size;
	int new_e_size1;
	int new_e_size2;
	//graph_vis_two(g, "vis.dot", c1, c2);
	compress_graph(comp_v_idx, comp_e_idx, comp_v_check, comp_v_deg, g->v_size, &new_v_size, &new_e_size1, &new_e_size2);
	FILE *dc = fopen("vis_com.dot", "w");
	fprintf(dc, "Graph{\n");
	for(i=0;i<new_v_size; i++){
		int start = i==0?0:comp_v_idx[0][i-1];
		int end = comp_v_idx[0][i];
		for(j=start; j<end; j++)
			fprintf(dc, "%d -- %d [color=red];\n", i, comp_e_idx[0][j]);
		start = i==0?0:comp_v_idx[1][i-1];
		end = comp_v_idx[1][i];
		for(j=start; j<end; j++)
			fprintf(dc, "%d -- %d [color=blue];\n", i, comp_e_idx[1][j]);
	}
	fprintf(dc, "}\n");
	fclose(dc);

	//this is the part to //compute the exemplar distance
	int possible_comb=1;
	int color;
	for(i=0; i<new_v_size; i++){
		if(comp_v_deg[0][i] != 0 && comp_v_deg[1][i]!=0){
			possible_comb = possible_comb*(comp_v_deg[0][i]*comp_v_deg[1][i]);
		}
	}
	int ***adj = (int***)malloc(sizeof(int**)*possible_comb);
	for(i=0;i<possible_comb;i++){
		adj[i] = (int**)malloc(sizeof(int*)*2);
		adj[i][0] = (int*)malloc(sizeof(int)*new_v_size);
		adj[i][1] = (int*)malloc(sizeof(int)*new_v_size);
	}

	int period =0;
	for(i=0;i<new_v_size;i++){
		period = 0;
		if(comp_v_deg[0][i]!=0 && comp_v_deg[1][i]!=0){
			period = comp_v_deg[0][i]*comp_v_deg[1][i];
			int **period_arr = (int**)malloc(sizeof(int*)*2);
			period_arr[0] = (int*)malloc(sizeof(int)*period);
			period_arr[1] = (int*)malloc(sizeof(int)*period);
			int start1 = i==0?0:comp_v_idx[0][i-1];
			int start2 = i==0?0:comp_v_idx[1][i-1];
			int end1 = comp_v_idx[0][i];
			int end2 = comp_v_idx[1][i];
			int pos=0;
			for(j=start1;j<end1;j++){
				for(k=start2;k<end2;k++){	
					period_arr[0][pos] = comp_e_idx[0][j];
					period_arr[1][pos] = comp_e_idx[1][k];
					pos++;
				}
			}
			for(j=0;j<possible_comb;j++){
				adj[j][0][i] = period_arr[0][j%period];
				adj[j][1][i] = period_arr[1][j%period];
			}
			free(period_arr[0]);
			free(period_arr[1]);
			free(period_arr);
		}
		else{
			int **period_arr = (int**)malloc(sizeof(int*)*2);
			if(comp_v_deg[0][i]==0 && comp_v_deg[1][i]==0)
				;
			else if(comp_v_deg[0][i]==0 && comp_v_deg[1][i]!=0){
				period = comp_v_deg[1][i];
				period_arr[0] = (int*)malloc(sizeof(int)*period);
				period_arr[1] = (int*)malloc(sizeof(int)*period);
				int start2 = i==0?0:comp_v_idx[1][i-1];
				int end2 = comp_v_idx[1][i];
				int pos=0;
				for(k=start2;k<end2;k++){
					period_arr[0][pos] = NUL;
					period_arr[1][pos] = comp_e_idx[1][k];
					pos++;
				}
			}
			else if(comp_v_deg[0][i]!=0 && comp_v_deg[1][i]==0){
				period = comp_v_deg[0][i];
				period_arr[0] = (int*)malloc(sizeof(int)*period);
				period_arr[1] = (int*)malloc(sizeof(int)*period);
				int start1 = i==0?0:comp_v_idx[0][i-1];
				int end1 = comp_v_idx[0][i];
				int pos=0;
				for(k=start1;k<end1;k++){
					period_arr[0][pos] = comp_e_idx[0][k];
					period_arr[1][pos] = NUL;
					pos++;
				}
			}
			if(period != 0){
				for(j=0;j<possible_comb;j++){
					adj[j][0][i] = period_arr[0][j%period];
					adj[j][1][i] = period_arr[1][j%period];
				}
			}
			else{
				for(j=0;j<possible_comb;j++){
					adj[j][0][i] = NUL;
					adj[j][1][i] = NUL;
				}
			}
			free(period_arr[0]);
			free(period_arr[1]);
			free(period_arr);
		}
	}
	//
	for(i=0;i<possible_comb;i++){
		for(j=0;j<new_v_size;j++){
			int to = adj[i][0][j];
			if(to!=NUL && to != CAP){
				if(adj[i][0][to] != j)
					adj[i][0][j]=NUL;
			}
			to = adj[i][1][j];
			//if(to > 100000){
			//	vis_adj("adj.dot", adj[i], new_v_size);
			//	FILE *w = fopen("comp.dot", "w");
			//	int m, n;
			//	fprintf(w, "Graph{\n");
			//	for(m=0; m<new_v_size; m++){
			//		int start = m==0?0:comp_v_idx[0][m-1];
			//		int end = comp_v_idx[0][m];
			//		for(n=start; n<end; n++)
			//			fprintf(w, "%d -- %d [color=%s];\n", i, comp_e_idx[0][n], "red");
			//		start = m==0?0:comp_v_idx[1][m-1];
			//		end = comp_v_idx[1][m];
			//		for(n=start; n<end; n++)
			//			fprintf(w, "%d -- %d [color=%s];\n", i, comp_e_idx[1][n], "blue");
			//	}
			//	fprintf(w, "}\n");
			//	fclose(w);
			//	graph_vis_two(g1, "g1.dot", c1, c2);
			//	graph_vis_two(g, "g.dot", c1, c2);
			//	HAS_PROBLEM = TRUE;
			//	return -1;
			//}
			if(to!=NUL && to != CAP){
				if(adj[i][1][to] != j)
					adj[i][1][j]=NUL;
			}
		}
	}
	//int c;
	//for(c=0; c<possible_comb; c++){	
	//	char f_buf[100];
	//	sprintf(f_buf, "vis_com%d.dot", c);
	//	FILE *writer = fopen(f_buf, "w");
	//	fprintf(writer, "Graph{\n");
	//	for(i=0;i<new_v_size; i++){
	//		int from = rev_tmp_idx[i];
	//		int to = adj[c][0][i];
	//		if(to != NUL && to !=CAP)
	//			to = rev_tmp_idx[to];
	//		if(to != NUL)
	//			fprintf(writer, "%d -- %d [color=%s];\n", from, to, "red");
	//		to = adj[c][1][i];
	//		if(to != NUL && to !=CAP)
	//			to = rev_tmp_idx[to];
	//		if(to != NUL)
	//			fprintf(writer, "%d -- %d [color=%s];\n", from, to, "blue");
	//	}
	//	fprintf(writer, "}\n");
	//	fclose(writer);
	//}
	//check conflicts
	int actual_comb =0;
	int *violate = (int*)malloc(sizeof(int)*possible_comb);
	for(i=0;i<possible_comb;i++){
		violate[i] = FALSE;
		for(j=0; j<new_v_size; j++){
			if((adj[i][0][j]!=NUL && adj[i][0][j]!=CAP && j!=adj[i][0][adj[i][0][j]])){
				printf("%d %d\n", j, adj[i][0][j]);
				violate[i] = TRUE;
				break;
			} 
			else if((adj[i][1][j]!=NUL && adj[i][1][j]!=CAP && j!=adj[i][1][adj[i][1][j]])){
				printf("%d %d\n", j, adj[i][0][j]);
				violate[i] = TRUE;
				break;
			}
		}
		if(violate[i] == FALSE)
			actual_comb++;
	}
	//for(i=0;i<possible_comb;i++)
	//	printf("%d %s\n", i, violate[0]==TRUE?"TRUE":"FALSE");
	//exit(1);
	//compute distances 
	int *vet_type = (int*)malloc(sizeof(int)*new_v_size);
	//final computation
	int min = 1000000000;
	int min_idx = 0;
	for(j=0;j<possible_comb;j++){
		if(violate[j]==FALSE){
			for(i=0;i<new_v_size; i++)
			{
				if(adj[j][1][i]!=NUL && adj[j][1][i]!=CAP && adj[j][0][i]==NUL)
					vet_type[i]=PI_OPEN;
				else if(adj[j][0][i]!=NUL && adj[j][0][i]!=CAP && adj[j][1][i]==NUL)
					vet_type[i]=GAMMA_OPEN;
				else if(adj[j][0][i]!=NUL && adj[j][1][i]!=NUL &&
						adj[j][0][i]!=CAP && adj[j][1][i]!=CAP)
					vet_type[i]=CYCLE;
				else if((adj[j][0][i]!=NUL && adj[j][0][i]==CAP) ||
						(adj[j][1][i]!=NUL && adj[j][1][i]==CAP))
					vet_type[i]=END;
				else
					vet_type[i]=NONE;
			}
			int results = compute_indel_dis(adj[j], vet_type, 0, 1, new_v_size, gene);
			if(results < min){
				min = results;
				min_idx = i;
			}
		}
	}
	//free graph
	if(IS_TWO == FALSE)
		free_graph(g);
	else
		free_graph_two(g);
	//free comp_***
	free(comp_v_idx[0]);
	free(comp_v_idx[1]);
	free(comp_v_idx);
	free(comp_v_deg[0]);
	free(comp_v_deg[1]);
	free(comp_v_deg);
	free(comp_e_idx[0]);
	free(comp_e_idx[1]);
	free(comp_e_idx);
	free(comp_v_check);
	//free adj
	for(i=0;i<possible_comb; i++){
		free(adj[i][0]);
		free(adj[i][1]);
		free(adj[i]);
	}
	free(adj);
	free(vet_type);
	free(violate);
	return min;
}

//this function is used to compact the graph
void 
compress_graph(int **comp_v_idx, int **comp_e_idx, 
		int *comp_v_check, int **v_deg,
		int size, int *new_v_size, 
		int *new_e_size1, int *new_e_size2){
	int i, j;

	//get the degrees
	for(i=0;i<size;i++){
		int start_one = i==0?0:comp_v_idx[0][i-1];
		int start_two = i==0?0:comp_v_idx[1][i-1];
		int end_one = comp_v_idx[0][i];
		int end_two = comp_v_idx[1][i];
		if((end_one-start_one) <= 0)
			v_deg[0][i]=0;
		else if(comp_v_check[i] == NUL)
			v_deg[0][i]=0;
		else
			v_deg[0][i]=(end_one - start_one);
		if((end_two-start_two) <= 0)
			v_deg[1][i]=0;
		else if(comp_v_check[i] == NUL)
			v_deg[1][i]=0;
		else
			v_deg[1][i]=(end_two - start_two);
		if(v_deg[0][i]==0 && v_deg[1][i]==0)
			comp_v_check[i]=NUL;
	}
	//compressing
	int has_inner = TRUE;
	while(has_inner == TRUE){
		has_inner = FALSE;
		for(i=0; i<size; i++){
			if(comp_v_check[i]==NUL)
				continue;
			int one = i;
			if(v_deg[0][one]==1 && v_deg[1][one]==1){
				int two = comp_e_idx[0][one==0?0:comp_v_idx[0][one-1]];
				if(two != CAP && v_deg[0][two]==1 && v_deg[1][two]==1){
					int three = comp_e_idx[1][one==0?0:comp_v_idx[1][one-1]];
					int four = comp_e_idx[1][two==0?0:comp_v_idx[1][two-1]];
					if(three != CAP && four != CAP && 
							three != two && four != one &&
							v_deg[1][three]==1 && v_deg[1][four]==1){
						int start_three = three==0?0:comp_v_idx[1][three-1];
						int end_three = comp_v_idx[1][three];
						int start_four = four==0?0:comp_v_idx[1][four-1];
						int end_four = comp_v_idx[1][four];
						for(j=start_three; j<end_three; j++)
							if(comp_e_idx[1][j]==one)
								comp_e_idx[1][j] = four;
						for(j=start_four; j<end_four; j++)
							if(comp_e_idx[1][j]==two)
								comp_e_idx[1][j] = three;
						has_inner = TRUE;
						v_deg[0][one]=0;
						v_deg[1][one]=0;
						comp_e_idx[0][one==0?0:comp_v_idx[0][one-1]] = NUL;
						comp_e_idx[1][one==0?0:comp_v_idx[1][one-1]] = NUL;
						v_deg[0][two]=0;
						v_deg[1][two]=0;
						comp_e_idx[0][two==0?0:comp_v_idx[0][two-1]] = NUL;
						comp_e_idx[1][two==0?0:comp_v_idx[1][two-1]] = NUL;
						comp_v_check[one] = NUL;
						comp_v_check[two] = NUL;
					}
				}
			}
		}
	}
	has_inner = TRUE;
	while(has_inner == TRUE){
		has_inner = FALSE;
		for(i=0; i<size; i++){
			if(comp_v_check[i]==NUL)
				continue;
			int one = i;
			if(v_deg[1][one]==1 && v_deg[0][one]==1){
				int two = comp_e_idx[1][one==0?0:comp_v_idx[1][one-1]];
				if(two != CAP && v_deg[1][two]==1 && v_deg[0][two]==1){
					int three = comp_e_idx[0][one==0?0:comp_v_idx[0][one-1]];
					int four = comp_e_idx[0][two==0?0:comp_v_idx[0][two-1]];
					if(three != CAP && four != CAP && 
							three != two && four != one &&
							v_deg[0][three]==1 && v_deg[0][four]==1){
						int start_three = three==0?0:comp_v_idx[0][three-1];
						int end_three = comp_v_idx[0][three];
						int start_four = four==0?0:comp_v_idx[0][four-1];
						int end_four = comp_v_idx[0][four];
						for(j=start_three; j<end_three; j++)
							if(comp_e_idx[0][j]==one)
								comp_e_idx[0][j] = four;
						for(j=start_four; j<end_four; j++)
							if(comp_e_idx[0][j]==two)
								comp_e_idx[0][j] = three;
						has_inner = TRUE;
						v_deg[1][one]=0;
						v_deg[0][one]=0;
						comp_e_idx[1][one==0?0:comp_v_idx[1][one-1]] = NUL;
						comp_e_idx[0][one==0?0:comp_v_idx[0][one-1]] = NUL;
						v_deg[1][two]=0;
						v_deg[0][two]=0;
						comp_e_idx[1][two==0?0:comp_v_idx[1][two-1]] = NUL;
						comp_e_idx[0][two==0?0:comp_v_idx[0][two-1]] = NUL;
						comp_v_check[one] = NUL;
						comp_v_check[two] = NUL;
					}
				}
			}
		}
	}
	//FILE *writer = fopen("vis1_com.dot", "w"); 
	//fprintf(writer, "Graph{\n");
	//for(i=0;i<size;i++){
	//	if(comp_v_check[i]==NUL)
	//		continue;
	//	int s = i==0?0:comp_v_idx[0][i-1];
	//	int end = comp_v_idx[0][i];
	//	for(j=s;j<end;j++){
	//		int to = comp_e_idx[0][j];
	//		//if(to != NUL)
	//			fprintf(writer, "%d -- %d [color=%s];\n", i, to, "red");
	//	}
	//	s = i==0?0:comp_v_idx[1][i-1];
	//	end = comp_v_idx[1][i];
	//	for(j=s;j<end;j++){
	//		int to = comp_e_idx[1][j];
	//		//if(to != NUL)
	//			fprintf(writer, "%d -- %d [color=%s];\n", i, to, "blue");
	//	}
	//}
	//fprintf(writer, "}\n");
	//fclose(writer);
	//exit(1);
	//renaming
	int *tmp_idx = (int*)malloc(sizeof(int)*size);
	int *rev_tmp_idx = (int*)malloc(sizeof(int)*size);
	int sum=0;
	for(i=0;i<size;i++){
		if(comp_v_check[i]==NUL)
			continue;
		tmp_idx[i] = sum++;
		rev_tmp_idx[tmp_idx[i]] = i;
	}
	int **tmp_v_idx = (int**)malloc(sizeof(int*)*2); 
	int **tmp_e_idx = (int**)malloc(sizeof(int*)*2); 
	int **tmp_v_deg = (int**)malloc(sizeof(int*)*2);
	tmp_v_idx[0] = (int*)malloc(sizeof(int)*size); 
	tmp_v_idx[1] = (int*)malloc(sizeof(int)*size); 
	tmp_v_deg[0] = (int*)malloc(sizeof(int)*size); 
	tmp_v_deg[1] = (int*)malloc(sizeof(int)*size); 
	tmp_e_idx[0] = (int*)malloc(sizeof(int)*comp_v_idx[0][size-1]); 
	tmp_e_idx[1] = (int*)malloc(sizeof(int)*comp_v_idx[1][size-1]); 
	int sum1 = 0;
	int sum2 = 0;
	int e_size1 = 0;
	int e_size2 = 0;
	for(i=0; i<size; i++){
		if(comp_v_check[i]!=NUL){
			sum1 += v_deg[0][i];
			tmp_v_idx[0][tmp_idx[i]] = sum1;
			tmp_v_deg[0][tmp_idx[i]] = v_deg[0][i];
			int start = i==0?0:comp_v_idx[0][i-1];
			int end = comp_v_idx[0][i];
			for(j=start; j<end; j++){
				tmp_e_idx[0][e_size1] = comp_e_idx[0][j]==CAP?CAP:tmp_idx[comp_e_idx[0][j]];
				e_size1++;
			}
			sum2 += v_deg[1][i];
			tmp_v_idx[1][tmp_idx[i]] = sum2;
			tmp_v_deg[1][tmp_idx[i]] = v_deg[1][i];
			start = i==0?0:comp_v_idx[1][i-1];
			end = comp_v_idx[1][i];
			for(j=start; j<end; j++){
				tmp_e_idx[1][e_size2] = comp_e_idx[1][j]==CAP?CAP:tmp_idx[comp_e_idx[1][j]];
				e_size2++;
			}
		}
	}
	*new_v_size = sum;
	*new_e_size1 = e_size1;
	*new_e_size2 = e_size2;
	//copy back
	for(i=0;i<*new_v_size; i++){
		comp_v_idx[0][i] = tmp_v_idx[0][i];
		comp_v_idx[1][i] = tmp_v_idx[1][i];
		v_deg[0][i] = tmp_v_deg[0][i];
		v_deg[1][i] = tmp_v_deg[1][i];
	}
	//for(i=0;i<*new_v_size; i++)
	//	printf("%d|%d ", v_deg[0][i], v_deg[1][i]);
	for(i=0;i<*new_e_size1; i++){
		comp_e_idx[0][i] = tmp_e_idx[0][i];
	}
	for(i=0;i<*new_e_size2; i++){
		comp_e_idx[1][i] = tmp_e_idx[1][i];
	}

	//FILE *writer = fopen("vis1_com.dot", "w"); 
	//fprintf(writer, "Graph{\n");
	//for(i=0;i<*new_v_size;i++){
	//	int s = i==0?0:tmp_v_idx[0][i-1];
	//	int end = tmp_v_idx[0][i];
	//	for(j=s;j<end;j++){
	//		int to = tmp_e_idx[0][j];
	//		if(to != NUL)
	//			fprintf(writer, "%d -- %d [color=%s];\n", i, to, "red");
	//	}
	//	s = i==0?0:tmp_v_idx[1][i-1];
	//	end = tmp_v_idx[1][i];
	//	for(j=s;j<end;j++){
	//		int to = tmp_e_idx[1][j];
	//		if(to != NUL)
	//			fprintf(writer, "%d -- %d [color=%s];\n", i, to, "blue");
	//	}
	//}
	//fprintf(writer, "}\n");
	//fclose(writer);
	//exit(1);
	//free functions
	free(tmp_v_idx[0]); 
	free(tmp_v_idx[1]); 
	free(tmp_v_deg[0]); 
	free(tmp_v_deg[1]); 
	free(tmp_e_idx[0]); 
	free(tmp_e_idx[1]); 
	free(tmp_v_idx); 
	free(tmp_e_idx); 
	free(tmp_v_deg);
	free(tmp_idx);
	free(rev_tmp_idx);
}

int 
compute_all_possible_exemplar_orders(po o, int c1, 
		int c2, int gene){
	int i, j, k;
	pg g1 = (pg)malloc(sizeof(tg));
	g1->v_size = o->v_size;
	g1->e_size = (int*)malloc(sizeof(int)*2);
	g1->v_check = (int*)malloc(sizeof(int)*g1->v_size);
	for(i=0; i<g1->v_size; i++)
		g1->v_check[i] = TRUE;

	g1->v_idx = (int**)malloc(sizeof(int*)*2); 
	g1->e_idx = (int**)malloc(sizeof(int*)*2);
	g1->e_check = (int**)malloc(sizeof(int*)*2);
	g1->v_deg = (int**)malloc(sizeof(int*)*2);

	for(i=0;i<2;i++){
		int c=i==0?c1:c2;
		g1->v_idx[i] = (int*)malloc(sizeof(int)*g1->v_size); 
		g1->v_deg[i] = (int*)malloc(sizeof(int)*g1->v_size);
		g1->e_size[i] = o->v_idx[c][g1->v_size-1];
		g1->e_idx[i] = (int*)malloc(sizeof(int)*g1->e_size[i]);
		g1->e_check[i] = (int*)malloc(sizeof(int)*g1->e_size[i]);
		for(j=0; j<g1->v_size; j++){
			g1->v_idx[i][j] = o->v_idx[c][j];
			int start = j==0?0:g1->v_idx[i][j-1];
			int end = g1->v_idx[i][j];
			g1->v_deg[i][j] = end - start;
		}
		for(j=0;j<g1->e_size[i];j++){
			g1->e_idx[i][j] = o->e_idx[c][j];
			g1->e_check[i][j] = TRUE;
		}
	}
	IS_TWO = TRUE;
	//graph_vis_two(g1, "0_23.dot", 0, 1);
	int min = g1->v_size/2 - cpec_mem_save(g1, 0, 1, g1->v_size);
	//int min = compute_all_possible_exemplar_compact(g1, 0, 1, g1->v_size);	
	free_graph_two(g1);
	return min;
}

int compute_all_possible_exemplar(pg g1, int c1, 
		int c2, int gene){
	pg g = (pg)malloc(sizeof(tg));
	if(IS_TWO == FALSE){
		prepare_graph_for_distance(g1, g);
	}
	else{
		prepare_graph_for_distance_two(g1, g);
	}
	int i,j,k;
	int possible_comb=1;
	int color;
	for(i=0; i<g->v_size; i++){
		if(g->v_check[i] != NUL && g->v_deg[c1][i] != 0 && g->v_deg[c2][i]!=0)
			possible_comb = possible_comb*(g->v_deg[c1][i]*g->v_deg[c2][i]);
	}
	graph_vis_two(g, "vis.dot", c1, c2);
	int *results = (int*)malloc(sizeof(int)*possible_comb);
	int ***adj = (int***)malloc(sizeof(int**)*possible_comb);
	for(i=0;i<possible_comb;i++){
		adj[i] = (int**)malloc(sizeof(int*)*2);
		adj[i][0] = (int*)malloc(sizeof(int)*g->v_size);
		adj[i][1] = (int*)malloc(sizeof(int)*g->v_size);
	}
	int period =0;
	for(i=0;i<g->v_size;i++){
		period = 0;
		if(g->v_check[i]!=NUL){
			if(g->v_deg[c1][i]!=0 &&g->v_deg[c2][i]!=0){
				period = g->v_deg[c1][i]*g->v_deg[c2][i];
				int **period_arr = (int**)malloc(sizeof(int*)*2);
				period_arr[0] = (int*)malloc(sizeof(int)*period);
				period_arr[1] = (int*)malloc(sizeof(int)*period);
				int start1 = i==0?0:g->v_idx[c1][i-1];
				int start2 = i==0?0:g->v_idx[c2][i-1];
				int end1 = g->v_idx[c1][i];
				int end2 = g->v_idx[c2][i];
				int pos=0;
				for(j=start1;j<end1;j++){
					for(k=start2;k<end2;k++){	
						period_arr[0][pos] = g->e_idx[c1][j];
						period_arr[1][pos] = g->e_idx[c2][k];
						pos++;
					}
				}
				for(j=0;j<possible_comb;j++){
					adj[j][0][i] = period_arr[0][j%period];
					adj[j][1][i] = period_arr[1][j%period];
				}
			}
			else{
				int **period_arr = (int**)malloc(sizeof(int*)*2);
				if(g->v_deg[c1][i]==0 && g->v_deg[c2][i]==0)
					;
				else if(g->v_deg[c1][i]==0 && g->v_deg[c2][i]!=0){
					period = g->v_deg[c2][i];
					period_arr[0] = (int*)malloc(sizeof(int)*period);
					period_arr[1] = (int*)malloc(sizeof(int)*period);
					int start2 = i==0?0:g->v_idx[c2][i-1];
					int end2 = g->v_idx[c2][i];
					int pos=0;
					for(k=start2;k<end2;k++){
						period_arr[0][pos] = NUL;
						period_arr[1][pos] = g->e_idx[c2][k];
						pos++;
					}
				}
				else if(g->v_deg[c1][i]!=0 && g->v_deg[c2][i]==0){
					period = g->v_deg[c1][i];
					period_arr[0] = (int*)malloc(sizeof(int)*period);
					period_arr[1] = (int*)malloc(sizeof(int)*period);
					int start1 = i==0?0:g->v_idx[c1][i-1];
					int end1 = g->v_idx[c1][i];
					int pos=0;
					for(k=start1;k<end1;k++){
						period_arr[0][pos] = g->e_idx[c1][k];
						period_arr[1][pos] = NUL;
						pos++;
					}
				}
				if(period != 0){
					for(j=0;j<possible_comb;j++){
						adj[j][0][i] = period_arr[0][j%period];
						adj[j][1][i] = period_arr[1][j%period];
					}
				}
				else{
					for(j=0;j<possible_comb;j++){
						adj[j][0][i] = NUL;
						adj[j][1][i] = NUL;
					}
				}
			}
		}
		else{
			for(j=0;j<possible_comb;j++){
				adj[j][0][i] = NUL;
				adj[j][1][i] = NUL;
			}
		}
	}
	//check
	for(i=0;i<possible_comb;i++){
		for(j=0;j<g->v_size;j++){
			int to = adj[i][0][j];
			if(to!=NUL && to != CAP){
				if(adj[i][0][to] != j)
					adj[i][0][j]=NUL;
			}
			to = adj[i][1][j];
			if(to!=NUL && to != CAP){
				if(adj[i][1][to] != j)
					adj[i][1][j]=NUL;
			}
		}
	}
	//check conflicts
	int actual_comb =0;
	int *violate = (int*)malloc(sizeof(int)*possible_comb);
	for(i=0;i<possible_comb;i++){
		violate[i] = FALSE;
		for(j=0; j<g->v_size; j++){
			if((adj[i][0][j]!=NUL && adj[i][0][j]!=CAP && j!=adj[i][0][adj[i][0][j]])){
				violate[i] = TRUE;
				break;
			} 
			else if((adj[i][1][j]!=NUL && adj[i][1][j]!=CAP && j!=adj[i][1][adj[i][1][j]])){
				violate[i] = TRUE;
				break;
			}
		}
		if(violate[i] == FALSE)
			actual_comb++;
	}
	//int c;
	//for(c=0; c<possible_comb; c++){
	//        char f_buf[100];
	//        sprintf(f_buf, "vis_c%d.dot", c);
	//        FILE *writer = fopen(f_buf, "w");
	//        fprintf(writer, "Graph{\n");
	//        for(i=0;i<g->v_size; i++){
	//		if(g->v_check[i] != NUL){
	//			int to = adj[c][0][i];
	//			if(to != NUL)
	//				fprintf(writer, "%d -- %d [color=%s];\n", i, to, "red");
	//			to = adj[c][1][i];
	//			if(to != NUL)
	//				fprintf(writer, "%d -- %d [color=%s];\n", i, to, "blue");
	//		}
	//	}
	//	fprintf(writer, "}\n");
	//	fclose(writer);
	//}
	//get vet type
	int **vet_type = (int**)malloc(sizeof(int*)*possible_comb);
	for(i=0;i<possible_comb;i++)
		vet_type[i] = (int*)malloc(sizeof(int)*g->v_size);
	for(j=0;j<possible_comb;j++){
		if(violate[j] == FALSE)
			for(i=0;i<g->v_size; i++)
			{
				if(adj[j][1][i]!=NUL && adj[j][1][i]!=CAP && adj[j][0][i]==NUL)
					vet_type[j][i]=PI_OPEN;
				else if(adj[j][0][i]!=NUL && adj[j][0][i]!=CAP && adj[j][1][i]==NUL)
					vet_type[j][i]=GAMMA_OPEN;
				else if(adj[j][0][i]!=NUL && adj[j][1][i]!=NUL && 
						adj[j][0][i]!=CAP && adj[j][1][i]!=CAP)
					vet_type[j][i]=CYCLE;
				else if((adj[j][0][i]!=NUL && adj[j][0][i]==CAP) ||
						(adj[j][1][i]!=NUL && adj[j][1][i]==CAP))
					vet_type[j][i]=END;
				else
					vet_type[j][i]=NONE;	
			}
	}
	//start computing
	int min = 1000000000;
	int min_idx = 0;
	for(i=0;i<possible_comb;i++){
		if(violate[i]==FALSE){
			results[i] = compute_indel_dis(adj[i], vet_type[i], 0, 1, g->v_size, gene);
			if(results[i] < min){
				min = results[i];
				min_idx = i;
			}
		}
	}
	//printf("the min is %d\n", min);
	return min;
}

void to_adj_arr(pg g, int **adj, int *vet_type, int c1, int c2){
	int i,j;
	int color;
	for(color=0;color<2;color++){
		int c = color==0?c1:c2;
		for(i=0;i<g->v_size; i++){
			if(g->v_check[i] != NUL){
				int start = i==0?0:g->v_idx[c][i-1];
				int end = g->v_idx[c][i];
				for(j=start;j<end;j++){
					if(g->e_check[c][j]==TRUE)
						adj[c][i] = g->e_idx[c][j];
				}
			}
			else
				adj[c1][i] = NUL;
		}
	}
	//trying to get the type
	for(i=0;i<g->v_size; i++)
	{
		if(adj[c2][i]!=NUL && adj[c2][i]!=CAP && adj[c1][i]==NUL)
			vet_type[i]=PI_OPEN;
		else if(adj[c1][i]!=NUL && adj[c1][i]!=CAP && adj[c2][i]==NUL)
			vet_type[i]=GAMMA_OPEN;
		else if(adj[c1][i]!=NUL && adj[c2][i]!=NUL && 
				adj[c1][i]!=CAP && adj[c2][i]!=CAP)
			vet_type[i]=CYCLE;
		else if(adj[c1][i]!=NUL && adj[c2][i]!=NUL &&
				adj[c1][i]==CAP && adj[c2][i]==CAP)
			vet_type[i]=END;
		else
			vet_type[i]=NONE;	
	}
}

int compute_indel_dis(int **adj, int *vet_type, int c1, int c2, int size, int gene){
	//vis_adj("vis_adj.dot", adj, size);
	int used[MAX_VET_SIZE][2];
	int left = 0;
	int right = 0;

	int num_c=0;
	int num_even=0;
	int num_pi_pi=0;
	int num_gamma_gamma=0;
	int num_pi_gamma=0;
	int num_reg_even=0;
	int num_pi_odd=0;
	int num_pi_even=0;
	int num_gamma_odd=0;
	int num_gamma_even=0;
	int delta=0;

	int i,j;
	//open path
	for(i=0;i<size;i++)
	{
		if(vet_type[i]!=NONE && vet_type[i]!=CYCLE && vet_type[i]!= END)
		{
			if(vet_type[i]==PI_OPEN)
			{
				left = i;
				vet_type[left]=NONE;
				set_visited(vet_type, left);
				while(1)
				{
					right = adj[c2][left];
					if(vet_type[right]==PI_OPEN){
						num_pi_pi++;
						vet_type[right]=NONE;
						set_visited(vet_type, right);
						break;
					}
					if(vet_type[right]==END){
						num_pi_odd++;
						vet_type[right]=NONE;
						set_visited(vet_type, right);
						break;
					}
					vet_type[right]=NONE;
					set_visited(vet_type, right);

					left=adj[c1][right];
					if(vet_type[left]==GAMMA_OPEN){
						num_pi_gamma++;
						vet_type[left] = NONE;
						set_visited(vet_type, left);
						break;
					}
					if(vet_type[left]==END){
						num_pi_even++;
						vet_type[left] = NONE;
						set_visited(vet_type, left);
						break;
					}
					vet_type[left] = NONE;
					set_visited(vet_type, left);
				}
			}
			else if(vet_type[i]==GAMMA_OPEN)
			{
				left = i;
				vet_type[left]=NONE;
				set_visited(vet_type, left);
				while(1)
				{
					right = adj[c1][left];
					if(vet_type[right]==GAMMA_OPEN){
						num_gamma_gamma++;
						vet_type[right]=NONE;
						set_visited(vet_type, right);
						break;
					}
					if(vet_type[right]==END){
						num_gamma_odd++;
						vet_type[right]=NONE;
						set_visited(vet_type, right);
						break;
					}
					vet_type[right]=NONE;
					set_visited(vet_type, right);

					left=adj[c2][right];
					if(vet_type[left]==PI_OPEN){
						num_pi_gamma++;
						vet_type[left] = NONE;
						set_visited(vet_type, left);
						break;
					}
					if(vet_type[left]==END){
						num_gamma_even++;
						vet_type[left] = NONE;
						set_visited(vet_type, left);
						break;
					}
					vet_type[left] = NONE;
					set_visited(vet_type, left);
				}
			}
		}
	}
	//calculate clean paths
	int tmp_c1;
	int tmp_c2;
	for(i=0;i<size;i++)
	{
		if(vet_type[i]!=NONE && vet_type[i]!=CYCLE)
		{
			left = i;
			//if(vet_type[left]==NONE){
			//	printf("%d has visited\n", left);
			//	ERROR_PRINT();
			//}
			vet_type[left]=NONE;
			set_visited(vet_type, left);
			if(adj[c1][left]==CAP && adj[c2][left]!=CAP)
			{
				tmp_c1=1;
				tmp_c2=0;
			}
			else if(adj[c1][left]!=CAP && adj[c2][left]==CAP)
			{
				tmp_c1=0;
				tmp_c2=1;
			}
			else if(adj[c1][left]==CAP && adj[c2][left]==CAP)
			{
				num_even++;
				continue;
			}
			while(1)
			{
				right = adj[tmp_c1][left];
				if(right == NUL){
					break;
				}
				if(vet_type[right]==END){
					vet_type[right]=NONE;
					set_visited(vet_type, right);
					break;
				}
				vet_type[right]=NONE;
				set_visited(vet_type, right);

				left=adj[tmp_c2][right];
				if(vet_type[left]==END){
					num_even++;
					vet_type[left] = NONE;
					set_visited(vet_type, left);
					break;
				}
				vet_type[left] = NONE;
				set_visited(vet_type, left);
			}
		}
	}

	//calculate cycles
	for(i=0;i<size;i++)
	{
		if(vet_type[i]!=NONE)
		{
			int start = i;
			left = i;
			vet_type[left]=NONE;
			set_visited(vet_type, left);
			while(1)
			{
				right = adj[c1][left];
				if(right==start){
					num_c++;
					vet_type[right]=NONE;
					set_visited(vet_type, right);
					break;
				}
				vet_type[right]=NONE;

				left=adj[c2][right];
				if(left==start){
					num_c++;
					vet_type[left] = NONE;
					set_visited(vet_type, left);
					break;
				}
				vet_type[left] = NONE;
				set_visited(vet_type, left);
			}
		}
	}
	//calculate delta
	if((num_pi_gamma%2==0)&&
			((num_pi_odd>num_pi_even && num_gamma_odd>num_gamma_even)||
			 (num_pi_odd<num_pi_even && num_gamma_odd<num_gamma_even)))
		delta=1;
	else
		delta=0;

	//calculate result
	int result = 
		((num_c+num_pi_pi+num_gamma_gamma+(int)floor((double)num_pi_gamma/2))
		 +(num_even+delta+
			 num_pi_odd>num_pi_even?num_pi_even:num_pi_odd+
			 num_gamma_odd>num_gamma_even?num_gamma_even:num_gamma_odd)/2);
	//printf("the min is %d\n", result);
	return result;	
}

void set_visited(int *vet_type, int pos){
	vet_type[pos] = NONE;
	//printf("%d visisted\n", pos);
}

int evaluate(pg g){
	int sum =0; 
	int dis1, dis2, dis3;
	int effective_size=0, i, j, color;
	for(i=0;i<g->v_size;i++)
		if(g->v_check[i] != NUL)
			effective_size++;
	//for(color=0;color<4;color++){
	//	for(i=0;i<g->v_size; i++){
	//		int start = i==0?0:g->v_idx[color][i-1];
	//		int end = g->v_idx[color][i];
	//		for(j=start;j<end;j++){
	//			printf("%d:%d->%d ", j, i, g->e_idx[color][j]);
	//		}
	//	}
	//	printf("\n");
	//}
	prepare_dis(g);
	IS_TWO = FALSE;
#ifdef USE_COMPRESS
	g->distance[0] = dis1 = compute_all_possible_exemplar_compact(g, 0, 3, effective_size);
	g->distance[1] = dis2 = compute_all_possible_exemplar_compact(g, 1, 3, effective_size);
	g->distance[2] = dis3 = compute_all_possible_exemplar_compact(g, 2, 3, effective_size);
	
#elif MEM_SAVE
	g->distance[0] = dis1 = effective_size/2 - cpec_mem_save(g, 0, 3, effective_size);
	g->distance[1] = dis2 = effective_size/2 - cpec_mem_save(g, 1, 3, effective_size);
	g->distance[2] = dis3 = effective_size/2 - cpec_mem_save(g, 2, 3, effective_size);
#else
	g->distance[0] = dis1 = compute_all_possible_exemplar(g, 0, 3, effective_size);
	g->distance[1] = dis2 = compute_all_possible_exemplar(g, 1, 3, effective_size);
	g->distance[2] = dis3 = compute_all_possible_exemplar(g, 2, 3, effective_size);
#endif
	//printf("dis1 %d dis2 %d dis3 %d gene %d cycle %d\n", dis1, dis2, dis3, effective_size/2, g->cycle);
	sum = dis1+dis2+dis3-g->cycle;
	return sum;
}

void vis_adj(char *file, int **adj, int size){
	int i;
	FILE *writer = fopen(file, "w");
	fprintf(writer, "Graph{\n");
	for(i=0;i<size; i++){
		if(adj[0][i]!=NUL)
			fprintf(writer, "%d -- %d [color=red];\n", i, adj[0][i]);
		if(adj[1][i]!=NUL)
			fprintf(writer, "%d -- %d [color=blue];\n", i, adj[1][i]);
	}
	fprintf(writer, "}\n");
	fclose(writer);
}

void vis_compact_adj(char *file, int **comp_v_idx,
		int **comp_e_idx, int v_size){
	int i, j;	
	FILE *writer = fopen(file, "w");
	if(writer == NULL){
		printf("the file %s does not exists\n", file);
		exit(1);
	}
	fprintf(writer, "Graph{\n");
	for(i=0; i<v_size; i++){
		int start = i==0?0:comp_v_idx[0][i-1];
		int end = comp_v_idx[0][i];
		for(j=start; j<end; j++){
			int v_to = comp_e_idx[0][j];
			//if(i<v_to)
				fprintf(writer, "%d -- %d [color=%s];\n", i, v_to, "red");
		}
		start = i==0?0:comp_v_idx[1][i-1];
		end = comp_v_idx[1][i];
		for(j=start; j<end; j++){
			int v_to = comp_e_idx[1][j];
			//if(i<v_to)
				fprintf(writer, "%d -- %d [color=%s];\n", i, v_to, "blue");
		}
	}
	fprintf(writer, "}\n");
	fclose(writer);
}

#include "sim.h"
#include "list.h"
#include <math.h>

void 
init_orders_sim(po o, int gene_len, int gnm_num);
void 
mutate(po o, int pos, double theta);
void
rmv(po o, int pos, double gamma);
void
duplicate(po o, int pos, double phi);
int
fork_sim(po o, int pos, 
		double theta, double gamma, double phi);
void 
set_state(po o);
void 
compute_num_e(po o, int pos);
void
csr_to_adj(po o, pa a, int pos);
void
adj_to_csr(po o, pa a, int pos);
void
vis_order(po o, int pos_n);
void 
vis_fork_adj(pa a, int pos_n);
void
vis_order_step(po o, int pos, char *file);

int total_score;
/* **********************************************************************
 * this is the main function to generate gene order files based on a tree
 * **********************************************************************/
void 
run_simulate(char *outfile, int gene_len, int gnm_num, 
		double theta, double gamma, double phi){
	int i, j, k;
	total_score = 0;
	int effective_edge =0;
	po o = (po) malloc(sizeof(to));
	init_orders_sim(o, gene_len, gnm_num);
	compute_num_e(o, 0);
	int fork_num = log2(gnm_num);
	printf("we need to fork %d iterations\n", fork_num);
	for(i=0; i<fork_num; i++){
		int num_g = o->num_genome;
		for(j=0; j<num_g; j++){
			if(o->valid[j] == VALID2){
				effective_edge += fork_sim(o, j, theta, gamma, phi);
			}
		}
		set_state(o);
	}
	//output the results
	effective_edge=0;
	for(i=0; i<o->num_genome; i++){
		if(o->valid[i] != VALID2)
			continue;
		for(j=0;j<o->v_size; j++){
			int start = j==0?0:o->v_idx[i][j-1];
			int end = o->v_idx[i][j];
			effective_edge+=(end-start);
		}
	}
	FILE *writer = fopen(outfile, "w");
	fprintf(writer, "%d %d %d\n", o->v_size, gnm_num, effective_edge);
	int gnm_id=1;
	for(i=0; i<o->num_genome; i++){
		if(o->valid[i] != VALID2)
			continue;
		for(j=0;j<o->v_size; j++){
			int start = j==0?0:o->v_idx[i][j-1];
			int end = o->v_idx[i][j];
			for(k=start; k<end; k++)
				if(o->e_idx[i][k]!=NUL)
					fprintf(writer, "%d %d %d\n", j, o->e_idx[i][k], gnm_id);
		}
		gnm_id++;
	}
	fprintf(writer, "score: %d\n", total_score);
	fclose(writer);
}

/* *******************************************************
 * this is to init a seed genome and malloc memory for the 
 * genomes in the branchs
 * *******************************************************/
void 
init_orders_sim(po o, int gene_len, int gnm_num){
	int i, j, k;
	//set up basic information
	o->num_genome = 1;
	o->num_med_genome = -1;
	int result = gnm_num;
	o->num_input_genome = 0;
	while(result>0){
		o->num_input_genome += result;
		result = result/2;
	}
	o->med_genome_idx = 1;
	o->v_size = gene_len*2;
	//allocate memory
	o->e_size = (int*)malloc(sizeof(int)*o->num_input_genome);
	o->valid = (int*)malloc(sizeof(int)*o->num_input_genome);
	o->v_idx = (int**)malloc(sizeof(int*)*o->num_input_genome);
	o->e_idx = (int**)malloc(sizeof(int*)*o->num_input_genome);
	for(i=0;i<o->num_input_genome; i++){
		o->v_idx[i] = (int*)malloc(sizeof(int)*o->v_size);
		o->e_idx[i] = (int*)malloc(sizeof(int)*o->v_size*2);
		o->valid[i] = -1;
	}
	//init seed genome
	o->valid[0] = VALID2;
	for(i=0;i<o->v_size; i++)
		o->v_idx[0][i] = i+1;
	for(i=1;i<o->v_size-1; i+=2){
		o->e_idx[0][i] = i+1;
		o->e_idx[0][i+1] = i;
	}
	o->e_idx[0][0] = CAP;
	o->e_idx[0][o->v_size-1] = CAP;
}

/* *******************************************************
 * after forking turning the current active node into dead 
 * and tunring new generated node into active.
 * *******************************************************/
void 
set_state(po o){
	int i;
	for(i=0; i<o->num_genome; i++){
		if(o->valid[i] == VALID2)
			o->valid[i] = VALID1;
		if(o->valid[i] == VALID3)
			o->valid[i] = VALID2;
	}
}


/* **************************************************
 * from a gvien node, fork two new nodes,
 * and perform operations such as DCJ, indel, and dup
 * on the new generated genomes.
 * **************************************************/
int
fork_sim(po o, int pos, 
		double theta, double gamma, double phi){
	int i, j, k;
	//copy
	for(i=0; i<2; i++){
		int pos_n = o->med_genome_idx;
		printf("processing %dth node\n", pos_n);
		for(j=0; j<o->v_size; j++)
			o->v_idx[pos_n][j] = o->v_idx[pos][j];
		for(j=0; j<o->e_size[pos]; j++)
			o->e_idx[pos_n][j] = o->e_idx[pos][j];
		//perform mutations
		mutate(o, pos_n, theta);
		rmv(o, pos_n, gamma);
		duplicate(o, pos_n, phi);
		compute_num_e(o, pos_n);
		total_score += compute_all_possible_exemplar_orders(o, pos, pos_n, o->v_size/2);
		//change stats
		o->valid[pos_n] = VALID3;
		o->med_genome_idx++;
		o->num_genome++;
	}
	int result =0;
	for(i=0; i<o->num_genome; i++){
		if(o->valid[i] == VALID3)
			result += o->e_size[i];
	}
	//for(i=0; i<o->num_genome; i++)
	//	printf("%d ", o->valid[i]);
	//printf("\n");
	return result;
}

/* ********************************
 * compute how many edges there are
 * after performing an operation.
 * ********************************/
void 
compute_num_e(po o, int pos){
	int i, j;
	int count =0;
	for(i=0; i<o->v_size; i++){
		int start = i==0?0:o->v_idx[pos][i-1];
		int end = o->v_idx[pos][i];
		count+= (end-start);
	}
	o->e_size[pos] = count;
}

/* *************************************
 * performing DCJ operations on a genome
 * *************************************/
void 
mutate(po o, int pos, double theta){
	int i, j, k;
	int num_mut = (o->v_size/2) * theta;
	for(j=0; j<num_mut; j++){
		int from_one = rand() % o->v_size;
		int start = from_one == 0?0:o->v_idx[pos][from_one-1];
		int end = o->v_idx[pos][from_one];
		for(i=start; i<end; i++){
			if(o->e_idx[pos][i] == CAP){
				end = 0;
				start = 0;
			}
		}
		while((end-start) == 0){
			from_one = rand() % o->v_size;
			start = from_one == 0?0:o->v_idx[pos][from_one-1];
			end = o->v_idx[pos][from_one];
		}
		int to_one;
		for(i=start; i<end; i++){
			if(o->e_idx[pos][i]!=NUL && o->e_idx[pos][i]!=CAP){
				to_one = o->e_idx[pos][i];
				break;
			}
		}
		int from_two = rand() % o->v_size;
		start = from_two == 0?0:o->v_idx[pos][from_two-1];
		end = o->v_idx[pos][from_two];
		for(i=start; i<end; i++){
			if(o->e_idx[pos][i] == CAP){
				end = 0;
				start = 0;
			}
		}
		while((end-start) == 0){
			from_two = rand() % o->v_size;
			start = from_two == 0?0:o->v_idx[pos][from_two-1];
			end = o->v_idx[pos][from_two];
		}
		int to_two;
		for(i=start; i<end; i++)
			if(o->e_idx[pos][i]!=NUL && o->e_idx[pos][i]!=CAP){
				to_two = o->e_idx[pos][i];
				break;
			}
		int conn_type = rand() % 2;
		if(from_one != to_one && from_two != to_two 
				&& from_one != from_two && to_one != to_two
				&& from_one != to_two && from_two != to_one){
			//make change
			//if(pos==7){
			//	printf("from_one %d to_one %d from_two %d to_two %d %d\n", from_one, to_one, from_two, to_two, o->e_idx[pos][o->v_idx[pos][198]]);
			//}
			int start_from_one = from_one==0?0:o->v_idx[pos][from_one-1];
			int end_from_one = o->v_idx[pos][from_one]; 
			int start_to_one = to_one==0?0:o->v_idx[pos][to_one-1];
			int end_to_one = o->v_idx[pos][to_one];
			int start_from_two = from_two==0?0:o->v_idx[pos][from_two-1];
			int end_from_two = o->v_idx[pos][from_two];
			int start_to_two = to_two==0?0:o->v_idx[pos][to_two-1]; 
			int end_to_two = o->v_idx[pos][to_two];
			if(conn_type == CONN_TYPE_ONE){
				for(i=start_from_one;i<end_from_one;i++)
					if(o->e_idx[pos][i] == to_one){
						o->e_idx[pos][i] = from_two;
						break;
					}
				for(i=start_to_one;i<end_to_one;i++)
					if(o->e_idx[pos][i] == from_one){
						o->e_idx[pos][i] = to_two;
						break;
					}
				for(i=start_from_two;i<end_from_two;i++)
					if(o->e_idx[pos][i] == to_two){
						o->e_idx[pos][i] = from_one;
						break;
					}
				for(i=start_to_two;i<end_to_two;i++)
					if(o->e_idx[pos][i] == from_two){
						o->e_idx[pos][i] = to_one;
						break;
					}
			}
			else if(conn_type == CONN_TYPE_TWO){
				for(i=start_from_one;i<end_from_one;i++)
					if(o->e_idx[pos][i] == to_one){
						o->e_idx[pos][i] = to_two;
						break;
					}
				for(i=start_to_one;i<end_to_one;i++)
					if(o->e_idx[pos][i] == from_one){
						o->e_idx[pos][i] = from_two;
						break;
					}
				for(i=start_from_two;i<end_from_two;i++)
					if(o->e_idx[pos][i] == to_two){
						o->e_idx[pos][i] = to_one;
						break;
					}
				for(i=start_to_two;i<end_to_two;i++)
					if(o->e_idx[pos][i] == from_two){
						o->e_idx[pos][i] = from_one;
						break;
					}

			}
		}
		//check violate here
		//for(i=0; i<o->v_size; i++){
		//	int start = i==0?0:o->v_idx[pos][i-1];
		//	int end = o->v_idx[pos][i];
		//	for(k=start; k<end; k++){
		//	}
		//}
		//if(pos == 7){
		//	char f_buf[100];
		//	sprintf(f_buf, "vis/step/v%d.dot", j);
		//	vis_order_step(o, pos, f_buf);
		//}
	}
}


/* ***************************
 * perform the indel operation
 * ***************************/
void
rmv(po o, int pos, double gamma){
	int i, j, k;
	int gene_size = o->v_size / 2;
	int num_rmv = gene_size * gamma;
	for(j=0; j<num_rmv; j++){
		int which_gene = rand() % gene_size + 1;
		int head = 2*(which_gene-1);
		int head_to = NUL;
		int start = head==0?0:o->v_idx[pos][head-1];
		int end = o->v_idx[pos][head];
		for(i=start; i<end; i++){
			int to = o->e_idx[pos][i];
			if(to != CAP && to != NUL){
				head_to = to;
				o->e_idx[pos][i] = NUL;
				break;
			}
		}
		int tail = 2*(which_gene-1)+1;
		int tail_to = NUL;
		start = tail==0?0:o->v_idx[pos][tail-1];
		end = o->v_idx[pos][tail];
		for(i=start; i<end; i++){
			int to = o->e_idx[pos][i];
			if(to != CAP && to != NUL){
				tail_to = to;
				o->e_idx[pos][i] = NUL;
				break;
			}
		}
		if(head == NUL || tail == NUL)
			continue; //already been deleted
		//bridge out the vertex
		int start_head_to = head_to==0?0:o->v_idx[pos][head_to-1];
		int end_head_to = o->v_idx[pos][head_to];
		int start_tail_to = tail_to==0?0:o->v_idx[pos][tail_to-1];
		int end_tail_to = o->v_idx[pos][tail_to];
		for(i=start_head_to; i<end_head_to; i++)
			if(o->e_idx[pos][i] == head)
				o->e_idx[pos][i] = tail_to;
		for(i=start_tail_to; i<end_tail_to; i++)
			if(o->e_idx[pos][i] == tail)
				o->e_idx[pos][i] = head_to;
	}
	pa a = (pa)malloc(sizeof(ta));
	csr_to_adj(o, a, pos);
	adj_to_csr(o, a, pos);
}

/* *****************************************
 * this is to perform duplication operations
 * *****************************************/
void
duplicate(po o, int pos, double phi){
	int i, j, k;
	pa a = (pa)malloc(sizeof(ta));
	csr_to_adj(o, a, pos);
	int gene_size = o->v_size / 2;
	int num_dup = o->v_size * phi;
	for(j=0; j<num_dup; j++){
		//select gene to duplicate
		int dup_gene = rand() % gene_size + 1;
		int head = 2*(dup_gene-1);
		int tail = 2*(dup_gene-1) + 1;
		int num_dup = a->idx_list[head];
		for(i=0; i<num_dup; i++){
			//select pos to insert
			int insert_pos = rand() % o->v_size;
			int insert_to = a->list[insert_pos][0];
			while(a->idx_list[insert_pos] == 0 ||
					insert_to == CAP ||
					insert_to == NUL){
				insert_pos = rand() % o->v_size;
				insert_to = a->list[insert_pos][0];
			}
			if(insert_pos != head && insert_pos != tail
					&& insert_to != head && insert_to != tail){
				for(k=0; k<a->idx_list[insert_pos]; k++)
					if(a->list[insert_pos][k] == insert_to)
						a->list[insert_pos][k] = head;
				for(k=0; k<a->idx_list[insert_to]; k++)
					if(a->list[insert_to][k] == insert_pos)
						a->list[insert_to][k] = tail;
				a->list[head][a->idx_list[head]] = insert_pos;
				a->list[tail][a->idx_list[tail]] = insert_to;
				a->idx_list[head]++;
				a->idx_list[tail]++;
			}
		}
	}
	adj_to_csr(o, a, pos);
}

/* ***********************************************
 * turn csr into adj for the convinience of coding
 * ***********************************************/
void
csr_to_adj(po o, pa a, int pos){
	int i, j, k;
	a->list = (int**)malloc(sizeof(int*)*o->v_size);
	a->idx_list = (int*)malloc(sizeof(int)*o->v_size);
	for(i=0; i<o->v_size; i++){
		a->list[i] = (int*)malloc(sizeof(int)*o->v_size);
		a->idx_list[i] = 0;
	}
	a->list_size = o->v_size;
	for(i=0; i<o->v_size; i++){
		int start = i==0?0:o->v_idx[pos][i-1];
		int end = o->v_idx[pos][i];
		for(j=start; j<end; j++){
			int from = i;
			int to = o->e_idx[pos][j];
			a->list[from][a->idx_list[from]] = to;
			a->idx_list[from]++;
		}
	}
}

/* *******************
 * turn back into csr
 * *******************/
void
adj_to_csr(po o, pa a, int pos){
	int i, j, k;
	//compute v_idx first
	int sum = 0;
	for(i=0; i<o->v_size; i++){
		if(a->idx_list[i] == 0){
			o->v_idx[pos][i] = sum;
			continue;
		}
		for(j=0; j<a->idx_list[i]; j++){
			if(a->list[i][j] != NUL)
				sum++;
		}
		o->v_idx[pos][i] = sum;
	}
	//compute e_idx
	k=0;
	for(i=0;i<o->v_size; i++){
		if(a->idx_list[i] == 0)
			continue;
		for(j=0; j<a->idx_list[i]; j++){
			if(a->list[i][j] != NUL){
				o->e_idx[pos][k] = a->list[i][j];
				k++;
			}
		}
	}
}

/* *************
 * visualization
 * *************/
void
vis_order(po o, int pos_n){
	int i, j, k;
	char f_buf[100];
	sprintf(f_buf, "vis_fork_%d.dot", pos_n);
	FILE *writer = fopen(f_buf, "w");
	fprintf(writer, "Graph{\n");
	for(j=0;j<o->v_size; j++){
		int start = j==0?0:o->v_idx[pos_n][j-1];
		int end = o->v_idx[pos_n][j];
		for(k=start; k<end; k++){
			int from = j;
			int to = o->e_idx[pos_n][k];
				fprintf(writer, "%d -- %d [color=red];\n", from, to);
		}
	}
	fprintf(writer, "}\n");
	fclose(writer);
}

void
vis_order_step(po o, int pos_n, char *file){
	int i, j, k;
	char f_buf[100];
	FILE *writer = fopen(file, "w");
	fprintf(writer, "Graph{\n");
	for(j=0;j<o->v_size; j++){
		int start = j==0?0:o->v_idx[pos_n][j-1];
		int end = o->v_idx[pos_n][j];
		for(k=start; k<end; k++){
			int from = j;
			int to = o->e_idx[pos_n][k];
				fprintf(writer, "%d -- %d [color=red];\n", from, to);
		}
	}
	fprintf(writer, "}\n");
	fclose(writer);
}

void 
vis_fork_adj(pa a, int pos_n){
	int i, j, k;
	char f_buf[100];
	sprintf(f_buf, "vis_adj_%d.dot", pos_n);
	FILE *writer = fopen(f_buf, "w");
	fprintf(writer, "Graph{\n");
	for(i=0; i<a->list_size; i++){
		for(j=0; j<a->idx_list[i]; j++){
			int from = i;
			int to = a->list[i][j];
			fprintf(writer, "%d -- %d [color=red];\n", from, to);
		}
	}
	fprintf(writer, "}\n");
	fclose(writer);
}

void
init_adj(pa a, int size1, int size2){
	a->list = (int**)malloc(sizeof(int*)*size1);
	for(i=0; i<size1; i++)
		a->list[i] = (int*)malloc(sizeof(int)*size2);
	a->idx_list = (int*)malloc(sizeof(int)*siz1);
	a->list_size = size1;
}

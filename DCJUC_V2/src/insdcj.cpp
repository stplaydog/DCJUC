#include "insdcj.h"
#include "utils.h"
#include <fstream>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>

#define MAX_VET_SIZE 32768
#define PI_OPEN 1
#define GAMMA_OPEN 2
#define CYCLE 3
#define END 4
#define NONE 5
#define NUL -1
#define CAP 32768
#define UNSET 65533
#define ede_a   0.595639
#define ede_c   0.457748

void 
set_visited(int *vet_type, int pos);
bool
check_terminate(int *pos, int b_1, int b_2);
int
factorial(int num1, int num2);
bool
check_term(int* one, int *two, int b_1, int b_2);
int
ede ( int invdist, int ngene );

/*
 * this is the constructor
 **/
InsDCJ::InsDCJ(char *file, int tid){
	//pre
	order = (int**)malloc(sizeof(int*)*2);
	renamed_order = (int**)malloc(sizeof(int*)*2);
	adj = (int**)malloc(sizeof(int*)*2);
	dup_idx = (int**)malloc(sizeof(int*)*2);
	dup_pos = (int**)malloc(sizeof(int*)*2);
	//set file reader
	ifstream reader;
	string sline;
	reader.open(file);
	reader>>gf_size>>g_size[0]>>g_size[1];
	int max_size = g_size[0]>g_size[1]?g_size[0]:g_size[1];
	c1=0;
	c2=1;
	size = (g_size[0]>g_size[1]?g_size[0]:g_size[1])*2;
	gene = size/2;
	getline(reader, sline);
	e_size[0] = g_size[0]*2;
	e_size[1] = g_size[1]*2;
	num_dup[0] = num_dup[1] = 0;
#ifdef DEBUG_CONST
	printf("g_size1 %d g_size2 %d\n", g_size[0], g_size[1]);
#endif
	//start reading
	int* exists = (int*)malloc(sizeof(int)*(gf_size+1));
	int *dup_start = (int*)malloc(sizeof(int)*(gf_size+1));
	for(int i=0; i<2; i++){
		//read two lines to get the real gene order
		getline(reader, sline);
		getline(reader, sline);
		string buf;
		stringstream ss(sline);
		vector<string> genes;
		while (ss >> buf)
			genes.push_back(buf);
		//allocation
		order[i] = (int*)malloc(sizeof(int)*g_size[i]);
		renamed_order[i] = (int*)malloc(sizeof(int)*g_size[i]);
		adj[i] = (int*)malloc(sizeof(int)*max_size*4);
		dup_idx[i] = (int*)malloc(sizeof(int)*(gf_size+1));
		dup_pos[i] = (int*)malloc(sizeof(int)*(g_size[i]+1));
		for(int j=0; j<gf_size+1; j++) exists[j] = 0;
		//read genes into order
		for(int j=0; j<g_size[i]; j++){
			int gene_id = atoi(genes[j].c_str());
			order[i][j] = gene_id;
			exists[abs(gene_id)]++;
			if(exists[abs(gene_id)] == 2){
				//printf("dup: %d\n", abs(gene_id));
				num_dup[i]++;
			}
		}
		//assign duplication indexes and positions
		int sum =0;
		for(int j=1; j<gf_size+1; j++){
			dup_start[j]=sum;
			sum += exists[j];
			dup_idx[i][j] = sum;
		}
		for(int j=0; j<g_size[i]; j++){
#ifdef DEBUG_CONST
			//printf(">%d\n", dup_start[abs(order[i][j])]);	
#endif
			dup_pos[i][dup_start[abs(order[i][j])]++] = j;
		}
	}
	//branch related variables
	branch_code = (int*)malloc(sizeof(int)*num_dup[1]*100);
	branch_pos = (int*)malloc(sizeof(int)*num_dup[1]*100);
	branch_idx = (int*)malloc(sizeof(int)*10000);
	num_branch=0;
	branch_selected =0;
	last_dup=-1;
	free(exists);
	//for distance computation
	tmp_second = (int*)malloc(sizeof(int)*g_size[1]);
	vet_type = (int*)malloc(sizeof(int)*(max_size*4));
	//rename orders
	int renamed_id = gf_size+1;
	for(int i=0; i<2; i++)
		for(int j=0; j<g_size[i]; j++)
			renamed_order[i][j] = UNSET;
	//rename the order 0 and the unduplicated part of order 1
	for(int i=1; i<gf_size+1; i++){
		int start = i==0?0:dup_idx[0][i-1]+1;
		int end = dup_idx[0][i];
		for(int j=start; j<end; j++){
			int pos = dup_pos[0][j];
			renamed_order[0][pos] = renamed_id++;
		}
	}
	for(int i=0; i<g_size[0]; i++){
		if(renamed_order[0][i] == UNSET){
			renamed_order[0][i] = order[0][i];
		}
	}
	rnm_id = renamed_id;
	old_rnm = renamed_id;
	//rename order 1
	for(int i=1; i<gf_size+1; i++){
		int start = i==0?0:dup_idx[1][i-1];
		int end = dup_idx[1][i];
		if(end-start>1)
			for(int j=start; j<end; j++){
				renamed_order[1][dup_pos[1][j]] = UNSET+1;
			}
	}
	for(int i=0; i<g_size[1]; i++)
		if(renamed_order[1][i] == UNSET)
			renamed_order[1][i] = order[1][i];
		else
			renamed_order[1][i] = UNSET;
#ifdef DEBUG_CONST
	for(int i=0; i<2; i++){
		for(int j=0; j<g_size[i]; j++){
			printf("%d:%d ", j, renamed_order[i][j]);
		}
		printf("\n");
		printf("\n");
	}
	//compute bound
	for(int k=0; k<2; k++){
		for(int i=1; i<gf_size+1; i++){
			int start = i==0?0:dup_idx[k][i-1];
			int end = dup_idx[k][i];
			if((end-start)>1){
				for(int j=start; j<end; j++)
					printf("%d:%d ", i, dup_pos[k][j]);
			}
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");
	printf("\n");
#endif
	compute_bound();
#ifdef DEBUG_CONST
	//print to check if correct
	for(int i=0; i<2; i++){
		for(int j=0; j<g_size[i]; j++){
			printf("%d:%d ", j, order[i][j]);
		}
		printf("\n");
	}
	printf("num_dup1 %d num_dup2 %d\n", num_dup[0], num_dup[1]);
	printf("g_size1 %d g_size2 %d gf_size %d \n", g_size[0], g_size[1], gf_size);
#endif
	compute_num_elem();
	//exit(1);
}


InsDCJ::InsDCJ(int **ord, int one, int two, 
		int size_one, int size_two,
		int gf_size1){
	//printf("one %d two %d size_one %d size_two %d\n", one, two, size_one, size_two);
	//pre
	c1=0;c2=1;
	num_dup[0] = num_dup[1] = 0;
	gf_size = gf_size1;
	order = (int**)malloc(sizeof(int*)*2);
	renamed_order = (int**)malloc(sizeof(int*)*2);
	adj = (int**)malloc(sizeof(int*)*2);
	dup_idx = (int**)malloc(sizeof(int*)*2);
	dup_pos = (int**)malloc(sizeof(int*)*2);
	int max_size = size_one>size_two?size_one:size_two;
	//set file reader
	//start reading
	int* exists = (int*)malloc(sizeof(int)*(gf_size+1));
	int *dup_start = (int*)malloc(sizeof(int)*(gf_size+1));
	for(int i=0; i<2; i++){
		//read two lines to get the real gene order
		for(int j=0; j<gf_size+1; j++) exists[j] = 0;
		int c=i==0?one:two;
		int size = i==0?size_one:size_two;
		g_size[i] = size; 
		int *genes = ord[c];
		//allocation
		order[i] = (int*)malloc(sizeof(int)*g_size[i]);
		renamed_order[i] = (int*)malloc(sizeof(int)*g_size[i]);
		adj[i] = (int*)malloc(sizeof(int)*max_size*4);
		dup_idx[i] = (int*)malloc(sizeof(int)*(gf_size+1));
		dup_pos[i] = (int*)malloc(sizeof(int)*(g_size[i]+1));
		for(int j=0; j<gf_size+1; j++) exists[j] = 0;
		//read genes into order
		for(int j=0; j<g_size[i]; j++){
			int gene_id = genes[j];
			order[i][j] = gene_id;
			exists[abs(gene_id)]++;
			if(exists[abs(gene_id)] == 2){
				num_dup[i]++;
			}
		}
		//assign duplication indexes and positions
		int sum =0;
		for(int j=1; j<gf_size+1; j++){
			dup_start[j]=sum;
			sum += exists[j];
			dup_idx[i][j] = sum;
		}
		for(int j=0; j<g_size[i]; j++){
			dup_pos[i][dup_start[abs(order[i][j])]++] = j;
		}
	}
	free(exists);
	free(dup_start);
	//branch related variables
	if(num_dup[1]!=0){
		branch_code = (int*)malloc(sizeof(int)*num_dup[1]*100);
		branch_pos = (int*)malloc(sizeof(int)*num_dup[1]*100);
		branch_idx = (int*)malloc(sizeof(int)*10000);
	}
	num_branch=0;
	branch_selected =0;
	last_dup=-1;
	//for distance computation
	tmp_second = (int*)malloc(sizeof(int)*(g_size[1]+1));
	vet_type = (int*)malloc(sizeof(int)*(max_size*4));
	//rename orders
	int renamed_id = gf_size+1;
	for(int i=0; i<2; i++)
		for(int j=0; j<g_size[i]; j++)
			renamed_order[i][j] = UNSET;
	//rename the order 0 and the unduplicated part of order 1
	dup_idx[0][0] =0; dup_idx[1][0]=0;
	for(int i=1; i<gf_size+1; i++){
		int start = i==0?0:dup_idx[0][i-1]+1;
		int end = dup_idx[0][i];
		for(int j=start; j<end; j++){
			int pos = dup_pos[0][j];
			renamed_order[0][pos] = renamed_id++;
		}
	}
	for(int i=0; i<g_size[0]; i++){
		if(renamed_order[0][i] == UNSET){
			renamed_order[0][i] = order[0][i];
		}
	}
	rnm_id = renamed_id;
	old_rnm = renamed_id;
	//rename order 1
	for(int i=1; i<gf_size+1; i++){
		int start = i==0?0:dup_idx[1][i-1];
		int end = dup_idx[1][i];
		if(end-start>1)
			for(int j=start; j<end; j++){
				//if(j<0){
				//	printf("start %d end %d\n", start, end);
				//	for(int k=0; k<g_size[1]; k++)
				//		printf("%d ", order[1][k]);
				//	printf("\n");
				//	for(int k=1; k<gf_size+1; k++)
				//		printf("%d ", dup_idx[1][k]);
				//	printf("\n");
				//}
				renamed_order[1][dup_pos[1][j]] = UNSET+1;
			}
	}
	for(int i=0; i<g_size[1]; i++)
		if(renamed_order[1][i] == UNSET)
			renamed_order[1][i] = order[1][i];
		else
			renamed_order[1][i] = UNSET;
#ifdef DEBUG_CONST
	for(int i=0; i<2; i++){
		for(int j=0; j<g_size[i]; j++){
			printf("%d:%d ", j, renamed_order[i][j]);
		}
		printf("\n");
		printf("\n");
	}
	//compute bound
	for(int k=0; k<2; k++){
		for(int i=1; i<gf_size+1; i++){
			int start = i==0?0:dup_idx[k][i-1];
			int end = dup_idx[k][i];
			if((end-start)>1){
				for(int j=start; j<end; j++)
					printf("%d:%d ", i, dup_pos[k][j]);
			}
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");
	printf("\n");
#endif
	compute_bound();
#ifdef DEBUG_CONST
	//print to check if correct
	for(int i=0; i<2; i++){
		for(int j=0; j<g_size[i]; j++){
			printf("%d:%d ", j, order[i][j]);
		}
		printf("\n");
	}
	printf("num_dup1 %d num_dup2 %d\n", num_dup[0], num_dup[1]);
	printf("g_size1 %d g_size2 %d gf_size %d \n", g_size[0], g_size[1], gf_size);
#endif
	compute_num_elem();
	//exit(1);
}

/*
 * this is the copy constructor
 **/
InsDCJ::InsDCJ(const InsDCJ &other){
}

InsDCJ::~InsDCJ(){
	for(int i=0; i<2; i++){
		free(order[i]);
		free(renamed_order[i]);
		free(dup_pos[i]);
		free(dup_idx[i]);
		free(adj[i]);
	}
	if(num_dup[1]>0){	
		free(branch_code); //use the first genome as the standard
		free(branch_idx);
		free(branch_pos);
	}
	free(tmp_second);
	free(vet_type);
	free(order);
	free(renamed_order);
	free(dup_pos);
	free(dup_idx);
	free(adj);
}

/*
 * let the instance transform into a specific branch
 **/
void InsDCJ::to_branch(int which_branch){
	//branch_code[0] is which gene, branch_code[1] is which gene it is renamed to.
	int start = which_branch==0?0:branch_idx[which_branch-1];
	int end = branch_idx[which_branch];
	//printf("start %d end %d\n", start, end);
	int renamed_count = 0;
	for(int i=start; i<end; i++){
		int pos = branch_pos[i];
		int code = branch_code[i];
		renamed_order[1][pos] = code;
		if(abs(code) >= old_rnm)
			renamed_count++;
	}
	rnm_id += renamed_count; //the reason why adding this is for the following bound computation
	branch_selected = which_branch;
	compute_bound();
}

/*
 * recover from a branch
 **/
void InsDCJ::from_branch(){
	int start = branch_selected==0?0:branch_idx[branch_selected-1];
	int end = branch_idx[branch_selected];
	int renamed_count = 0;
	for(int i=start; i<end; i++){
		int pos = branch_pos[i];
		int code = branch_code[i];
		renamed_order[1][pos] = UNSET;
		if(abs(code) >= old_rnm)
			renamed_count++;
	}
	rnm_id -= renamed_count; //the reason why adding this is for the following bound computation
	branch_selected = -1;
}


int  
InsDCJ::get_num_branches(){
	//select all possible combinations in a gene family as branches
	int i;
	//find the first dup position
	for(i=1; i<gf_size+1; i++){
		//int start_1 = i==0?0:dup_idx[0][i-1];
		//int end_1 = dup_idx[0][i];
		int start_2 = i==0?0:dup_idx[1][i-1];
		int end_2 = dup_idx[1][i];
		if(end_2-start_2<2 || 
				renamed_order[1][dup_pos[1][start_2]] != UNSET)
			continue;
		last_dup = i;
		break;
	}
	//printf("last dup %d\n", last_dup);
	//print_encode();
	//get the number of occurrence of duplications for each genome
	int start_1 = i==0?0:dup_idx[0][last_dup-1];
	int end_1 = dup_idx[0][last_dup];
	int start_2 = i==0?0:dup_idx[1][last_dup-1];
	int end_2 = dup_idx[1][last_dup];
	int b_1 = (end_1-start_1);
	int b_2 = (end_2-start_2);
	b_1 = b_1==0?1:b_1;
	b_2 = b_2==0?1:b_2;
	num_branch=b_1>b_2?factorial(b_1, b_2):factorial(b_2, b_1);
	branch_offset = b_2*2;
	//for(int i=0; i<b_1; i++)
	//	num_branch = num_branch * (b_2 - i);
	int branch_id=0;
#ifdef DEBUG
	//printf("start_1 %d end_1 %d start_2 %d end_2 %d b_1 %d b_2 %d num_branch %d\n", start_1, end_1, start_2, end_2, b_1, b_2, num_branch);
	//for(int i=start_2; i<end_2; i++){
	//	printf("%d ", order[1][dup_pos[1][i]]);
	//}
	//printf("\n");
#endif
	//cope with situation when (end_1-start_1)=0
	//if((end_1-start_1)==0 && reverse_idx[last_dup]==NUL){
	if((end_1-start_1)==0){
		int start_id = rnm_id;
		num_branch = 1;
		branch_idx[0] = b_2;
		branch_code[0] = order[1][dup_pos[1][start_2]];
		for(int i=1, j=start_2+1; i<b_2; i++, j++){
			branch_code[i] = start_id++;
			branch_pos[i] = dup_pos[1][j];
		}
#ifdef DEBUG_GNUM
		for(int i=0; i<num_branch; i++){
			//printf("branch: %d: ", i);
			int start = i==0?0:branch_idx[i-1];
			int end = branch_idx[i];
			//for(int j=start; j<end; j++)
			//	printf("%d ", branch_code[j]);
			//printf("\n");
		}
#endif
		return num_branch;
	}

#ifdef USE_BFS
	//int **q = (int**)malloc(sizeof(int*)*2);
	//q[0] = (int*)malloc(sizeof(int*)*10000);
	//q[1] = (int*)malloc(sizeof(int*)*10000);
	int q[2][10000];
	int q_num = 1;
	//put the first null element
	for(int i=0; i<10000; i++) q[0][i] = q[1][i] = NUL;
	while(q_num>0){
		//get one from the queue
		int one[b_1];
		int two[b_2];
		int s_1 = (q_num-1)*b_1;
		int s_2 = (q_num-1)*b_2;
		for(int i=0; i<b_1; i++)
			one[i] = q[0][s_1+i];
		for(int i=0; i<b_2; i++)
			two[i] = q[1][s_2+i];
		q_num--;
		//check if terminate
		if(check_term(one, two, b_1, b_2)==true){
			int start_id = rnm_id;
			for(int i=0; i<b_2; i++){
				branch_pos[branch_id*b_2 + i] = dup_pos[1][start_2+i];
				if(two[i] == NUL){
					if(order[1][dup_pos[1][start_2+i]]<0) branch_code[branch_id*b_2 + i] = start_id * (-1);
					else branch_code[branch_id*b_2 + i] = start_id;
					start_id++;
				}
				else{
					int one_pos = two[i];
					if(order[1][dup_pos[1][start_2+i]]<0) branch_code[branch_id*b_2 + i] = renamed_order[0][dup_pos[0][start_1+one_pos]]*-1;
					else branch_code[branch_id*b_2 + i] = renamed_order[0][dup_pos[0][start_1+one_pos]];
				}
			}
			branch_id++;
			branch_idx[branch_id-1] = branch_id * b_2;
			continue;
		}
		if(b_1<=b_2){
			for(int i=0; i<b_1; i++){
				if(one[i] == NUL){
					for(int j=0; j<b_2; j++){
						if(two[j] == NUL){
							one[i] = j;
							two[j] = i;
							//insert
							for(int k=0; k<b_1; k++)
								q[0][q_num*b_1+k] = one[k];
							for(int k=0; k<b_2; k++)
								q[1][q_num*b_2+k] = two[k];
							q_num++;
							//printf("%d \n", q_num);
							two[j] = NUL;
						}
					}
					break;
				}
			}
		}
		else{
			for(int i=0; i<b_2; i++){
				if(two[i] == NUL){
					for(int j=0; j<b_1; j++){
						if(one[j] == NUL){
							two[i] = j;
							one[j] = i;
							//insert
							for(int k=0; k<b_1; k++)
								q[0][q_num*b_1+k] = one[k];
							for(int k=0; k<b_2; k++)
								q[1][q_num*b_2+k] = two[k];
							q_num++;
							one[j] = NUL;
						}
					}
					break;
				}
			}
		}
	}
#ifdef DEBUG
	//for(int i=0; i<num_branch; i++){
	//	printf("branch: %d: ", i);
	//	int start = i==0?0:branch_idx[i-1];
	//	int end = branch_idx[i];
	//	for(int j=start; j<end; j++)
	//		printf("%d->%d ", branch_pos[j], branch_code[j]);
	//	printf("\n");
	//}
#endif
	//free(q[0]);
	//free(q[1]);
	//free(q);
#endif


#ifdef USE_DFS
	//using recursive dfs to get the branch code
	int **stack = (int**)malloc(sizeof(int*)*3);
	stack[0] = (int*)malloc(sizeof(int)*1000); //stack[0] store the position of second genome
	stack[1] = (int*)malloc(sizeof(int)*1000); //stack[1] store the according mapping
	stack[2] = (int*)malloc(sizeof(int)*1000); //stack[2] store the position of first genome
	int stack_idx = 0;
	//select by is the indicator for which position in genome 2 has been reached
	int *select_by = (int*)malloc(sizeof(int)*b_1);
	for(int i=0; i<b_1; i++)
		select_by[i] = -1;
	//chosen is indicating which copy of the gene family in genome 1 is chosen
	bool *chosen = (bool*)malloc(sizeof(bool)*b_1);
	bool *chosen_g2 = (bool*)malloc(sizeof(bool)*b_2);
	for(int i=0; i<b_1; i++) chosen[i] = false;
	for(int i=0; i<b_2; i++) chosen_g2[i] = false;
	//push into the stack
	while(true){
		//if there is nothing to map do the rest
		if(stack_idx==b_1 || stack_idx==b_2){
			int start_id = rnm_id;
			int start = branch_id==0?0:branch_idx[branch_id-1];
			for(int i=0; i<b_2; i++)
				branch_code[start+i] = NUL;
			if(stack_idx==b_1 && b_1<=b_2){
				//copy from the stack
				for(int i=0; i<b_1; i++){
					branch_pos[start+i] = stack[0][i];
					branch_code[start+i] = stack[1][i]; 
				}
				int b_start = start + b_1;
				//copy the rest by rename
				for(int i=0; i<b_2; i++){
					if(chosen_g2[i] == false){
						int idx = (last_dup==0?0:dup_idx[1][last_dup-1])+i;
						branch_pos[b_start] = dup_pos[1][idx]; //
						branch_code[b_start] = order[1][branch_pos[b_start]]>0?start_id:start_id*(-1); //need to consider sign
						b_start++;
						start_id++;
					}
				}
			}
			else{
				for(int i=0; i<b_2; i++){
					//stack[0] store the position
					//stack[1] store the rename
					branch_pos[start+i] = stack[0][i];
					branch_code[start+i] = stack[1][i];
				}
			}
			branch_idx[branch_id] = start+b_2;
			branch_id+=1;
			//pop an element and increase the count
			chosen[stack[2][stack_idx-1]]=false;
			chosen_g2[select_by[stack_idx-1]] = false;
			stack_idx--;
			if(branch_id>=num_branch) break;
		}
		//move to the next combination
		if(select_by[stack_idx]==(b_2-1)){ //for the case when already reach the limit
			select_by[stack_idx] = -1;
			//change chosen information
			chosen_g2[select_by[stack_idx-1]] = false;
			chosen[stack[2][stack_idx-1]] = false;
			stack_idx--;
			continue;
		}
		int g1_pos, g2_pos, mapping;
		for(int i=0; i<b_1; i++){
			if(chosen[i]==false){
				g1_pos = i;
				mapping = renamed_order[0][dup_pos[0][start_1+i]];
				chosen[i] = true;
				for(int j=select_by[stack_idx]+1; j<b_2; j++){
					if(chosen_g2[j] == false){
						g2_pos = dup_pos[1][start_2+j];
						select_by[stack_idx] = j;
						chosen_g2[j] = true;
						break;
					}
				}
				break;
			}
		}
		stack[0][stack_idx] = g2_pos;
		stack[1][stack_idx] = order[1][g2_pos]>0?abs(mapping):abs(mapping)*(-1); //need to consider the sign
		stack[2][stack_idx] = g1_pos;
#ifdef DEBUG 
		//for(int i=0; i<stack_idx; i++)
		//	printf("%d:%d:%d >>>", stack[0][stack_idx], stack[1][stack_idx], stack[2][stack_idx]);
		//printf("\n");
#endif
		stack_idx++;
	}
	free(stack[0]);
	free(stack[1]);
	free(stack[2]);
	free(stack);
#endif
	return num_branch;
}

/*
 * shows which gene is renamed to which
 ***/
int  
InsDCJ::get_encode(int* encode){
	int encode_idx =0;
	for(int i=1; i<gf_size+1; i++){
		int start = i==1?0:dup_idx[1][i-1];
		int end = dup_idx[1][i];
		if(end-start>=2){
			for(int j=start; j<end; j++){
				int pos = dup_pos[1][j];
				if(renamed_order[1][pos] != UNSET) { //rename happend
					encode[encode_idx++] = pos;
					encode[encode_idx++] = renamed_order[1][pos];
				}
			}
		}
	}
	return encode_idx;
}


void 
InsDCJ::to_encode(int* encode, int size){
	//rename order 1
	for(int i=0; i<g_size[1]; i++)
		renamed_order[1][i] = UNSET;
	for(int i=1; i<gf_size+1; i++){
		int start = i==0?0:dup_idx[1][i-1];
		int end = dup_idx[1][i];
		if(end-start>1)
			for(int j=start; j<end; j++){
				renamed_order[1][dup_pos[1][j]] = UNSET+1;
			}
	}
	for(int i=0; i<g_size[1]; i++)
		if(renamed_order[1][i] == UNSET)
			renamed_order[1][i] = order[1][i];
		else
			renamed_order[1][i] = UNSET;
	rnm_id = old_rnm; //remember to change the encode
	int num_encode = size/2;
	for(int i=0; i<num_encode; i++){
		int pos = encode[i*2];
		int code = encode[i*2+1];
		renamed_order[1][pos] = code;
		//printf("%d:%d ", pos, code);
		if(abs(code)>=old_rnm){
			//printf("rnm_id %d\n", rnm_id);
			rnm_id++;
		}
	}
}


void InsDCJ::print_encode(){
	//int encode[100];
	//int num = get_encode(encode);
	//for(int i=0; i<num; i++)
	//	printf("%d ", encode[i]);
	//printf("\n");
	//fflush(stdout);
	//for(int i=0; i<g_size[1]; i++)
	//	if(order[1][i] != renamed_order[1][i])
	//		printf("%d~%d ", order[1][i], renamed_order[1][i]);
	//printf("\n");
	//for(int i=0; i<num_dup[1]; i++){
	//	int idx = dup_idx[1][i];
	//	if(order[1][idx] != renamed_order[1][idx]) { //rename happend
	//		printf("%d~%d ", idx, renamed_order[1][idx]);
	//	}
	//}
	//for(int i=0; i<2; i++){
	//	for(int j=0; j<g_size[i]; j++){
	//		printf("%d ", renamed_order[i][j]);
	//	}
	//	printf("\n");
	//}
	//printf("%d distance is %d\n", upper_bound, ede (upper_bound, 200));
}


void InsDCJ::copy_code(int *encode, int size){
	int j=0;
	int k;
	to_encode(encode, size);
	//
	int new_name = rnm_id;
	//printf("rnm_id %d\n", rnm_id);
	//rename the dup first
	for(int i=0; i<g_size[1]; i++)
		tmp_second[i] = renamed_order[1][i];
	//for(int i=0; i<g_size[1]; i++)
	//	printf("%d ", order[1][i]);
	//printf("\n");
	//for(int i=0; i<g_size[1]; i++)
	//	if(tmp_second[i]==UNSET)
	//	printf("%d:%d ", i, tmp_second[i]);
	//printf("\n");
	//printf("\n");
	for(int i=1; i<gf_size+1; i++){
		int start_1 = i==0?0:dup_idx[0][i-1];
		int end_1 = dup_idx[0][i];
		int start_2 = i==0?0:dup_idx[1][i-1];
		int end_2 = dup_idx[1][i];
		if((end_2-start_2)>1){ //there is unset renames
			if((end_2-start_2)<=(end_1-start_1)){
				for(j=start_1, k=start_2; j<end_1; j++,k++){
					int pos_1 = dup_pos[0][j];
					int pos_2 = dup_pos[1][k];
					if(tmp_second[pos_2] == UNSET){
						//check sign
						if(order[1][pos_2]>0)
							tmp_second[pos_2] = abs(renamed_order[0][pos_1]);
						else
							tmp_second[pos_2] = abs(renamed_order[0][pos_1])*(-1);
					}
				}
			}
			else{
				for(j=start_1, k=start_2; j<end_1; j++,k++){
					int pos_1 = dup_pos[0][j];
					int pos_2 = dup_pos[1][k];
					if(tmp_second[pos_2] == UNSET){
						//check sign
						if(order[1][pos_2]>0)
							tmp_second[pos_2] = abs(renamed_order[0][pos_1]);
						else
							tmp_second[pos_2] = abs(renamed_order[0][pos_1])*(-1);
					}
				}
				if((end_1-start_1)==0)
					tmp_second[dup_pos[1][start_2]] = order[1][dup_pos[1][start_2]];
				k=((end_1-start_1)==0?start_2+1:k);
				for(;k<end_2; k++){
					int pos_2 = dup_pos[1][k];
					if(tmp_second[pos_2] == UNSET){
						tmp_second[pos_2] = new_name++;
					}
				}
			}
			
		}
	}
	for(int i=0; i<g_size[1]; i++)
		renamed_order[1][i] = tmp_second[i];
}


int  InsDCJ::get_value(){
	return upper_bound;
}


int  
InsDCJ::compute_indel_dis(){ 
	//vis_adj("vis_adj.dot", adj, size);
	int left = 0;
	int right = 0;

	int num_c=0;
	int num_even=0;
	int num_pi_pi=0;
	int num_gamma_gamma=0;
	int num_pi_gamma=0;
	//int num_reg_even=0;
	int num_pi_odd=0;
	int num_pi_even=0;
	int num_gamma_odd=0;
	int num_gamma_even=0;
	int delta=0;

	int i;
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
	int tmp_c1 = 0;
	int tmp_c2 = 0;
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
				//printf("%d \n", right);
				if(right==start){
					num_c++;
					vet_type[right]=NONE;
					set_visited(vet_type, right);
					break;
				}
				vet_type[right]=NONE;

				left=adj[c2][right];
				//printf("%d \n", left);
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
	return result;	
}


void 
InsDCJ::compute_bound(){
	//randomly rename to compute the upper bound
	int j=0;
	int k;
	//the reason why new_name is equal to g_size[0] is that 
	//duplications in order1 and order2 is to neccessarily having the same
	//number of genes for a specific gene family
	//questionable 
	int new_name = rnm_id;
	//printf("rnm_id %d\n", rnm_id);
	//rename the dup first
	for(int i=0; i<g_size[1]; i++)
		tmp_second[i] = renamed_order[1][i];
	//for(int i=0; i<g_size[1]; i++)
	//	printf("%d ", order[1][i]);
	//printf("\n");
	//for(int i=0; i<g_size[1]; i++)
	//	if(tmp_second[i]==UNSET)
	//	printf("%d:%d ", i, tmp_second[i]);
	//printf("\n");
	//printf("\n");
	if(num_dup[1]>0){
		for(int i=1; i<gf_size+1; i++){
			int start_1 = i==0?0:dup_idx[0][i-1];
			int end_1 = dup_idx[0][i];
			int start_2 = i==0?0:dup_idx[1][i-1];
			int end_2 = dup_idx[1][i];
			if((end_2-start_2)>1){ //there is unset renames
				if((end_2-start_2)<=(end_1-start_1)){
					for(j=start_1, k=start_2; j<end_1; j++,k++){
						int pos_1 = dup_pos[0][j];
						int pos_2 = dup_pos[1][k];
						if(tmp_second[pos_2] == UNSET){
							//check sign
							if(order[1][pos_2]>0)
								tmp_second[pos_2] = abs(renamed_order[0][pos_1]);
							else
								tmp_second[pos_2] = abs(renamed_order[0][pos_1])*(-1);
						}
					}
				}
				else{
					for(j=start_1, k=start_2; j<end_1; j++,k++){
						int pos_1 = dup_pos[0][j];
						int pos_2 = dup_pos[1][k];
						if(tmp_second[pos_2] == UNSET){
							//check sign
							if(order[1][pos_2]>0)
								tmp_second[pos_2] = abs(renamed_order[0][pos_1]);
							else
								tmp_second[pos_2] = abs(renamed_order[0][pos_1])*(-1);
						}
					}
					if((end_1-start_1)==0)
						tmp_second[dup_pos[1][start_2]] = order[1][dup_pos[1][start_2]];
					k=((end_1-start_1)==0?start_2+1:k);
					for(;k<end_2; k++){
						int pos_2 = dup_pos[1][k];
						if(tmp_second[pos_2] == UNSET){
							tmp_second[pos_2] = new_name++;
						}
					}
				}
			}
		}
	}
	//for(int i=0; i<g_size[1]; i++)
	//	printf("%d:%d ", order[1][i], tmp_second[i]);
	int tmp_size = size;
	size = (new_name-1)*2;
	from_order_to_adj();
	//printf("come here! g_size[1] %d\n", g_size[1]);
	//if(last_dup==22)
	//printf("old_rnm %d rnm_id %d new name %d, size %d\n", old_rnm, rnm_id, new_name, size);
#ifdef DEBUG_BOUND
	//print to check
	//for(int i=0; i<2; i++){
	//	for(int j=0; j<size; j++){
	//		printf("%d->%d ", j, adj[i][j]);
	//	}
	//	printf("\n");
	//}
	//printf("\n");
	//vis_graph();
#endif
	upper_bound = size/2 - compute_indel_dis();
	size = tmp_size;
	//printf("size c1 %d size %d rest %d\n", g_size[0], size, size/2 - upper_bound);
	//print_encode();
	//print_adj();
	////////////////////////////////////////////////
	//
	//delete all the copy to compute the lower bound
	j=0;
	for(int i=0; i<g_size[1]; i++)
		if(renamed_order[1][i] != UNSET)
			tmp_second[j++] = renamed_order[1][i];
	//for (int i=0; i<j; i++)
	//	printf("%d ", tmp_second[i]);
	//printf("\n");
	//printf("size of second genome %d\n", j);
	tmp_size = size;
	size = (g_size[0]>rnm_id?g_size[0]:rnm_id)*2; //change the size
	int tmp_g_size = g_size[1];
	g_size[1] = j;
	//compute
	from_order_to_adj();
	lower_bound = size/2 - compute_indel_dis();
	//resume
	size = tmp_size; //resume the size
	g_size[1] = tmp_g_size;
	score = ede(lower_bound, 50);
#ifdef DEBUG_BOUND
	//vis_graph();
	//printf("upper_bound %d lower_bound %d\n", upper_bound, lower_bound);
#endif
}

/*
 * turning orders into adjacencies
 * **/
void
InsDCJ::from_order_to_adj(){
	for(int i=0; i<2; i++)
		for(int j=0; j<size; j++)
			adj[i][j] = NUL;
	//cope with first genome
	int head, pre_tail;
	if(renamed_order[c1][0]>0)
		head = (abs(renamed_order[c1][0])-1)*2; 
	else
		head = (abs(renamed_order[c1][0])-1)*2+1;
	if(renamed_order[c1][g_size[c1]-1]>0)
		pre_tail = (abs(renamed_order[c1][g_size[c1]-1])-1)*2+1; 
	else
		pre_tail = (abs(renamed_order[c1][g_size[c1]-1])-1)*2;
	adj[c1][pre_tail] = head;
	adj[c1][head] = pre_tail;
	for(int i=1; i<g_size[c1]; i++){
		if(renamed_order[c1][i]>0)
			head = (abs(renamed_order[c1][i])-1)*2;
		else
			head = (abs(renamed_order[c1][i])-1)*2+1;
		if(renamed_order[c1][i-1]>0)
			pre_tail = (abs(renamed_order[c1][i-1])-1)*2+1;
		else
			pre_tail = (abs(renamed_order[c1][i-1])-1)*2;
		//printf("order %d preorder %d head %d tail %d\n", renamed_order[c1][i], renamed_order[c1][i-1], head, pre_tail);
		adj[c1][pre_tail] = head;
		adj[c1][head] = pre_tail;
	}
	//cope with second genome
	if(tmp_second[0]>0)
		head = (abs(tmp_second[0])-1)*2; 
	else
		head = (abs(tmp_second[0])-1)*2+1;
	if(tmp_second[g_size[c2]-1]>0)
		pre_tail = (abs(tmp_second[g_size[c2]-1])-1)*2+1; 
	else
		pre_tail = (abs(tmp_second[g_size[c2]-1])-1)*2;
	adj[c2][pre_tail] = head;
	adj[c2][head] = pre_tail;
	for(int i=1; i<g_size[c2]; i++){
		if(tmp_second[i]>0)
			head = (abs(tmp_second[i])-1)*2;
		else
			head = (abs(tmp_second[i])-1)*2+1;
		if(tmp_second[i-1]>0)
			pre_tail = (abs(tmp_second[i-1])-1)*2+1;
		else
			pre_tail = (abs(tmp_second[i-1])-1)*2;
#ifdef DEBUG_ADJ
		printf("pre %d: now %d\n", tmp_second[i-1], tmp_second[i]);
		printf("head %d pre_tail %d\n", head, pre_tail);
#endif
		adj[c2][pre_tail] = head;
		adj[c2][head] = pre_tail;
	}
	//set vet_type
	for(int i=0;i<size; i++)
	{
		if(adj[1][i]!=NUL && adj[1][i]!=CAP && adj[0][i]==NUL)
			vet_type[i]=PI_OPEN;
		else if(adj[0][i]!=NUL && adj[0][i]!=CAP && adj[1][i]==NUL)
			vet_type[i]=GAMMA_OPEN;
		else if(adj[0][i]!=NUL && adj[1][i]!=NUL &&
				adj[0][i]!=CAP && adj[1][i]!=CAP)
			vet_type[i]=CYCLE;
		else if((adj[0][i]!=NUL && adj[0][i]==CAP) ||
				(adj[1][i]!=NUL && adj[1][i]==CAP))
			vet_type[i]=END;
		else
			vet_type[i]=NONE;
	}
}

bool
check_terminate(int *pos, int b_1, int b_2){
	int end = b_1>b_2?b_2:b_1;
	for(int i=0; i<end; i++){
		//if there is still one roon to choose
		//continue
		if(pos[i]<(b_2-1))
			return false;
	}
	return true;
}

void 
set_visited(int *vet_type, int pos){
	vet_type[pos] = NONE;
	//printf("%d visisted\n", pos);
}

void
InsDCJ::vis_graph(){
	FILE *writer = fopen("graph.dot", "w");
	fprintf(writer, "graph{\n");
	for(int i=0; i<2; i++){
		const char *cstr = i==0?"red":"blue";
		for(int j=0; j<size; j++){
			int from = j;
			int to = adj[i][from];
			if(to != NUL)
				fprintf(writer, "%d -- %d [color=%s];\n", from, to, cstr);
		}
	}
	fprintf(writer, "}\n");
	fclose(writer);
}

int
factorial(int num1, int num2){
	int result = 1;
	for(int i=num1, j=0; j<num2; i--, j++)
		result = result * i;
	//for(int i=num2; i>1; i--)
	//	result = result * i;
	return result;
}

bool
check_term(int* one, int *two, int b_1, int b_2){
	int num_1=0;
	int num_2=0;
	for(int i=0; i<b_1; i++)
		if(one[i] == NUL)
			break;
		else
			num_1++;
	for(int i=0; i<b_2; i++)
		if(two[i] == NUL)
			break;
		else
			num_2++;
	
	if(num_1==b_1 || num_2==b_2)
		return true;	
	else
		return false;
}

void
InsDCJ::print_adj(){
	for(int c=0; c<2; c++){
		for(int i=0; i<size; i++){
			printf("%d ", adj[c][i]);
		}
		printf("\n");
	}
	for(int i=0; i<size; i++)
		printf("%d ", vet_type[i]);
	printf("\n");
}


void 
InsDCJ::compute_num_elem(){
	num_count = 0;
	for(int i=1; i<gf_size; i++){
		int start_1 = i==1?0:dup_idx[0][i-1];
		int start_2 = i==1?0:dup_idx[1][i-1];
		int end_1 = dup_idx[0][i];
		int end_2 = dup_idx[1][i];
		int b_1 = end_1 -start_1;
		int b_2 = end_2 - start_2;
		int fac = b_1>b_2?factorial(b_1, b_2):factorial(b_2, b_1);	
		if(fac > num_count)
			num_count = fac;
	}
}

int
ede ( int invdist, int ngene )
{
    double ll, tt, kk, pp, dval;
    int newvalue;

    kk = invdist / ( ngene + 0.0 );

    if ( kk >= 0.999999999999 )
    {                           /* the distance correction has singularity at 1 */
        kk = 0.999999999999;
    }
    if ( kk <= 1 - ede_c )
        return invdist;

    ll = ede_c * kk - ede_a;
    tt = 4 * ede_a * ( 1 - kk ) * kk + ll * ll;
    tt = ll + sqrt ( tt );
    pp = tt / ( 2 * ( 1 - kk ) );
    pp *= ngene;

    dval = pp;
    newvalue = ( int ) ceil ( dval );
    /*if (newvalue-dval > 0) return newvalue-1; */
    return newvalue;
}


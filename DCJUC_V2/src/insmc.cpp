#include "insmc.h"
#include "utils.h"
#include <string>
#include <iostream>
#include <cstring>
#include <string>

using namespace std;

#define DEG_VAL 0
#define DEG_IDX 1

/***************************************************
 * these are public functions
 * for branching operations
 * *************************************************/

/**
 * constructor and destructors
 * **/
InsMC::InsMC(char *f) : InsMC(f, "", ""){
}

/**
 * constructor and destructors
 * **/
InsMC::InsMC(char *f, char *of) : InsMC(f, of, ""){
}

/**
 * constructor and destructors
 * **/
InsMC::InsMC(char *f, char *of, 
		char *od){
	//***some basic info
	file = (char*)malloc(sizeof(char)*500);
	file = f;
	out_file = (char*)malloc(sizeof(char)*500);
	out_file = of;
	writer = fopen(out_file, "w");
	CC_folder = od;
	upper_bound = NUL;
	lower_bound = NUL;
	branch_id = NUL;
	num_branch = NUL;

	//***init graph and init some related data structures
	init_graph();
	printf("finished reading graph!\n");
	//here if we have the folder information just perform the CC algorithm
	if(CC_folder[0] != '\0'){
		separate_CC();
		get_num_triangles();
		return;
	}

	//***get initial R, and perform the color assignment function
	get_full_R(); 
	get_deg_sorted();
	assign_C();

	//compute upper and lower bound now
	upper_bound = max_C;
	lower_bound = 2; //this is important, smaller than this will not be considered
	printf("initial upper bound %d lower bound %d\n", upper_bound, lower_bound);
	get_num_elem();
}

/**
 * destructor to free all the things
 * **/
InsMC::~InsMC(){
	free(file);
	free(out_file);
	free(v_idx);
	free(e_idx);
	free(Q);
	free(R);
	free(tmp_R);
	free(R_deg);
	free(R_avail);
	free(C);
	free(C_type);
	free(C_map);
	fclose(writer);
}

/**
 * one thing to notice is, one branch can generate multiple Qs
 * we compute Rp := R \intersect p
 * **/
void 
InsMC::to_branch(int which_branch){
	//do nothing here
	branch_id = which_branch;
	Q[Q_size++] = branch_v[branch_id];
	//then update upper and lower bound
	compute_bound();
}

/**
 * resume from branch
 * **/
void 
InsMC::from_branch(){
	//no need to do anything
	branch_id = -1;
	Q_size--;
}

/**
 * compute the upper and lower bound
 * **/
void 
InsMC::compute_bound(){
	//upper bound is easy
	int which_vet = tmp_R[branch_id];
	int color = C_map[which_vet];
	upper_bound = Q_size + color;
	//lower bound is itself 
	lower_bound = Q_size;
	//printf("branch_id %d which_vet %d color %d upper_bound %d lower_bound %d\n", branch_id, which_vet, color, upper_bound, lower_bound);
}

/**
 * calculate the combinations 
 * and get number of branches
 * this is actually the core algorithm
 * **/
int 
InsMC::get_num_branches(){
	//number of possible branches is based on R and the adjacent vertices of v
	num_branch = 0;
	for(int i=0; i<R_size; i++)
		R_avail[tmp_R[i]] = true;
	//main body, check every vertex in R and its adjacencies
	for(int i=0; i<R_size; i++){
		int v = tmp_R[i];
		int start = v==0?0:v_idx[v-1]; 
		int end = v_idx[v];
		for(int j=start; j<end; j++){
			int to = e_idx[j];
			if(R_avail[to] == true){
				branch_v[num_branch++] = v;
				break;
			}
		}
	}
	//resume	
	for(int i=0; i<R_size; i++)
		R_avail[tmp_R[i]] = false;
	return num_branch;
}

/**
 * from adj to encode
 * for the storage of the graph
 * **/
int 
InsMC::get_encode(int *encode){
	for(int i=0; i<Q_size; i++)
		encode[i] = Q[i];
	return Q_size;
}

/**
 * encode is stored in CSR format then transformed into 
 * adj format
 * **/
void 
InsMC::to_encode(int *encode, int size){
	for(int i=0; i<size; i++) Q[i] = encode[i];
	Q_size = size;
	//for(int i=0; i<Q_size; i++)
	//	printf("%d ", Q[i]);
	//printf("Q_size %d\n", Q_size);
	get_full_R();
	for(int i=0; i<Q_size; i++)
		join_R(Q[i]);
	get_deg_sorted();
	assign_C();
}

/**
 * print the content of the encode
 * **/
void 
InsMC::print_encode(){
	printf("encode: ");
	for(int i=0; i<Q_size; i++) printf("%d ", Q[i]);
	printf("\n");
}

/**
 * copy the 'best' encode into current code
 * **/
void 
InsMC::copy_code(int *encode, int size){
	for(int i=0; i<size; i++) Q[i] = encode[i];
	Q_size = size;
}

/**
 * get the best value
 * **/
int 
InsMC::get_value(){
	return upper_bound;
}

/**
 * get number of possible elements in the search list
 * **/
void
InsMC::get_num_elem(){
	num_count = upper_bound;
}

/**
 * private information is for debugging purposes
 * **/
int 
InsMC::return_private_info(){
	if(Q_size<3)
		return 0;
	//Utils::q_sort(Q, 0, Q_size-1);
	//bool is_clique = true;
	//for(int i=0; i<Q_size; i++){
	//	int v = Q[i];
	//	int start = v==0?0:v_idx[v-1];
	//	int end = v_idx[v];
	//	int m=0;
	//	int n=start;
	//	int num_common = 0;
	//	while(m<Q_size && n<end){
	//		if(Q[m]==e_idx[n]){
	//			m++;
	//			n++;
	//			num_common++;
	//		}
	//		else if(Q[m]<e_idx[n])
	//			m++;
	//		else if(Q[m]>e_idx[n])
	//			n++;
	//	}
	//	if(num_common<Q_size-1)
	//		is_clique = false;
	//}
	//if(is_clique==false){
	for(int i=0; i<Q_size; i++)
		fprintf(writer, "%d ", Q[i]);
	fprintf(writer, "\n");
	fflush(writer);
	//}
	return 0;
}

/***************************************************
 * these are private functions
 * for graph operations
 * *************************************************/

/**
 * Initiate the Ins_Dis first
 * **/
void
InsMC::init_graph(){
	//scan the graph to get the basic information
	scan_graph();
	read_graph();
}

/**
 * scan the graph to get some basic information
 * **/
void
InsMC::scan_graph(){
	//declare variables
	FILE *stream;
	char file_buf[500];
	int v_num;
	int e_num;
	int g_num;
	int v_id=0;
	int v_to=0;
	int sum=0;
	int color;

	//read the head information, and allocate according variables
	//printf("file %s\n", file);
	sprintf(file_buf, "%s", file);
	if((stream=fopen(file_buf,"r"))==NULL){
		printf("the file %s you input does not exist!\n", file_buf);
		exit(1);
	}
	if(fscanf(stream, "%d %d %d\n", &v_num, &g_num, &e_num)==EOF)
		printf("error here\n");
	v_size = v_num;
	e_size = e_num;
	v_idx = (int*)malloc(sizeof(int)*(v_size+1));
	e_idx = (int*)malloc(sizeof(int)*(e_size+1));
	Q = (int*)malloc(sizeof(int)*(v_size+1));
	R = (int*)malloc(sizeof(int)*(v_size+1));
	tmp_R = (int*)malloc(sizeof(int)*(v_size+1));
	R_deg = (int*)malloc(sizeof(int)*(v_size+1));
	R_avail = (bool*)malloc(sizeof(bool)*(v_size+1));
	C = (int*)malloc(sizeof(int)*(v_size+1));
	C_type = (int*)malloc(sizeof(int)*(v_size+1));
	C_map = (int*)malloc(sizeof(int)*(v_size+1));
	branch_v = (int*)malloc(sizeof(int)*(v_size+1));
	//scan the real content
	for(int i=0; i<v_size; i++) v_idx[i] = 0;
	for(int i=0;i<e_num;i++)
	{
		if(fscanf(stream, "%d %d %d\n",&v_id , &v_to, &color)==EOF)
			printf("error here\n");
		v_idx[v_id] +=1;
	}
	sum=0;
	for(int i=0; i<v_num; i++){
		sum += v_idx[i];
		v_idx[i]=sum;
	}
	//allocate adjs
	fclose(stream);
}

/**
 * read the real content  into CSR
 * **/
void
InsMC::read_graph(){
	FILE *stream;
	char file_buf[500];
	int v_num;
	int e_num;
	int g_num;
	int v_id=0;
	int v_to=0;
	int C=0;
	int i,j;
	//read basic information
	sprintf(file_buf, "%s", file);
	if((stream=fopen(file_buf,"r"))==NULL){
		printf("the file %s you input does not exist!\n", file_buf);
		exit(1);
	}
	if(fscanf(stream, "%d %d %d\n", &v_num, &g_num, &e_num)==EOF)
		printf("error here\n");
	//this is for the start of the different postitions
	int *idx=(int*)malloc(sizeof(int)*(v_size+1));
	for(j=0;j<v_size;j++){
		idx[j] = j==0?0:v_idx[j-1];
	}
	//real read part
	for(i=0;i<e_num;i++)
	{
		if(fscanf(stream, "%d %d %d\n",&v_id , &v_to, &C)==EOF)
			printf("error here\n");
		int pos = idx[v_id];
		e_idx[pos] = v_to;
		idx[v_id]++;
	}
	//e_idx also needs to be sorted, for the purpose of join operations
	for(int i=0; i<v_size; i++){
		int start = i==0?0:v_idx[i-1];
		int end = v_idx[i];
		Utils::q_sort(e_idx, start, end-1);
	}
	//frees
	free(idx);
	fclose(stream);
}

/**
 * detect connected components
 * and output them to file
 * use BFS algorithms
 * **/
void 
InsMC::separate_CC(){
	bool *visited = (bool*)malloc(sizeof(bool)*(v_size+1));
	for(int i=0; i<v_size; i++) visited[i] = false;
	int num_visited =0;
	int *queue = (int*)malloc(sizeof(int)*(v_size+1));
	int q_idx = 0;
	char f_buf[100];
	int count =0;
	while(true){
		count++;
		//detect a un-visited vertex
		int seed = -1;
		for(int i=0; i<v_size; i++){
			if(visited[i]==false){
				int start = i==0?0:v_idx[i-1];
				int end = v_idx[i];
				int deg = end-start;
				if(deg>=1){
					seed = i;
					visited[seed] = true;
					break;
				}
				else
					visited[i] = true;
			}
		}
		if(seed==-1)
			break;
		//perform bfs
		sprintf(f_buf,"%s%d.gr", CC_folder, count);
		FILE *writer = fopen(f_buf, "w");
		q_idx=0;
		queue[q_idx++] = seed;
		int edge_count =0;
		while(q_idx>0){
			int v = queue[--q_idx];
			int start = v==0?0:v_idx[v-1];
			int end = v_idx[v];
			for(int i=start; i<end; i++){
				if(visited[e_idx[i]]==false){
					visited[e_idx[i]] = true;
					queue[q_idx++] = e_idx[i];  
					fprintf(writer, "%d %d\n", v, e_idx[i]);
					fprintf(writer, "%d %d\n", e_idx[i], v);
					edge_count += 2;
				}
			}
		}
		fclose(writer);
	}
	printf("there are %d number of CCs in total\n", count);
}

/**
 * initialize the full R
 * **/
void
InsMC::get_full_R(){
	for(int i=0; i<v_size; i++)
		R[i] = i;
	R_size = v_size;
}

/**
 * at the very begining R is the full list of all vertices
 * when a v is added into Q
 * Q = Q \union v
 * R = R \join \gamma(v)
 * \gamma(v) is the set of vertices that adjacent to v
 * **/
void
InsMC::join_R(int v){
	int m = 0;
	int n = v==0?0:v_idx[v-1];
	int v_end = v_idx[v];
	R_size =0;
	while(m<v_size && n<v_end){
		if(R[m] == e_idx[n]){
			R[R_size++] = R[m];
			m++;
			n++;
		}
		else if(R[m] > e_idx[n])
			n++;
		else if(R[m] < e_idx[n])
			m++;
	}
	//keep R sorted
	Utils::q_sort(R, 0, R_size-1);
}

/**
 * this is to get a sorted list of each vertex in the graph
 * **/
void
InsMC::get_deg_sorted(){
	//initialize some info
	for(int i=0; i<R_size; i++){ 
		tmp_R[i] = R[i];
		R_avail[R[i]] = true;
	}
	//get degrees first
	for(int i=0; i<R_size; i++){
		int v = R[i];
		int start = v==0?0:v_idx[v-1];
		int end = v_idx[v];
		R_deg[i] = (end - start)*-1;
	}
	//then sort
	Utils::q_sort_two(R_deg, tmp_R, 0, R_size-1);
	for(int i=0; i<R_size; i++)
		R_avail[R[i]] = false;
}

/**
 * assign C to each vertex used for the upper bound update
 * using heuristic mentioned in Janes Konc's paper:
 * An improved branch and bound algorithm for the maximum clique problem
 * some improvement might be possible using vertex cover
 * **/
void
InsMC::assign_C(){
	int total_added=0;
	int now_C = 1;
	for(int i=0; i<R_size; i++){
		R_avail[tmp_R[i]] = true;
		C_map[tmp_R[i]] = -1;
	}
	//list each tree
	while(total_added < R_size){
		//***you have to put a first vertex in the C 
		//find the first, it's neccessary to have this
		int first = -1;
		for(int i=0; i<R_size; i++) {
			if(C_map[tmp_R[i]] == -1){
				first = i;
				break;
			}
		}
		C[first] = tmp_R[first];
		C_map[tmp_R[first]] = now_C;
		//printf("add %d color %d\n", tmp_R[first], now_C);
		total_added++;

		//***for each round, scan the vertex list
		for(int i=0; i<R_size; i++){
			if(i == first) continue;
			//each vertex will try to see if it's adjacet is already in
			int v = tmp_R[i];
			int start = v==0?0:v_idx[v-1];
			int end = v_idx[v];
			bool is_independent = true;
			for(int j=start; j<end; j++){
				//printf("---from %d to %d\n", v, e_idx[j]);
				if(R_avail[e_idx[j]]==true && C_map[e_idx[j]] == now_C){
					//printf("v %d c_map %d\n", v, C_map[v]);
					//printf("break !\n\n");
					is_independent = false;
					break;
				}
			}
			if(is_independent == true && C_map[v] == -1){
				C[total_added] = v;
				C_type[total_added] = now_C;
				C_map[v] = now_C;
				//printf("*add %d color %d tmp_R[i] %d\n", v, now_C, tmp_R[i]);
				total_added++;
			}
		}
		now_C++;
	}
	for(int i=0; i<R_size; i++) R_avail[tmp_R[i]] = false;
	max_C = now_C-1;
}


void 
InsMC::get_num_triangles(){
	int num_triangle = 0;
	int *pos_larger = (int*)malloc(sizeof(int)*v_size);
	for(int i=0; i<v_size; i++){
		int start = i==0?0:v_idx[i-1];
		int end = v_idx[i];
		for(int j=start; j<end; j++){
			if(e_idx[j]>i){
				pos_larger[i] = j;
				break;
			}
		}
	}
	for(int i=0; i<v_size; i++){
		if(i%10==0)
			printf("process %dth vertex\n", i);
		int start = i==0?0:v_idx[i-1];
		int end = v_idx[i];
		for(int j=start; j<end; j++){
			int to = e_idx[j];
			int to_start = to==0?0:v_idx[to-1];
			int to_end = v_idx[to];
			int m = start+j; 
			int n = pos_larger[to];
			while(m<end && n<to_end){
				if(e_idx[m]==e_idx[n]){
					m++;
					n++;
					num_triangle++;
				}
				else if(e_idx[m]<e_idx[n])
					m++;
				else if(e_idx[m]>e_idx[n])
					n++;
			}
		}
	}
	printf("total number of triangles %d\n", num_triangle);
}

/****************************************************
 * for debugging
 * ***************************************************/
/**
 * visualize csr graph
 * **/
void
InsMC::vis_csr(){
	FILE *writer = fopen("vis/vis_csr.dot", "w");
	fprintf(writer, "graph{\n");
	for(int i=0; i<v_size; i++){
		int start = i==0?0:v_idx[i-1];
		int end = v_idx[i];
		for(int j=start; j<end; j++){
			fprintf(writer, "%d -- %d [color = red]\n", i, e_idx[j]);
		}
	}
	fprintf(writer, "}\n");
}

/**
 * print degree information
 * **/
void 
InsMC::print_deg(){
	for(int i=0; i<R_size; i++){
		printf("%d:%d ", tmp_R[i], R_deg[i]);
	}
	printf("\n");
}


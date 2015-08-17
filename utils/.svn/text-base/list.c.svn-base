#include "list.h"
#include "graph.h"

void copy_tle(tle from, tle to);
void copy_tbe(tbe from, tbe to);
void
double_list_child(psl l, int buck_id);
void
double_list_parent(psl l, int buck_id);
void 
from_file(psl l, int buck_id);
void 
to_file(psl l, int buck_id);

extern int DEBUG_RETRIEVE;

void 
reset_list(psl l){
	int i, j;
	for(i=0;i<l->buck_size;i++){
		l->idx_child[i] = 0;
		l->idx_parent[i] = 0;
		for(j=0; j<l->list_size[i]; j++)
			l->parent_size[i][j]=0;
	}
	//also need to remove list stored in the file
}

void free_list(psl l){
	int i, j, k;
	for(i=0; i<l->buck_size; i++){
		for(j=0; j<l->list_size[i]; j++){
			free(l->parent[i][j]);
		}
		free(l->parent[i]);
		free(l->parent_size[i]);
		free(l->child[i]);
	}
	free(l->idx_child);
	free(l->child_size);
	free(l->child);
	free(l->parent);
	free(l->parent_size);
	free(l->idx_parent);
	free(l->f_check);
	free(l->list_size);
	free(l);
}

void
create_search_list(int buck_size, int list_size, 
		int p_size, psl l){
	int i,j,k;
	l->buck_size = buck_size;
	l->list_size = (int*)malloc(sizeof(int)*buck_size);	
	for(i=0;i<buck_size; i++)
		l->list_size[i] = list_size;	
	l->p_size = p_size;
	l->idx_parent = (int*)malloc(sizeof(int)*buck_size);
	l->idx_child = (int*)malloc(sizeof(int)*buck_size);
	l->parent = (int***)malloc(sizeof(int**)*buck_size);
	l->parent_size = (int**)malloc(sizeof(int*)*buck_size);
	l->child_size = (int*)malloc(sizeof(int)*buck_size);
	l->child = (tle**)malloc(sizeof(tle*)*buck_size);	
	l->f_check = (int*)malloc(sizeof(int)*buck_size); 
	for(i=0; i<buck_size; i++){
		l->idx_child[i] = 0;
		l->idx_parent[i] = 0;
		l->child_size[i] = list_size*AVG_CHILD_PER_PARENT;
		l->parent[i] = (int**)malloc(sizeof(int*)*list_size);
		l->parent_size[i] = (int*)malloc(sizeof(int)*list_size);
		l->child[i] = (tle*)malloc(sizeof(tle)*list_size*AVG_CHILD_PER_PARENT);
		for(j=0;j<list_size;j++)
			l->parent[i][j] = (int*)malloc(sizeof(int)*p_size);
		l->f_check[i] = 0;
	}
	
}

int
add_parent(pg g, psl l, 
		int buck_id){
	int i;
	int buck_pos = l->idx_parent[buck_id];
	l->idx_parent[buck_id]++;
	if(l->idx_parent[buck_id] == l->list_size[buck_id]){
		printf("buckpos %d list idx is %d list size is %d\n", buck_id, l->idx_parent[buck_id], l->list_size[buck_id]);
		if(l->list_size[buck_id]<1000000)
			double_list_parent(l, buck_id);
		else{
			to_file(l, buck_id);
			buck_pos = l->idx_parent[buck_id];
			l->idx_parent[buck_id]++;
		}
	}
	for(i=0;i<g->e_size[3];i++){
		l->parent[buck_id][buck_pos][i] = g->e_idx[3][i];
	}
	l->parent_size[buck_id][buck_pos] = 0;
	return buck_pos;
}

void
add_child(int from_one, int to_one, 
		int from_two, int to_two, 
		int conn_type, int idx_parent,
		int buck_id, psl l){
	int buck_pos_child = l->idx_child[buck_id];
	l->idx_child[buck_id]++;
	//it's possible to add some code here to increment the size of child bucket if overflowed
	if(l->idx_child[buck_id]>=l->child_size[buck_id]){
		double_list_child(l, buck_id);
	}
	//set child
	l->child[buck_id][buck_pos_child].from_one = from_one;
	l->child[buck_id][buck_pos_child].to_one = to_one;
	l->child[buck_id][buck_pos_child].from_two = from_two;
	l->child[buck_id][buck_pos_child].to_two = to_two;
	l->child[buck_id][buck_pos_child].new_connect_type = conn_type;
	l->child[buck_id][buck_pos_child].idx_parent = idx_parent;
#ifdef USE_DBG
	//if(buck_id == 0)
	//printf("buck_id %d buck_pos %d parent size %d from_one %d, to_one %d, from_two %d, to_two %d, type %s\n", 
	//		buck_id, idx_parent, l->parent_size[buck_id][idx_parent], from_one, to_one, from_two, to_two, conn_type==CONN_TYPE_TWO?"one":"two");
#endif
	//increment parent
	l->parent_size[buck_id][idx_parent]++;
}

/* ***********************************************
 * it need to be doubled by both child and parent.
 * ***********************************************/
void
double_list_child(psl l, int buck_id){
	printf("doubling child list\n");
	int i, j;
	//for child
	tle *tmp_child = (tle*)malloc(sizeof(tle)*l->child_size[buck_id]);
	for(i=0;i<l->child_size[buck_id];i++){
		copy_tle(l->child[buck_id][i], tmp_child[i]);
	}
	l->child[buck_id] = (tle*)realloc(l->child[buck_id], sizeof(tle)*l->child_size[buck_id]*2); //double the size
	//copy back
	for(i=0;i<l->child_size[buck_id];i++){
		copy_tle(tmp_child[i], l->child[buck_id][i]);
	}
	l->child_size[buck_id] = l->child_size[buck_id]*2;
	free(tmp_child);
}

/*
 *
 * **/
void 
to_file(psl l, int buck_id){
	int i, j;
	char f_buf[100];
	sprintf(f_buf, "./tmp/local_%d_%d", buck_id, l->f_check[buck_id]);
	printf("start writing file to %s, there are %d child and %d parent\n", f_buf, l->idx_child[buck_id], l->idx_parent[buck_id]);
	FILE *writer = fopen(f_buf, "wb");
	fwrite(&l->idx_child[buck_id], sizeof(int), 1, writer);
	fwrite(&l->idx_parent[buck_id], sizeof(int), 1, writer);
	fwrite(l->child[buck_id], sizeof(tle), l->idx_child[buck_id], writer);
	fwrite(l->parent_size[buck_id], sizeof(int), l->idx_parent[buck_id], writer);
	for(i=0;i<l->idx_parent[buck_id]; i++)
		fwrite(l->parent[buck_id][i], sizeof(int), l->p_size, writer);
	fclose(writer);
	//restore info
	l->idx_child[buck_id] = 0;
	l->idx_parent[buck_id] = 0;
	for(j=0; j<l->list_size[buck_id]; j++)
		l->parent_size[buck_id][j]=0;
	//update file info
	l->f_check[buck_id]++;
}

/*
 *
 * **/
void 
from_file(psl l, int buck_id){
	//update file info
	l->f_check[buck_id]--;
	int i, j;
	char f_buf[100];
	sprintf(f_buf, "./tmp/local_%d_%d", buck_id, l->f_check[buck_id]);
	FILE *reader = fopen(f_buf, "rb");
	if(fread(&l->idx_child[buck_id], sizeof(int), 1, reader)==EOF)
		printf("error here\n");
	if(fread(&l->idx_parent[buck_id], sizeof(int), 1, reader)==EOF)
		printf("error here\n");
	if(fread(l->child[buck_id], sizeof(tle), l->idx_child[buck_id], reader)==EOF)
		printf("error here\n");
	if(fread(l->parent_size[buck_id], sizeof(int), l->idx_parent[buck_id], reader)==EOF)
		printf("error here\n");
	for(i=0;i<l->idx_parent[buck_id]; i++)
		if(fread(l->parent[buck_id][i], sizeof(int), l->p_size, reader)==EOF)
			printf("error here\n");
	printf("getting info from file %s there are %d child and %d parent\n", f_buf, l->idx_child[buck_id], l->idx_parent[buck_id]);
	fclose(reader);
}

void
double_list_parent(psl l, int buck_id){
	printf("doubling parent list\n");
	int i, j, k;
	//for parent
	int list_size = l->list_size[buck_id];
	int **tmp_parent = (int**)malloc(sizeof(int*)*list_size);
	for(i=0;i<list_size; i++)
		tmp_parent[i] = (int*)malloc(sizeof(int)*l->p_size);
	for(i=0; i<list_size;i++)
		for(j=0;j<l->p_size; j++)
			tmp_parent[i][j] = l->parent[buck_id][i][j];
	for(i=0;i<list_size; i++)
		free(l->parent[buck_id][i]);
	l->parent[buck_id] = (int**)realloc(l->parent[buck_id], sizeof(int*)*list_size*2);
	for(i=0;i<list_size*2; i++)
		l->parent[buck_id][i] = (int*)malloc(sizeof(int)*l->p_size);
	for(i=0;i<list_size; i++){
		for(j=0;j<l->p_size; j++)
			l->parent[buck_id][i][j] = tmp_parent[i][j];
		free(tmp_parent[i]);
	}
	free(tmp_parent);
	int *tmp_parent_size = (int*)malloc(sizeof(int)*list_size);
	for(i=0;i<list_size; i++)
		tmp_parent_size[i] = l->parent_size[buck_id][i];
	l->parent_size[buck_id] = (int*)realloc(l->parent_size[buck_id], sizeof(int)*list_size*2);
	for(i=0;i<list_size; i++)
		l->parent_size[buck_id][i] = tmp_parent_size[i];
	free(tmp_parent_size);
	l->list_size[buck_id] = list_size*2;
}

void copy_tle(tle from, tle to){
	from.from_one = to.from_one;
	from.to_one = to.to_one;
	from.from_two = to.from_two;
	from.to_two = to.to_two;
	from.new_connect_type = to.new_connect_type;
	from.idx_parent = to.idx_parent;
}

int 
retrive(pg g, psl l, 
		int buck_id){
	int i,j,k;
	int buck_pos_child = l->idx_child[buck_id]-1;
	if(buck_pos_child<=0){
		printf("buck_id %d file num %d %d\n", buck_id, l->f_check[buck_id], l->f_check[buck_id+1]);
		if(l->f_check[buck_id]>0){
			from_file(l, buck_id);
			buck_pos_child = l->idx_child[buck_id]-1;
		}
		else
			return FALSE;
	}
	l->idx_child[buck_id]--;
	int buck_pos_parent = l->child[buck_id][buck_pos_child].idx_parent;
	//get child and parent and do the operation on graph
	for(i=0; i<g->e_size[3];i++)
		g->e_idx[3][i] = l->parent[buck_id][buck_pos_parent][i];
	//perform the kl switch
	int from_one = l->child[buck_id][buck_pos_child].from_one;
	int to_one = l->child[buck_id][buck_pos_child].to_one;
	int from_two = l->child[buck_id][buck_pos_child].from_two;
	int to_two = l->child[buck_id][buck_pos_child].to_two;
	int conn_type = l->child[buck_id][buck_pos_child].new_connect_type;
//#ifdef USE_DBG
	//if(buck_pos_child==3483){
	//printf("buck_id %d buck_pos %d num_parent %d from_one %d, to_one %d, from_two %d, to_two %d type %s ", 
	//	buck_id, buck_pos_parent, l->parent_size[buck_id][buck_pos_parent], from_one, to_one, from_two, to_two, conn_type==CONN_TYPE_ONE?"one":"two");	
	//graph_vis_two(g, "vis0.dot", 1, 3);
	//}
//	printf("buck_pos_child %d\n", buck_pos_child);
//#endif
	//printf("pos_child %d from_one %d to_one %d from_two %d to_two%d conn type %d\n", 
	//	buck_pos_child, from_one, from_two, to_one, to_two, conn_type);
	int start_from_one = from_one==0?0:g->v_idx[3][from_one-1];
	int end_from_one = g->v_idx[3][from_one];
	int start_to_one = to_one==0?0:g->v_idx[3][to_one-1];
	int end_to_one = g->v_idx[3][to_one];
	int start_from_two = from_two==0?0:g->v_idx[3][from_two-1];
	int end_from_two = g->v_idx[3][from_two];
	int start_to_two = to_two==0?0:g->v_idx[3][to_two-1];
	int end_to_two = g->v_idx[3][to_two];
	if(conn_type == CONN_TYPE_ONE){
		for(i=start_from_one;i<end_from_one;i++){
			if(g->e_idx[3][i] == to_one){
				g->e_idx[3][i] = from_two;
				if(DEBUG_RETRIEVE == TRUE)
				printf("now %d is connected to %d\n", from_one, from_two);
				break;
			}
		}
		for(i=start_to_one;i<end_to_one;i++){
			if(g->e_idx[3][i] == from_one){
				g->e_idx[3][i] = to_two;
				if(DEBUG_RETRIEVE == TRUE)
				printf("now %d is connected to %d\n", to_one, to_two);
				break;
			}
		}
		for(i=start_from_two;i<end_from_two;i++){
			if(g->e_idx[3][i] == to_two){
				g->e_idx[3][i] = from_one;
				if(DEBUG_RETRIEVE == TRUE)
				printf("now %d is connected to %d\n", from_two, from_one);
				break;
			}
		}
		for(i=start_to_two;i<end_to_two;i++){
			if(g->e_idx[3][i] == from_two){
				g->e_idx[3][i] = to_one;
				if(DEBUG_RETRIEVE == TRUE)
				printf("now %d is connected to %d\n", to_two, to_one);
				break;
			}
		}
	}
	else if(conn_type == CONN_TYPE_TWO){
		for(i=start_from_one;i<end_from_one;i++)
			if(g->e_idx[3][i] == to_one)
				g->e_idx[3][i] = to_two;
		for(i=start_to_one;i<end_to_one;i++)
			if(g->e_idx[3][i] == from_one)
				g->e_idx[3][i] = from_two;
		for(i=start_from_two;i<end_from_two;i++)
			if(g->e_idx[3][i] == to_two)
				g->e_idx[3][i] = to_one;
		for(i=start_to_two;i<end_to_two;i++)
			if(g->e_idx[3][i] == from_two)
				g->e_idx[3][i] = from_one;
		
	}	
	//decrease parent
	l->parent_size[buck_id][buck_pos_parent]--;
	if(l->parent_size[buck_id][buck_pos_parent]==0){
		//there might be no need to do this step
		for(i=0; i<l->p_size;i++)
			l->parent[buck_id][buck_pos_parent][i] = -1;
		l->idx_parent[buck_id]--;
	}
	return TRUE;
}

int 
get_list_size(psl l){
	int i;
	int result = 0;
	for(i=0; i<l->buck_size; i++){
		result += l->idx_child[i];
	}
	return result;
}

void 
shrink(pg g, int from_one,
		int to_one, int from_two,
		int to_two, int conn_type){
	int i,j,k;
	//make change
	int start_from_one = from_one==0?0:g->v_idx[3][from_one-1];
        int end_from_one = g->v_idx[3][from_one]; 
        int start_to_one = to_one==0?0:g->v_idx[3][to_one-1];
        int end_to_one = g->v_idx[3][to_one];
        int start_from_two = from_two==0?0:g->v_idx[3][from_two-1];
        int end_from_two = g->v_idx[3][from_two];
        int start_to_two = to_two==0?0:g->v_idx[3][to_two-1]; 
        int end_to_two = g->v_idx[3][to_two];
        if(conn_type == CONN_TYPE_ONE){
                for(i=start_from_one;i<end_from_one;i++)
                        if(g->e_idx[3][i] == to_one)
                                g->e_idx[3][i] = from_two;
                for(i=start_to_one;i<end_to_one;i++)
                        if(g->e_idx[3][i] == from_one)
                                g->e_idx[3][i] = to_two;
                for(i=start_from_two;i<end_from_two;i++)
                        if(g->e_idx[3][i] == to_two)
                                g->e_idx[3][i] = from_one;
                for(i=start_to_two;i<end_to_two;i++)
                        if(g->e_idx[3][i] == from_two)
                                g->e_idx[3][i] = to_one;
        }
        else if(conn_type == CONN_TYPE_TWO){
                for(i=start_from_one;i<end_from_one;i++)
                        if(g->e_idx[3][i] == to_one)
                                g->e_idx[3][i] = to_two;
                for(i=start_to_one;i<end_to_one;i++)
                        if(g->e_idx[3][i] == from_one)
                                g->e_idx[3][i] = from_two;
                for(i=start_from_two;i<end_from_two;i++)
                        if(g->e_idx[3][i] == to_two)
                                g->e_idx[3][i] = to_one;
                for(i=start_to_two;i<end_to_two;i++)
                        if(g->e_idx[3][i] == from_two)
                                g->e_idx[3][i] = from_one;

        }
}

void 
resume(pg g, int from_one,
		int to_one, int from_two,
		int to_two, int conn_type){
	int i,j,k;
	int start_from_one = from_one==0?0:g->v_idx[3][from_one-1];
        int end_from_one = g->v_idx[3][from_one];
        int start_to_one = to_one==0?0:g->v_idx[3][to_one-1];
        int end_to_one = g->v_idx[3][to_one];
        int start_from_two = from_two==0?0:g->v_idx[3][from_two-1];
        int end_from_two = g->v_idx[3][from_two];
        int start_to_two = to_two==0?0:g->v_idx[3][to_two-1];
        int end_to_two = g->v_idx[3][to_two];
        if(conn_type == CONN_TYPE_ONE){
                for(i=start_from_one;i<end_from_one;i++)
                        if(g->e_idx[3][i] == from_two)
                                g->e_idx[3][i] = to_one;
                for(i=start_to_one;i<end_to_one;i++)
                        if(g->e_idx[3][i] == to_two)
                                g->e_idx[3][i] = from_one;
                for(i=start_from_two;i<end_from_two;i++)
                        if(g->e_idx[3][i] == from_one)
                                g->e_idx[3][i] = to_two;
                for(i=start_to_two;i<end_to_two;i++)
                        if(g->e_idx[3][i] == to_one)
                                g->e_idx[3][i] = from_two;
        } 
        else if(conn_type == CONN_TYPE_TWO){
                for(i=start_from_one;i<end_from_one;i++)
                        if(g->e_idx[3][i] == to_two)
                                g->e_idx[3][i] = to_one;
                for(i=start_to_one;i<end_to_one;i++)
                        if(g->e_idx[3][i] == from_two)
                                g->e_idx[3][i] = from_one;
                for(i=start_from_two;i<end_from_two;i++)
                        if(g->e_idx[3][i] == to_one)
                                g->e_idx[3][i] = to_two;
                for(i=start_to_two;i<end_to_two;i++)
                        if(g->e_idx[3][i] == from_one)
                                g->e_idx[3][i] = from_two;
        }	
}

void
create_search_list_bnb(int buck_size, int list_size, 
		int p_size, pbl l){
	int i,j,k;
	l->buck_size = buck_size;
	l->list_size = list_size;	
	l->p_size = p_size;
	l->idx_parent = (int*)malloc(sizeof(int)*buck_size);
	l->idx_child = (int*)malloc(sizeof(int)*buck_size);
	l->parent = (int***)malloc(sizeof(int**)*buck_size);
	l->parent_size = (int**)malloc(sizeof(int*)*buck_size);
	l->child_size = (int*)malloc(sizeof(int)*buck_size);
	l->child = (tbe**)malloc(sizeof(tbe*)*buck_size);	
	for(i=0; i<buck_size; i++){
		l->idx_child[i] = 0;
		l->idx_parent[i] = 0;
		l->child_size[i] = list_size*AVG_CHILD_PER_PARENT;
		l->parent[i] = (int**)malloc(sizeof(int*)*list_size);
		l->parent_size[i] = (int*)malloc(sizeof(int)*list_size);
		l->child[i] = (tbe*)malloc(sizeof(tbe)*list_size*AVG_CHILD_PER_PARENT);
		for(j=0;j<list_size;j++)
			l->parent[i][j] = (int*)malloc(sizeof(int)*p_size);
	}
}

void
add_nodes_bnb(pg g, pbl l, int base,
		int *upper_bound, int* lower_bound){
	int i,j;
	for(i=0;i<g->v_size; i++){
		if(g->v_check[i]==NUL)
			continue;
		for(j=0;j<g->v_size; j++){
			if(g->v_check[j]==NUL)
				continue;
			int from = i;
			int to = j;
			if(from != to){
				//then each from and to combinnation
			}

		}
	}
}

int 
retrive_bnb(pg g, pbl l, int buck_id){
	int i,j,k;
	int buck_pos_child = l->idx_child[buck_id]-1;
	if(buck_pos_child==0)
		return FALSE;
	l->idx_child[buck_id]--;
	int buck_pos_parent = l->child[buck_id][buck_pos_child].idx_parent;
	//get child and parent and do the operation on graph
	for(i=0; i<g->e_size[3];i++)
		g->e_idx[3][i] = l->parent[buck_id][buck_pos_parent][i];
	//perform the kl switch
	int from = l->child[buck_id][buck_pos_child].from;
	int to = l->child[buck_id][buck_pos_child].to;
	int start_from = from==0?0:g->v_idx[3][from-1];
	int end_from = g->v_idx[3][from];
	int start_to = to==0?0:g->v_idx[3][to-1];
	int end_to = g->v_idx[3][to];
	//might need to be changed to cope with duplications
	for(i=start_from; i<end_from; i++)
		if(g->e_idx[3][i] == to)
			g->e_idx[3][i] = from;
	for(i=start_to;i<end_to;i++)
		if(g->e_idx[3][i] == from)
			g->e_idx[3][i] = to;
	//decrease parent
	l->parent_size[buck_id][buck_pos_parent]--;
	if(l->parent_size[buck_id][buck_pos_parent]==0){
		//there might be no need to do this step
		for(i=0; i<l->p_size;i++)
			l->parent[buck_id][buck_pos_parent][i] = -1;
		l->idx_parent[buck_id]--;
	}
	return TRUE;
}

void
double_list_bnb(pbl l, int buck_id){
	printf("I am doubling the list size on level %d\n", buck_id);
	tbe *tmp_child = (tbe*)malloc(sizeof(tbe)*l->child_size[buck_id]);
	int i;
	for(i=0;i<l->child_size[buck_id];i++){
		copy_tbe(l->child[buck_id][i], tmp_child[i]);
	}
	l->child[buck_id] = (tbe*)realloc(l->child[buck_id], sizeof(tbe)*l->child_size[buck_id]*2); //double the size
	//copy back
	for(i=0;i<l->child_size[buck_id];i++){
		copy_tbe(tmp_child[i], l->child[buck_id][i]);
	}
	l->child_size[buck_id] = l->child_size[buck_id]*2;
	free(tmp_child);
}

void copy_tbe(tbe from, tbe to){
	from.from = to.from;
	from.to = to.to;
}

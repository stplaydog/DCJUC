#include "listdis.h"
#include "instance.h"

#ifndef LIST_SIZE
#define LIST_SIZE 1000
#endif

extern int g_count;

DISList::DISList(myint b_size, myint list_size, 
		myint b, myint num_t, 
		bool is_ub, int num_elem) : List(b_size, list_size, b, num_t, is_ub){
	int i;
	content = new int*[b_size];
	for(i=0; i<b_size; i++)
    {
		content[i] = new int[LIST_SIZE*num_elem];
    }
	c_idx = new int*[b_size];
	for(i=0; i<b_size; i++)
    {
		c_idx[i] = new int[LIST_SIZE];
    }
	c_pos = new int[b_size];	
}

DISList::~DISList(){
	for(int i=0; i<buck_size; i++){
		free(content[i]);
		free(c_idx[i]);
	}
	free(content);
	free(c_idx);
	free(c_pos);
	//free original stuff
	for(int i=0; i<num_threads; i++){
		free(count[i]);
		free(idx_sum[i]);
		free(end[i]);
		free(start[i]);
		free(eliminate[i]);
		free(nei_score[i]);
	}
	free(count);
	free(idx_sum);
	free(end);
	free(start);
	free(eliminate);
	free(nei_score);
	//
	free(content_size);
    free(num);
    free(elim_idx);
}

void 
DISList::add(lint pos, myint buck_id, 
		int num_code, int *encode, int ins_id){
	int i;
	int stt = pos==0?0:c_idx[buck_id][pos-1];
	for(i=0; i<num_code; i++){
		content[buck_id][stt++] = encode[i];
	}
	//printf("add ");
	//for(int i=0; i<num_code; i+=2)
	//	printf("%d:%d ", encode[i], encode[i+1]);
	//printf("\n");
	c_idx[buck_id][pos] =stt;
#ifdef USE_DEBUG
    //add_g_count_list(buck_id, g_count);
#endif
}

void 
DISList::get(lint pos, myint buck_id, 
		int *encode, int *number, int ins_id){
    assert(buck_id>=0 && pos >=0);
	int i, j;
	int start = pos==0?0:c_idx[buck_id][pos-1];
	int end = c_idx[buck_id][pos];
	for(i=0,j=start;j<end; i++,j++)
		encode[i] = content[buck_id][j];
	*number = (end-start);
#ifdef USE_DEBUG
    //printf("######the g_count of this item is %d\n", g_count_list[buck_id][pos]);
#endif
}

void
DISList::prepare_parallel_list(int buck_id){
	int i, j;
	//get start position for each thread
	for (i=0; i<buck_size; i++){
		int sum = i==buck_id?0:num[i];
		for(j=0; j<num_threads; j++){
			start[j][i] = sum;
			sum += count[j][i];
			end[j][i] = sum;
		}
	}
	//compute array position
	for(i=0; i<buck_size; i++){
		if(i==buck_id)
			continue;
		int sum=c_idx[i][num[i]==0?0:num[i]-1];
		sum += idx_sum[0][i];
		for(j=1; j<num_threads; j++){
			int pre_pos = start[j][i]==0?0:start[j][i]-1;
			c_idx[i][pre_pos] = pre_pos==0?0:sum;
			sum += idx_sum[j][i];
		}
	}
	//cope with zero bucket
	if(is_ub==true){
		c_idx[0][0] = 0;
		int sum = idx_sum[0][buck_id];
		for(j=1; j<num_threads; j++){
			int pre_pos = start[j][buck_id]==0?0:start[j][buck_id]-1;
			c_idx[0][pre_pos] = pre_pos==0?0:sum;
			sum += idx_sum[j][buck_id];
		}
	}
	else{
		c_idx[buck_size-1][0] = 0;
		int sum = idx_sum[0][buck_id];
		for(j=1; j<num_threads; j++){
			int pre_pos = start[j][buck_id]==0?0:start[j][buck_id]-1;
			c_idx[buck_size-1][pre_pos] = pre_pos==0?0:sum;
			sum += idx_sum[j][buck_id];
		}
	}
}

int
DISList::copy_zero(int buck_id){
	int i;
	if(is_ub==true){
		num[0] = end[num_threads-1][buck_id];
		//copy back
		int copy_size = num[0]==0?0:c_idx[0][num[0]-1]; 
#pragma omp parallel for num_threads(num_threads) private(i)
		for(i=0; i<copy_size; i++){
			content[buck_id][i] = content[0][i];
		}
#pragma omp parallel for num_threads(num_threads) private(i)
		for(i=0; i<num[0]; i++){
			c_idx[buck_id][i] = c_idx[0][i];
			c_idx[0][i] = 0;
		}
		num[buck_id] = num[0];
		num[0]=0;
		return copy_size;
	}
	else{
		num[buck_size-1] = end[num_threads-1][buck_id];
		//copy back
		int copy_size = num[buck_size-1]==0?0:c_idx[buck_size-1][num[buck_size-1]-1]; 
#pragma omp parallel for num_threads(num_threads) private(i)
		for(i=0; i<copy_size; i++){
			content[buck_id][i] = content[buck_size-1][i];
		}
#pragma omp parallel for num_threads(num_threads) private(i)
		for(i=0; i<num[buck_size-1]; i++){
			c_idx[buck_id][i] = c_idx[buck_size-1][i];
			c_idx[buck_size-1][i] = 0;
		}
		num[buck_id] = num[buck_size-1];
		num[buck_size-1]=0;
		return copy_size;
	}
}

void
DISList::reset_num(){
	int i;
	//if there is anything, reset the number
	for(i=lower_bound; i<upper_bound; i++){
		if(start[0][i-base] < end[num_threads-1][i-base])
			num[i-base] = end[num_threads-1][i-base];
	}
}

void
DISList::copy(List *other, int buck_id, int size){
	DISList *l = (DISList*)other;
	int start = num[buck_id]-size;	
	int offset = c_idx[buck_id][start-1];
	//copy c_idx
	for(int i=start, j=0; i<num[buck_id]; i++, j++)
		l->c_idx[buck_id][j] = c_idx[buck_id][i]-offset;
	//copy content
	int pos_start = c_idx[buck_id][start-1];
	int pos_end = c_idx[buck_id][num[buck_id]-1];
	for(int i=pos_start, j=0; i<pos_end; i++, j++)
		l->content[buck_id][j] = content[buck_id][i];
	l->num[buck_id] = size;
	num[buck_id]-=size;
}


void
DISList::print_bucket(int buck_id){
	printf("buck %d num %lld:", buck_id, num[buck_id]);
	for (int i=0; i<num[buck_id]; i++){
		printf("|");
		int start = i==0?0:c_idx[buck_id][i-1];
		int end = c_idx[buck_id][i];
		for(int j=start; j<end; j++)
			printf("%d ", content[buck_id][j]);
	}
	printf("\n");
}


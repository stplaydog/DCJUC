#include "listmed.h"
#include "instance.h"

#ifndef LIST_SIZE
#define LIST_SIZE 1000000
#endif

MedList::MedList(myint b_size, myint list_size, 
		myint b, myint num_t, 
		bool is_ub, int num_elem) : List(b_size, list_size, b, num_t, is_ub){
    is_parallel=false;
	int i;
	content = (int**)malloc(sizeof(int*)*(b_size));
	for(i=0; i<b_size; i++)
		content[i] = (int*)malloc(sizeof(int)*LIST_SIZE*num_elem);
	c_idx = (int**)malloc(sizeof(int*)*(b_size));
	for(i=0; i<b_size; i++)
		c_idx[i] = (int*)malloc(sizeof(int)*LIST_SIZE);
	c_pos = (int*)malloc(sizeof(int)*b_size);	
}

MedList::~MedList(){
	for(int i=0; i<buck_size; i++){
                free(content[i]);
                free(c_idx[i]);
        }       
        free(content);
        free(c_idx);
        free(c_pos);
}

void 
MedList::add(lint pos, myint buck_id, 
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
}

void 
MedList::get(lint pos, myint buck_id, 
		int *encode, int *number, int ins_id){
	int i, j;
	int start = pos==0?0:c_idx[buck_id][pos-1];
	int end = c_idx[buck_id][pos];
	for(i=0,j=start;j<end; i++,j++)
		encode[i] = content[buck_id][j];
	*number = (end-start);
}

/**
 *this is the only function that needs to be changed 
 * **/
void
MedList::prepare_parallel_list(int buck_id){
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
	// cope with zero bucket
	// there is no need to cope with zero bucket
}

int
MedList::copy_zero(int buck_id){
	//no need to do anything
	return -1;
}

void
MedList::reset_num(){
	int i;
	//if there is anything, reset the number
	for(i=lower_bound; i<upper_bound; i++){
		if(start[0][i-base] < end[num_threads-1][i-base])
			num[i-base] = end[num_threads-1][i-base];
	}
}

void
MedList::copy(List *other, int buck_id, int size){
	MedList *l = (MedList*)other;
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
MedList::print_bucket(int buck_id){
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


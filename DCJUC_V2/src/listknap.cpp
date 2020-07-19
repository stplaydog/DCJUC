#include "listknap.h"
#include "instance.h"


KnapList::KnapList(myint b_size, myint list_size, 
		myint b, myint num_t, 
		bool is_ub, int num_elem) : List(b_size, list_size, b, num_t, is_ub)
{
	int i;
	content = (int**)malloc(sizeof(int*)*(b_size));
	for(i=0; i<b_size; i++)
    {
        content[i] = nullptr;
    }
	c_idx = (int**)malloc(sizeof(int*)*(b_size));
	for(i=0; i<b_size; i++)
    {
        c_idx[i] = nullptr;
    }
	c_pos = (int*)malloc(sizeof(int)*b_size);	
}

KnapList::~KnapList(){
}

void 
KnapList::add(lint pos, myint buck_id, 
		int num_code, int *encode, int ins_id)
{
	int i;
    /* if this array at this pos is not allocated, allocate the memory */
    if(content[buck_id] == nullptr && c_idx[buck_id] == nullptr)
    {
		content[buck_id] = (int*)malloc(sizeof(int)*LIST_SIZE*100);
		c_idx[buck_id] = (int*)malloc(sizeof(int)*LIST_SIZE);
        content_size[buck_id] = LIST_SIZE*100;
    }
	int stt = pos==0?0:c_idx[buck_id][pos-1];
    if(stt>=content_size[buck_id]-100)
    {
        reallocate_buck(buck_id);
    }
	for(i=0; i<num_code; i++)
    {
		content[buck_id][stt++] = encode[i];
    }
	c_idx[buck_id][pos] =stt;
	//printf("pos %d buck_id %d ins %d num_code %d\n", pos, buck_id, ins_id, num_code);
	//false sharing should be changed
		//num[buck_id]++;
	//if(encode[0]==7 && encode[6]==7){
	//	printf("duplicate here!\n");
	//	exit(1);
	//}
}

void 
KnapList::get(lint pos, myint buck_id, 
		int *encode, int *number, int ins_id)
{
	int i, j;
	int start = pos==0?0:c_idx[buck_id][pos-1];
	int end = c_idx[buck_id][pos];
	for(i=0,j=start;j<end; i++,j++)
    {
		encode[i] = content[buck_id][j];
    }
	*number = (end-start);
	//printf("%p\n", content);
	//fflush(stdout);
	//ins[ins_id]->to_encode(encode, (end-start));
}

void
KnapList::prepare_parallel_list(int buck_id){
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
	c_idx[0][0] = 0;
	int sum = idx_sum[0][buck_id];
	for(j=1; j<num_threads; j++){
		int pre_pos = start[j][buck_id]==0?0:start[j][buck_id]-1;
		c_idx[0][pre_pos] = pre_pos==0?0:sum;
		sum += idx_sum[j][buck_id];
	}
	//print
	//for(i=0; i<buck_size; i++)
	//	for(j=0; j<num_threads; j++){
	//		int pre_pos = start[j][i]==0?0:start[j][i]-1;
	//		printf("thread %d on bucket %d start %d end %d pre_pos %d arr_start %d occupied %d\n", j, i+base, start[j][i], end[j][i], pre_pos, i==buck_id?c_idx[0][pre_pos]:c_idx[i][pre_pos], idx_sum[j][i]);
	//	}
}

int
KnapList::copy_zero(int buck_id){
	int i;
	num[0] = end[num_threads-1][buck_id];
	//copy back
	int copy_size = num[0]==0?0:c_idx[0][num[0]-1]; 
#pragma omp parallel for num_threads(num_threads) private(i)
	for(i=0; i<copy_size; i++){
		content[buck_id][i] = content[0][i];
		//printf("%d ", content[0][i]);
	}
	//printf("\n");
#pragma omp parallel for num_threads(num_threads) private(i)
	for(i=0; i<num[0]; i++){
		c_idx[buck_id][i] = c_idx[0][i];
		c_idx[0][i] = 0;
	}
	num[buck_id] = num[0];
	num[0]=0;
	return copy_size;
}

void
KnapList::reset_num(){
	int i;
	for(i=lower_bound; i<upper_bound; i++){
		if(start[0][i-base] < end[num_threads-1][i-base])
		num[i-base] = end[num_threads-1][i-base];
	}
}

void
KnapList::copy(List *other, int buck_id, int size){
	KnapList *l = (KnapList*)other;
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
KnapList::print_bucket(int buck_id){
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

/**
 *
**/
void KnapList::free_buck(int buck_id)
{
    free(content[buck_id]); 
    free(c_idx[buck_id]); 
}

/**
 *
**/
void KnapList::reallocate_buck(int buck_id)
{
    int *tmp_cont = (int*) malloc(sizeof(int)*content_size[buck_id]*2);
    int *tmp_c_idx = (int*) malloc(sizeof(int)*content_size[buck_id]*2/100000);
    /* copy back */
    for(int i=0; i<content_size[buck_id]; i++)
    {
        tmp_cont[i] = content[buck_id][i];
    }
    for(int i=0; i<content_size[buck_id]/100000; i++)
    {
        tmp_c_idx[i] = c_idx[buck_id][i];
    }
    /* copy the pointer back */
    free(content[buck_id]);
    free(c_idx[buck_id]);
    content[buck_id] = tmp_cont;
    c_idx[buck_id] = tmp_c_idx;
    content_size[buck_id] = content_size[buck_id] * 2;
}

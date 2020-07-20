#include "list.h"
#include "insknap.h"

extern int parallel_mode; 
extern int which_med;
int g_count;
bool control;

#define FILTER 4

/**
 * The constructor for list.
 *
 * @param[in]       b_size          bucket size
 * @param[in]       list_size       size of each bucket
 * @param[in]       b               base index
 * @param[in]       num_t           number of threads
 * @param[in]       is_u            compute based on upper bound
 * @param[in]       is_enumerate_a  unumerate all possible solutions
 * @param[in]       log             place to keep the log file
 *
**/
List::List(int32_t b_size, int32_t list_size, int32_t b, int32_t num_t, bool is_u, bool *is_enumerate_a, const char* log)
{
    count            = new int*[num_threads];
	eliminate        = new bool*[num_threads];
	end              = new int*[num_threads];
    idx_sum          = new int*[num_threads];
	nei_score        = new int*[num_t];
	start            = new int*[num_threads];

    content_size     = new int64_t[b_size];
	elim_idx         = new int[num_t*CACHE_FILL];
    num              = new int64_t[b_size];

    base             = b;
    buck_size        = b_size; 
    is_ub            = is_u;
    is_enumerate_all = *is_enumerate_a;
	max_val          = 0;
	min_val          = 100000000;
	num_threads      = num_t;
	read_num         = 0;
	read_cnt         = 0;
	search_space     = 0;
	write_num        = 0;
	write_cnt        = 0;

    for(int i =0; i<num_threads; i++)
    {
        count[i]     = new int[b_size];
		eliminate[i] = new bool[LIST_SIZE*10]; 
		end[i]       = new int[b_size];
		idx_sum[i]   = new int[b_size+2];
		nei_score[i] = new int[LIST_SIZE*10]; 
		start[i]     = new int[b_size];
    }
	/* init values */
	for(int i=0; i<b_size; i++)
    {
		content_size[i] = list_size;
		num[i]          = 0;
		for(int j=0; j<num_t; j++)
        {
			count[j][i]   = 0;
			idx_sum[j][i] = 0;
			end[j][i]     = 0;
			start[j][i]   = 0;
		}
	}

	for(int i=0; i<10; i++)
    {
		time_elem e;
		t.list.push_back(e);
		t.init_timer(i);
	}
}

List::~List(){
}

bool
List::compute_partition(int32_t buck_id, Instance** ins){
//#ifdef USE_BW
	omp_lock_t writelock;
        omp_init_lock(&writelock);
//#endif
	int32_t i, j;
	//init counts
	for (i=0; i<num_threads; i++){
		elim_idx[i*CACHE_FILL] = 0;
		for(j=0; j<buck_size; j++){
			count[i][j] = 0;
			idx_sum[i][j] = 0;
		}
	}
	//parallel computing for partition
	bool is_lk_terminate=false;
#pragma omp parallel for num_threads(num_threads) private(i,j) 
	for (i=0; i<num[buck_id]; i++){
		int32_t tid = omp_get_thread_num();
		Instance *instance = ins[tid];
		int encode[30000];
		int number;
		get(i, buck_id, encode, &number, tid);
		
#ifdef USE_BW
		omp_set_lock(&writelock);
		read_num += number;
		read_cnt += 1;
		omp_unset_lock(&writelock);
#endif
		instance->to_encode(encode, number);	
		//printf("get from %d |", buck_id);
		//instance->print_encode();
		if(is_lk_terminate==true)
			continue; // find a better solution no need to continue search
		int possible_branch = instance->get_num_branches();
		//this part is for parallel lin-kernighan algorithm
		if(is_lk == true){
			instance->is_computing_bound=true;
			int start = elim_idx[tid*CACHE_FILL];
			int end = elim_idx[tid*CACHE_FILL]+possible_branch;
			instance->compute_bound();
			int o_score = instance->upper_bound;
			//for(int j=0; j<instance->o_size[3]; j++)
			for(j=0; j<possible_branch; j++){
				//omp_set_lock(&writelock);
				//instance->print_encode();
				//fflush(stdout);
				//omp_unset_lock(&writelock);
				instance->to_branch(j);
				//printf("%d %d\n", instance->upper_bound, upper_bound);
				if(instance->upper_bound < upper_bound){ //find the best solution
					//printf("parallel find a better %d<%d\n", instance->upper_bound , upper_bound);
					is_lk_terminate = true;
					//copy the best result
					omp_set_lock(&writelock);
					int encode[30000];
					int number = instance->get_encode(encode);
					ins[num_threads]->copy_code(encode, number);
					omp_unset_lock(&writelock);
				}
				else{
					if(is_lk_max==true)
						eliminate[tid][elim_idx[tid*CACHE_FILL]++] = true;
					else{
						if(instance->upper_bound<=o_score+FILTER){
							eliminate[tid][elim_idx[tid*CACHE_FILL]++] = false;
							//nei_score[tid][elim_idx[tid*CACHE_FILL]++];
							count[tid][buck_id+1]++;
							idx_sum[tid][buck_id+1] += instance->branch_offset;
						}
						else{
							eliminate[tid][elim_idx[tid*CACHE_FILL]++] = true;
						}
					}
				}
				instance->from_branch();
			}
			//decide the max one
			if(is_lk_max==true){
				int min = 1000000000;
				int min_id=-1;
				for(j=start; j<end; j++){
					if(nei_score[tid][j]<min){
						min = nei_score[tid][j];
						min_id = j;
					}
				}
				//eliminate[min_id] = false;
				instance->to_branch(j);
				//nei_score[tid][elim_idx[tid*CACHE_FILL]++];
				count[tid][instance->lower_bound]++;
				idx_sum[tid][instance->lower_bound] += instance->branch_offset;
				instance->from_branch();
			}
		}
		//this part is for branch and bound algoirhtm
		else if(is_lk==false){
			for(j=0; j<possible_branch; j++){
				instance->to_branch(j);
				if(is_ub==true){
					if(instance->upper_bound<lower_bound)
						eliminate[tid][elim_idx[tid*CACHE_FILL]++] = true;
					else{
						eliminate[tid][elim_idx[tid*CACHE_FILL]++] = false;
						count[tid][instance->upper_bound-base]++;
						//different problem has different branch_offset
						idx_sum[tid][instance->upper_bound-base] += (number + instance->branch_offset);
					}
#ifdef USE_BW
					omp_set_lock(&writelock);
					write_cnt += 1;
					omp_unset_lock(&writelock);
#endif
				}
				else{
					if(instance->lower_bound>upper_bound)
						eliminate[tid][elim_idx[tid*CACHE_FILL]++] = true;
					else{
						eliminate[tid][elim_idx[tid*CACHE_FILL]++] = false;
						count[tid][instance->lower_bound-base]++;
						//different problem has different branch_offset
						idx_sum[tid][instance->lower_bound-base] += (number + instance->branch_offset);
					}
#ifdef USE_BW
					omp_set_lock(&writelock);
					write_cnt += 1;
					omp_unset_lock(&writelock);
#endif
				}
				instance->from_branch();
			}
		}
	}
	//parallel list and sequential list are different 
	prepare_parallel_list(buck_id);
	if(is_lk_terminate == true)
		return false;
	else
		return true;
}

void
List::reorder_list(int32_t buck_id, Instance** ins){
	int i, j;
	omp_lock_t writelock;
	omp_init_lock(&writelock);
	for(int i=0; i<num_threads; i++)
		elim_idx[i*CACHE_FILL] = 0;
	//printf("num %d\n", num[buck_id]);
#pragma omp parallel for num_threads(num_threads) private(i,j) 
	for (i=0; i<num[buck_id]; i++){
		int32_t tid = omp_get_thread_num();
		Instance *instance = ins[tid];
		int encode[30000];
		int number;
		get(i, buck_id, encode, &number, tid);
#ifdef USE_BW
		omp_set_lock(&writelock);
		read_num += number;
		read_cnt += 1;
		omp_unset_lock(&writelock);
#endif
		//if(tid==0)
		//	t.tic_ser(0);
		instance->to_encode(encode, number);
		//if(tid==0)
		//	t.toc_ser(0);
		int possible_branch = instance->get_num_branches();
		//this part is for lin-kernighan algorithm
		if(is_lk == true){
			instance->is_computing_bound = false;
			for(j=0; j<possible_branch; j++){
				if(eliminate[tid][elim_idx[tid*CACHE_FILL]++]==true)
					continue;
				instance->to_branch(j);
				int pos = start[tid][buck_id+1]++;
				int encode[30000];
				int num_code = instance->get_encode(encode);
				add(pos, buck_id+1, num_code, encode, tid);	
				instance->from_branch();
			}	
		}
		else{
			for(j=0; j<possible_branch; j++){
				if(eliminate[tid][elim_idx[tid*CACHE_FILL]++]==true)
					continue;
				//if(tid==0)
				//	t.tic_ser(1);
				instance->to_branch(j);
				//if(tid==0)
				//	t.toc_ser(1);
				if(is_ub==true){
					//if(instance->upper_bound < lower_bound)
					//	continue;
					if(instance->lower_bound > lower_bound){
						//lock here
						omp_set_lock(&writelock);
						lower_bound = instance->lower_bound;
						omp_unset_lock(&writelock);
					}
					else if (instance->lower_bound >= upper_bound){
						omp_set_lock(&writelock);
						int encode[30000];
						int number = instance->get_encode(encode);
						ins[num_threads]->copy_code(encode, number);
						omp_unset_lock(&writelock);
					}
					if(instance->get_value()>max_val){
						//last round get the result
						omp_set_lock(&writelock);
						max_val = instance->get_value();
						int encode[30000];
						int number = instance->get_encode(encode);
						ins[num_threads]->copy_code(encode, number);
						omp_unset_lock(&writelock);
					}
					int pos = start[tid][instance->upper_bound-base]++;
					if(instance->upper_bound==buck_id+base){
						int encode[30000];
						int num_code = instance->get_encode(encode);
						//printf("thread %d add at buck %d at pos %d of size %d\n", tid, 0, pos, num_code);
						add(pos, 0, num_code, encode, tid);
#ifdef USE_BW
						omp_set_lock(&writelock);
						write_num += num_code;
						write_cnt += 1;
						omp_unset_lock(&writelock);
#endif
#ifdef USE_SPC
						omp_set_lock(&writelock);
						search_space++;
						omp_unset_lock(&writelock);
#endif
					}
					else{
						int bu_id = instance->upper_bound-base;
						int encode[30000];
						//if(tid==0)
						//	t.tic_ser(2);
						int num_code = instance->get_encode(encode);
						//if(tid==0)
						//	t.toc_ser(2);
						//printf("thread %d add at buck %d at pos %d of size %d\n", tid, bu_id, pos, num_code);
						add(pos, bu_id, num_code, encode, tid);
#ifdef USE_BW
						omp_set_lock(&writelock);
						write_num += num_code;
						write_cnt += 1;
						omp_unset_lock(&writelock);
#endif
#ifdef USE_SPC
						omp_set_lock(&writelock);
						search_space++;
						omp_unset_lock(&writelock);
#endif
					}
				}
				else{
					//if(instance->upper_bound < lower_bound)
					//	continue;
					if(instance->upper_bound < upper_bound){
						//lock here
						omp_set_lock(&writelock);
						upper_bound = instance->upper_bound;
						omp_unset_lock(&writelock);
					}
					else if (instance->upper_bound <= lower_bound){
						omp_set_lock(&writelock);
						int encode[30000];
						int number = instance->get_encode(encode);
						ins[num_threads]->copy_code(encode, number);
						omp_unset_lock(&writelock);
					}
					if(instance->get_value()<min_val){
						//last round get the result
						omp_set_lock(&writelock);
						min_val = instance->get_value();
						int encode[30000];
						int number = instance->get_encode(encode);
						ins[num_threads]->copy_code(encode, number);
						omp_unset_lock(&writelock);
					}
					int pos = start[tid][instance->lower_bound-base]++;
					if(instance->lower_bound==buck_id+base){
						int encode[30000];
						int num_code = instance->get_encode(encode);
						add(pos, buck_size-1, num_code, encode, tid);
#ifdef USE_BW
						omp_set_lock(&writelock);
						write_num += num_code;
						write_cnt += 1;
						omp_unset_lock(&writelock);
#endif
#ifdef USE_SPC
						omp_set_lock(&writelock);
						search_space++;
						omp_unset_lock(&writelock);
#endif
					}
					else{
						int bu_id = instance->lower_bound-base;
						int encode[30000];
						//if(tid==0)
						//	t.tic_ser(2);
						int num_code = instance->get_encode(encode);
						//if(tid==0)
						//	t.toc_ser(2);
						add(pos, bu_id, num_code, encode, tid);
						//printf("add ");
						//instance->print_encode();
						//if(number+instance->branch_offset != num_code)
						//	printf("wrong!\n");
						//printf("thread %d add at buck %d at pos %d of size %d\n", tid, bu_id, pos, num_code);
						//printf("number %d actual %d\n", number+instance->branch_offset, num_code);
#ifdef USE_BW
						omp_set_lock(&writelock);
						write_num += num_code;
						write_cnt += 1;
						omp_unset_lock(&writelock);
#endif
#ifdef USE_SPC
						omp_set_lock(&writelock);
						search_space++;
						//printf("search_space %d\n", search_space);
						omp_unset_lock(&writelock);
#endif
					}
				}
				//if(tid==0)
				//	t.tic_ser(3);
				instance->from_branch();
				//if(tid==0)
				//	t.toc_ser(3);
			}
		}
	}
	//copy temporary back
	int copy_size = copy_zero(buck_id);
#ifdef USE_BW
	write_num += num[buck_id];
	write_num += copy_size;
#endif
	//reset_numbers
	reset_num();
	//print_bucket(buck_id);
}

void
List::bnb_parallel_bucket(Instance **ins){
	is_parallel=true;
	if(is_ub==true){
		//sequential search first
		bnb(ins, 0);
		printf("finished sequential search!\n");
		while(upper_bound>lower_bound){
			int buck_id = upper_bound - base;
			if(num[buck_id]>0){
				//parallel version
				t.tic_ser(0);
				compute_partition(buck_id, ins);
				reorder_list(buck_id, ins);
				t.toc_ser(0);
			}
			else{
				upper_bound--;
			}
		}
	}
	else if(is_ub==false){
		//sequential search first
		bnb(ins, 0);
		printf("finished sequential bnb\n");
		while(upper_bound>lower_bound){
			int buck_id = lower_bound - base;
			if(num[buck_id]>0){
				//parallel version
				t.tic_ser(0);
				compute_partition(buck_id, ins);
				reorder_list(buck_id, ins);
				t.toc_ser(0);
			}
			else{
				lower_bound++;
				//printf("increase lower bound to %d\n", lower_bound);
			}
		}
	}
	ins[num_threads]->print_encode();
	//printf("time taken %f  %f %f %f %f %f \n", t.list[0].time, t.list[1].time, t.list[2].time, t.list[3].time, t.list[4].time, t.list[5].time);
}

void
List::bnb(Instance** ins, int ins_id){
	int tid = 0;
	upper_bound = ins[ins_id]->upper_bound;
	lower_bound = ins[ins_id]->lower_bound;
	if(is_ub==true)
    {
		int encode[30000];
		int number = ins[ins_id]->get_encode(encode);
        int buck_id = upper_bound -base;
		add(0, buck_id, number, encode, 0);
		num[buck_id]++;
	}
	else
    {
		int encode[30000];
		int number = ins[ins_id]->get_encode(encode);
		add(0, lower_bound-base, number, encode, 0);
		num[lower_bound-base]++;
	}
	int count =0;
	omp_lock_t writelock;
	omp_init_lock(&writelock);
	if(parallel_mode == MODE_SEQ)
		t.tic_ser(0);
	g_count =0;
	while(upper_bound > lower_bound){
		count++;
		/*******************for parallel********************************/
		if(parallel_mode == MODE_PAR){
			if(num[upper_bound - base]>parallel_base*num_threads)
				return;
		}
		else if(parallel_mode == MODE_TPAR){
			if(num[upper_bound - base]>parallel_base*num_threads){
				return;
			}
		}
		/******************end for parallel*********************************/
		if(is_ub==true)
        {
			int buck_id = upper_bound - base;
			while(num[buck_id]<=0 && upper_bound > lower_bound)
            {
                if(content_size[buck_id]>0)
                {
                    free_buck(buck_id);
                }
				buck_id--;
				upper_bound--;
				//printf("decrease to %d upper_bound %d lower_bound %d\n", buck_id+base, upper_bound, lower_bound);
			}
            if(num[buck_id] == 0)
            {
                //printf("come here the lower bound is %d upper bound is %d\n", upper_bound, lower_bound);
                ins[num_threads]->print_encode();
                return;
            }
			control=false;
			expand_ub(ins, tid, buck_id, writelock);
			//if(g_count % 100 == 0)
			//	printf("%d processed!\n", g_count);
            /**
             * TODO if g_count>=1000000
             * when processing gorilla_zebrafish
             * g_count=724866, there is a bug
             * */
			//if(g_count>=500)
			//	return;
			//print_buck();
			//exit(1);
			//printf("upper_bound %d, lower_bound %d\n", upper_bound, lower_bound);
		}
		else{
			int buck_id = lower_bound - base;
			while(num[buck_id]<=0 && upper_bound>lower_bound){
				buck_id++;
				lower_bound++;
				//printf("increase to %d\n", buck_id+base);
			}
			//printf("come to search %d with %lld\n", buck_id+base, num[buck_id]);
			expand_lb(ins, tid, buck_id, writelock);
		}

	}
	if(parallel_mode == MODE_SEQ)
		t.toc_ser(0);

    ins[num_threads]->print_encode();
    /* finalize bnb */
	//printf("search space %d\n", g_count);
}



void
List::expand_ub(Instance **ins, int tid, 
		int buck_id, omp_lock_t writelock){
	int j;
	int encode[30000];
	int number;
	Instance *instance = ins[tid];
	//print_buck();
	get(num[buck_id]-1, buck_id, encode, &number, 0);
	instance->to_encode(encode, number);
	num[buck_id]--;
#ifdef USE_BW
	read_num+=number;
	read_cnt+=1;
#endif
	int possible_branch = instance->get_num_branches();
	//printf("num possible branch %d \n", possible_branch);
	for(j=0; j<possible_branch; j++){
		g_count++;
		//printf("j %d upper_bound %d lower_bound %d g_count %d\n", j, instance->upper_bound, instance->lower_bound, g_count);
		instance->to_branch(j);
		if(instance->upper_bound < lower_bound)
			;
		else if(instance->upper_bound > upper_bound)
			instance->upper_bound = upper_bound;
		else{
			if(instance->lower_bound > lower_bound){
				if(is_enumerate_all==true)
					; //do nothing there
				else //change lower bound
					lower_bound = instance->lower_bound;
			}
			if(instance->lower_bound >= upper_bound){
				if(parallel_mode == MODE_TPAR){
					omp_set_lock(&writelock);
				}
				int encode[30000];
				int number = instance->get_encode(encode);
				ins[num_threads]->copy_code(encode, number);
				//ins[num_threads]->print_encode();
				//printf("best value is %d\n", lower_bound);
				if(parallel_mode == MODE_TPAR){
					omp_unset_lock(&writelock);
				}
				return;
			}
#ifdef USE_BW
			write_num+=number;
			write_cnt+=1;
#endif
			if(instance->upper_bound>lower_bound){
				int encode[30000];
				int number = instance->get_encode(encode);
				int buck_id = instance->upper_bound-base;
				int pos = num[buck_id];
				add(pos, buck_id, number, encode, 0);
				num[buck_id]++;
			}
			if(instance->get_value() > max_val){
				if(parallel_mode == MODE_TPAR){
					omp_set_lock(&writelock);
				}
				int encode[30000];
				int size = instance->get_encode(encode);
				ins[num_threads]->copy_code(encode, size);
				max_val = instance->get_value();
				if(parallel_mode == MODE_TPAR){
					omp_unset_lock(&writelock);
				}
				//printf("value %d\n", max_lb);
			}
#ifdef USE_SPC
			search_space++;
#endif
		}
		instance->from_branch();
	}
}

void
List::expand_lb(Instance **ins, int tid, 
		int buck_id, omp_lock_t writelock){
	int j;
	int encode[30000];
	int number;
	Instance *instance = ins[tid];
	get(num[buck_id]-1, buck_id, encode, &number, 0);
	instance->to_encode(encode, number);
	num[buck_id]--;
#ifdef USE_BW
	read_num+=number;
	read_cnt+=1;
#endif
	int possible_branch = instance->get_num_branches();
	for(j=0; j<possible_branch; j++){
		instance->to_branch(j);
		//printf("computing branch %d lower_bound %d upper_bound %d\n", j, instance->lower_bound, instance->upper_bound);
		if(instance->lower_bound > upper_bound)
			;
		else if(instance->lower_bound < lower_bound)
			instance->lower_bound = lower_bound;
		else{
			if(instance->upper_bound < upper_bound){
				upper_bound = instance->upper_bound;
				//printf("decreasing upper bound to %d\n", upper_bound);
			}
			if(instance->upper_bound <= lower_bound){
				if(parallel_mode == MODE_TPAR){
					omp_set_lock(&writelock);
				}
				int encode[30000];
				int number = instance->get_encode(encode);
				ins[num_threads]->copy_code(encode, number);
				ins[num_threads]->upper_bound = instance->upper_bound;
				ins[num_threads]->lower_bound = instance->lower_bound;
				//ins[num_threads]->print_encode();
				//printf("best value is %d\n", upper_bound);
				if(parallel_mode == MODE_TPAR){
					omp_unset_lock(&writelock);
				}
				return;
			}
			int buck_id = instance->lower_bound-base;
			int pos = num[buck_id];
			int encode[30000];
			int number = instance->get_encode(encode);
			//printf("pos %d buck_id %d, number %d\n", pos, buck_id, number);
			//printf("add to %d lower %d g upper %d\n", buck_id+base, instance->lower_bound, upper_bound);
#ifdef USE_BW
			write_num+=number;
			write_cnt+=1;
#endif
			add(pos, buck_id, number, encode, 0);
			num[buck_id]++;
			if(instance->get_value() < min_val){
				if(parallel_mode == MODE_TPAR){
					omp_set_lock(&writelock);
				}
				int encode[30000];
				int size = instance->get_encode(encode);
				ins[num_threads]->copy_code(encode, size);
				ins[num_threads]->upper_bound = instance->upper_bound;
				ins[num_threads]->lower_bound = instance->lower_bound;
				min_val = instance->get_value();
				if(parallel_mode == MODE_TPAR){
					omp_unset_lock(&writelock);
				}
				//printf("value %d\n", max_lb);
			}
#ifdef USE_SPC
			search_space++;
#endif
		}
		instance->from_branch();
	}
}

/**
 * this is the implementation of lin-kernighan algorithm
 * which is parallelizable 
 *
 * @param[in]       ins             instances
 * @param[in]       heu_level       heuristic level to run
 * @param[in]       term_move       which level to stop running             
 * @param[in]       is_opt          is using k-opt method or not
**/
void
List::lk(Instance** ins, int heu_level, 
		int term_move, bool is_opt){
	/* init the list */
	bool improved=true;
	bool accept=false;	
	int best_score = ins[0]->upper_bound;
	lk_terminate = false;
	is_lk=true;
	/* upperbound(best_score) is the current score */
	/* lowerbound is (d_12 + d_23 + d_13)/2 */
	int use_heu = true;
	while(improved==true)
    {
		int current_level=0;
		accept=false;
		count = 0;
		improved = false;
		bool expansion_result = false;
		while(accept==false && current_level<term_move)
        {
			if(is_opt==false && current_level < heu_level)
            {
				expansion_result = add_child_combinations(ins, current_level, use_heu);
			}
			else if(is_opt==false && current_level >= heu_level)
            {
				expansion_result = add_child_max(ins, current_level, term_move);
			}
			/* the following part is for k-opt algorithm */
			else if(is_opt==true && current_level < heu_level)
            {
				expansion_result = add_child_combinations(ins, current_level, use_heu);
			}
			else if(is_opt==true && current_level >= heu_level)
            {
				;
			}
			current_level += 1;
			
			if(expansion_result == true)
            {
				accept=true;
				best_score = upper_bound;
				improved = true;
				reset_list();
				/* copy code to instance 0 */
				int encode[30000];
				int number = ins[num_threads]->get_encode(encode);
				ins[0]->copy_code(encode, number);
			}
		}
	}
}



bool 
List::add_child_combinations(Instance **ins, int current_level, int use_heu)
{
	bool result=false;
	if(current_level == 0  || false == is_parallel)
    {
		ins[0]->is_computing_bound=true;
		int possible_branch = ins[0]->get_num_branches();
		for(int j=0; j<possible_branch; j++)
        {
            printf("to %dth branch!\n", j);
			ins[0]->to_branch(j);
            exit(1);
            /* if find a better result */
			if(ins[0]->upper_bound < upper_bound)
            {
				upper_bound = ins[0]->upper_bound;
				int encode[30000];
				int number = ins[0]->get_encode(encode);
				ins[num_threads]->copy_code(encode, number);
				result = true;
				break;
			}
			int pos = j;
			int encode[30000];
			int num_code = ins[0]->get_encode(encode);
			if(ins[0]->upper_bound <= upper_bound+FILTER)
            {
				add(pos, current_level+1, num_code, encode, 0);	
				num[current_level+1]++;
			}
			ins[0]->from_branch();
            exit(1);
		}
	}
	else
    {
		//
		is_lk_max = false;
		printf("come to partition!\n");
		result = compute_partition(current_level, ins);
		printf("finished partition!\n");
		if (result == true)
			reorder_list(current_level, ins);
		if(result==true) result=false; else result = true;
	}
	return result;
}


bool 
List::add_child_max(Instance **ins, int current_level, int use_heu){
    bool result = false;
    if(false == is_parallel)
    {
		ins[0]->is_computing_bound=true;
		int possible_branch = ins[0]->get_num_branches();
        int min_bound = 10000000;
        int min_pos =-1;
		for(int j=0; j<possible_branch; j++)
        {
			ins[0]->to_branch(j);
            /* if find a better result */
			if(ins[0]->upper_bound < upper_bound)
            {
				upper_bound = ins[0]->upper_bound;
				int encode[30000];
				int number = ins[0]->get_encode(encode);
				ins[num_threads]->copy_code(encode, number);
				result = true;
				break;
			}
            if(ins[0]->upper_bound < min_bound)
            {
                min_bound = ins[0]->upper_bound;
                min_pos = j;
            }
			ins[0]->from_branch();
		}
        /* perform insert */
        int pos = min_pos;
        int encode[30000];
        ins[0]->to_branch(pos);
        int num_code = ins[0]->get_encode(encode);
        if(ins[0]->upper_bound <= upper_bound+FILTER)
        {
            add(pos, current_level+1, num_code, encode, 0);	
            num[current_level+1]++;
        }
        ins[0]->from_branch();
    }
    else
    {
        printf("come to recursive max!\n");
        is_lk_max = true;
        result = compute_partition(current_level, ins);
        printf("finished partition!\n");
        if (result == true)
            reorder_list(current_level, ins);
        if(result==true) result=false; else result = true;
    }
    return result;

}


void 
List::reset_list(){
    for(int i=0; i<buck_size; i++)
        num[i] = 0;
}


//fclose(dbg);
//printf("the best score is %d\n", best_score);
void
free_list(){
}

void 
List::add(int64_t pos, int32_t buck_id, 
        int num_code, int *encode, int ins_id){
}

void 
List::get(int64_t pos, int32_t buck_id, 
        int *encode, int *number, int ins_id){
}

void
List::prepare_parallel_list(int buck_id){
}

int 
List::copy_zero(int buck_id){
    return -1;
}

void
List::reset_num(){
}

void 
List::copy(List *other, int buck_id, int size){
}

void
List::print_bucket(int buck_id){
}

void
List::print_buck(){
    for(int i=lower_bound; i<=upper_bound; i++)
        if(num[i-base] != 0)
            printf("%d(%d): %lld\n", i, i-base, num[i-base]);
}

void 
List::add_g_count_list(int buck_id, int g_count){
#ifdef USE_DEBUG
    g_count_list[buck_id][g_count_idx[buck_id]++] = g_count;
#endif
}

void List::free_buck(int buck_id)
{
}

void List::reallocate_buck(int buck_id)
{
}

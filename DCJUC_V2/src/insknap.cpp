#include "insknap.h"
#include "utils.h"
#include <cmath>


InsKnapsack::InsKnapsack(const char *file, int id, const char *lf) : 
    Instance(lf)
{
    int i;
    FILE *reader = fopen(file, "r");
    if(reader==NULL){
        printf("file %s not exits!\n", file);
        exit(1);
    }
    int total_item, space_limit; 
    fscanf(reader, "%d %d\n", &total_item, &space_limit);
    items_space = (int*)malloc(sizeof(int)*total_item);
    items_value = (int*)malloc(sizeof(int)*total_item);
    branch_id = (int*)malloc(sizeof(int)*total_item);
    sorted = (int*)malloc(sizeof(int)*total_item);
    map = (int*)malloc(sizeof(int)*total_item);
    in_bag = (bool*)malloc(sizeof(bool)*total_item);
    for(i=0; i<total_item; i++){
        fscanf(reader, "%d %d\n", &items_space[i], &items_value[i]);
        in_bag[i] = false;    
        sorted[i] = (int)ceil((double)items_value[i]/(double)items_space[i]);
        map[i] = i;
    }
    Utils::q_sort_two(sorted, map, 0, total_item-1);
    total_space = space_limit;
    value = 0;
    used_space = 0;
    num_item = total_item;
    num_count = total_item;
    chosen_now = -1;
    num_branch = 0;
    compute_min_max();
    compute_bound();
    tid = id;
    time_elem e;
    t.list.push_back(e);
    t.init_timer(0);
}

InsKnapsack::InsKnapsack(const InsKnapsack &other) : 
    Instance(other)
{
    int i;
    num_item = other.num_item;
    items_space = (int*)malloc(sizeof(int)*num_item);
    items_value = (int*)malloc(sizeof(int)*num_item);
    branch_id = (int*)malloc(sizeof(int)*num_item);
    sorted = (int*)malloc(sizeof(int)*num_item);
    map = (int*)malloc(sizeof(int)*num_item);
    in_bag = (bool*)malloc(sizeof(bool)*num_item);
    for (i=0; i<num_item; i++){
        items_value[i] = other.items_value[i];
        items_space[i] = other.items_space[i];
        in_bag[i] = false;      
        sorted[i] = (int)((double)items_value[i]/(double)items_space[i]);
        map[i] = i;
    }   
    Utils::q_sort_two(sorted, map, 0, num_item-1);
    total_space = other.total_space;
    value = 0;
    used_space = 0;
    num_count = num_item;
    chosen_now = -1; 
    num_branch = 0;
    max = other.max;
    min = other.min;
}

void 
InsKnapsack::to_branch(int which_branch){
    //bool enable = false;
    //if(used_space==656 && value==1277){
    //    enable = true;
    //}
    int b_id = branch_id[which_branch];
    in_bag[b_id] = true;
    used_space += items_space[b_id];
    value += items_value[b_id];
    chosen_now = b_id;
    if(chosen_now==max_now || chosen_now==min_now)
        compute_min_max();
    compute_bound();
    //if(enable ==true){
    //print_encode();
    //printf("used_space %d value %d\n", used_space, value);
    //}
    //printf("branch: %d space %d value %d max_now %d min_now %d max %3.1f min %3.1f ", 
    //    chosen_now, used_space,  
    //    value, 
    //    max_now, min_now, max, min);
}

void 
InsKnapsack::from_branch(){
    in_bag[chosen_now] = false;
    used_space -= items_space[chosen_now];
    value -= items_value[chosen_now];
    double val = (double)items_value[chosen_now]/(double)items_space[chosen_now];
    if(val>max || val<min){
        compute_min_max();
    }
    chosen_now = -1;
}

void 
InsKnapsack::compute_bound(){
    int i;
    int space_left = total_space-used_space;
    upper_bound = value + space_left*max;
    int tmp_val=value;
    int tmp_spc=space_left;
    //if(control==true){
    //    printf("space %d val %d\n", tmp_spc, tmp_val);
    //    int t_v=0;
    //    int t_s=0;
    //    for(i=0; i<num_item; i++)
    //        if(in_bag[i] ==true){
    //            t_v += items_value[i];
    //            t_s += items_space[i];
    //            printf("%d %d\n", t_s, t_v);
    //        }
    //}
    bool tmp_in_bag[100];
    for (i=0; i<num_item; i++)
        tmp_in_bag[i] = in_bag[i];
    for(i=num_item-1; i>=0; i--){
        int space = items_space[map[i]];    
        if(tmp_spc-space <0 || tmp_in_bag[map[i]]==true)
            continue;
        tmp_spc-=space;
        tmp_val+=items_value[map[i]];
        tmp_in_bag[map[i]] = true;
        if(control==true)
            printf("space %d value %d\n", tmp_spc, tmp_val);
    }
    lower_bound = tmp_val;
}

int 
InsKnapsack::get_num_branches(){
    int i;
    int result =0;
    for(i=0; i<num_item; i++){
        if(in_bag[i]==false && items_space[i]+used_space <= total_space){
            branch_id[result++] = i;
        }
    }
    num_branch = result;
    branch_offset = 1;
    return result;
}

void
InsKnapsack::compute_min_max(){
    int i;
    max =0;
    for(i=0; i<num_item; i++){
        if(in_bag[map[i]]==true)
            continue;
        int val = sorted[i];
        if(val>max){
            max = val;
            max_now = map[i];
        }
    }
    //find min and compute lower_bound
    min =100000000;
    for(i=0; i<num_item; i++){
        if(in_bag[map[i]]==true)
            continue;
        int val = sorted[i];
        if(val<min){
            min = val;
            min_now = map[i];
        }
    }
}


int  
InsKnapsack::get_encode(int* encode){
    int i, j;
    //printf("encode: ");
    //if (tid==0)
    //    t.tic_ser(0);
    for(i=0,j=0; i<num_item; i++){
        if(in_bag[i]==true){
            encode[j++] = i;
            //printf("%d ", i);
        }
    }
    //if (tid==0)
    //    t.toc_ser(0);
    //printf("\n");
    return j;
}

void
InsKnapsack::print_encode(){
    int i, j;
    for(i=0,j=0; i<num_item; i++){
        if(in_bag[i]==true){
            printf("1 ");
        }
        else 
            printf("0 ");
    }
    printf("\n");
}

void 
InsKnapsack::to_encode(int* encode, int size){
    int i;
    used_space = 0;    
    value = 0;
    for(i=0; i<num_item; i++){
        in_bag[i] = false;
    }
    for(i=0; i<size; i++){
        //printf("%d ", encode[i]);
        in_bag[encode[i]] = true;
        used_space += items_space[encode[i]];
        value += items_value[encode[i]];
        //printf("%d:%d ", used_space,value);
    }
    //for(i=0; i<num_item; i++)
    //    printf("%s ", in_bag[i]?"1":"0");
    //printf("space %d value %d\n", used_space, value);
    //printf("\n");
    compute_min_max();
}

void 
InsKnapsack::copy_code(int *encode, int size){
    int i;
    to_encode(encode, size);
    int space_left = total_space-used_space;
    int tmp_val=value;
    int tmp_spc=space_left;
    for(i=num_item-1; i>=0; i--){
        int space = items_space[map[i]];
        if(tmp_spc-space <0 || in_bag[i]==true)
            continue;
        tmp_spc-=space;
        tmp_val+=items_value[map[i]];
        in_bag[map[i]] = true;
    }
}


int
InsKnapsack::get_value(){
    return value;
}

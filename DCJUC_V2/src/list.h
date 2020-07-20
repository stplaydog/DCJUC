#pragma once


#ifndef _LIST_H
#define _LIST_H

#include <cstdio>
#include <cstdlib>
#include "omp.h"
#include "instance.h"
#include <vector>
#include <assert.h>
//#include "timer.h"
#define MODE_SEQ 1
#define MODE_PAR 2
#define MODE_TPAR 3
#define LIST_SIZE 1000000
#define CACHE_FILL 64

using namespace std;

class List{
public:
    time_list t;           ///<

    int32_t **count;       ///<
    int32_t **end;         ///<
    int32_t **idx_sum;     ///<
    int32_t **l_lb;        ///<
    int32_t **l_ub;        ///<
    int32_t **nei_score;   ///<
    int32_t **start;       ///<

    bool    **eliminate;   ///<

    int64_t *content_size; ///<
    int64_t *num;          ///<

    int32_t *elim_idx;     ///<

    int64_t read_num;      ///<
    int64_t write_num;     ///<
    int64_t read_cnt;      ///<
    int64_t write_cnt;     ///<

    int32_t upper_bound;   ///<
    int32_t lower_bound;   ///<
    int32_t num_threads;   ///<
    int32_t buck_size;     ///< 
    int32_t base;          ///<
    int32_t max_val;       ///<
    int32_t min_val;       ///<
    int32_t search_space;  ///<
    int32_t parallel_base; ///<

    bool is_ub;            ///<
    bool is_enumerate_all; ///<
    bool is_parallel;      ///<
    bool lk_terminate;     ///<
    bool is_lk_max;        ///<
    bool is_lk;            ///<

    FILE *logger;

    //this for debug purpose only
#ifdef USE_DEBUG
    int32_t **g_count_list;
    int32_t *g_count_idx;
#endif

public:
    /* main functions */
    List(int32_t buck_size, int32_t list_size, int32_t base, int32_t num_t, bool is_ub, bool *is_enumerate_all = NULL, const char* log = NULL);
    virtual ~List();

    bool compute_partition(int32_t buck_id, Instance** ins);
    void reorder_list(int32_t buck_id, Instance** ins);
    void bnb_parallel_bucket(Instance** ins);
    void bnb(Instance** ins, int32_t ins_id);
    void expand_ub(Instance **ins, int32_t tid, int32_t buck_id, omp_lock_t writelock);
    void expand_lb(Instance **ins, int32_t tid, int32_t buck_id, omp_lock_t writelock);
    void print_buck();

    //lin kernighan algorithm implementation
    void lk(Instance** ins, int32_t heu_level, int32_t term_move, bool is_opt);
    bool add_child_combinations(Instance **ins, int32_t current_level, int32_t use_heu);
    bool add_child_max(Instance **ins, int32_t current_level, int32_t use_heu);
    void reset_list();
    void add_g_count_list(int32_t pos, int32_t g_count);

    //virtual functions
    virtual void add(int64_t pos, int32_t buck_id, int32_t num_code, int32_t *encode, int32_t ins_id);
    virtual void get(int64_t pos, int32_t buck_id, int32_t *encode, int32_t *num, int32_t ins_id);
    virtual void prepare_parallel_list(int32_t buck_id);
    virtual int32_t copy_zero(int32_t buck_id);
    virtual void reset_num();
    virtual void copy(List *other, int32_t buck_id, int32_t size);
    virtual void print_bucket(int32_t buck_id);
    virtual void free_buck(int32_t buck_id);
    virtual void reallocate_buck(int32_t buck_id);
};

#endif

#pragma once
#ifndef _H_INSTANCE_
#define _H_INSTANCE_

#include "timer.h"
#include "utils.h"
#include <assert.h>

using namespace std;

class Instance{
public:
    int upper_bound;
    int lower_bound;
    int num_count;
    bool control;
    int tid;
    time_list t;
    int branch_offset;
    bool is_computing_bound;
    int score;
    FILE *logger;
public:
    Instance();
    Instance(const Instance &other);
    Instance(const char* lf = NULL);
    virtual ~Instance();
    virtual void to_branch(int which_branch);
    virtual void from_branch();
    virtual void compute_bound();
    virtual int get_num_branches();
    virtual int get_encode(int *encode);
    virtual void to_encode(int *encode, int size);
    virtual void print_encode();
    virtual void copy_code(int *encode, int size);
    virtual int get_value();
    virtual int return_private_info();
};


#endif

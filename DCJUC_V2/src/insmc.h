/* Copyright (c) 2015, YMSys and/or its affiliates. All rights reserved.*/

/**
  @file     insmc.cpp

  @brief    Instance for solving maximum clique problem.

  MODIFIED   (MM/DD/YY)
  zyin    07/15/15 - creation
 **/


#pragma once
#ifndef _H_INSMC
#define _H_INSMC

#include "instance.h"
#include <cstdio>
#include <cstdlib>

class InsMC : public Instance {
public:
    /* where to read the graph */
    char *in_file;
    char *CC_folder;

    //***this is for the storage of the graph
    int v_size;
    int e_size;
    int *v_idx;
    int *e_idx;

    //***this is specific to the implemetation of Cing algorithm
    //Q 
    int *Q; //Q should be sorted
    int Q_size; //Q size
    //R
    int *R; 
    int *tmp_R;
    int *R_deg;
    int R_size;
    bool *R_avail;
    //C
    int *C_map;
    int max_C;

    //***these are for branching functions
    int branch_id;
    int *branch_v; //which vertex in this branch
    int num_branch;

public:
    //***constructors and destructors
    InsMC(const char *f = NULL, const char *od = NULL, const char *lf = NULL);
    ~InsMC();

    //***branch and bound functions
    void to_branch(int which_branch);
    void from_branch();
    void compute_bound();
    int get_num_branches();

    //***storage related functions
    int get_encode(int *encode);
    void to_encode(int *encode, int size);
    void print_encode();
    void copy_code(int *encode, int size);
    int get_value();
    int return_private_info();
    void get_num_elem();
private:
    //for initialization of graph
    void init_graph();
    void separate_CC();
    //Three function need to be implemented for the basic BnB algo
    void join_R(int v);
    void get_deg_sorted();
    void assign_C();
    //something for debugging
    void vis_csr();
    void print_deg();
    void get_num_triangles();
    void allocate_data_structure(const int v_num, 
            const int e_num);

    /**
     * initialize the full R
     *
     * @return      N/A
     **/
    void get_full_R()
    {
        for(int i=0; i<v_size; i++)
        {
            R[i] = i;
        }
        /* All the vertices are in R */
        R_size = v_size;
    }

    /**
     * Find the first available vertex. 
     *
     * @return      the first available vertex. 
     **/
    int find_first_avail_pos()
    {
        int first = OPTKIT_NULL;
        for(int i=0; i<R_size; i++) 
        {
            if(C_map[tmp_R[i]] == OPTKIT_NULL)
            {
                first = i;
                break;
            }
        }
        return first;
    }
};

#endif

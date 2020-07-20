#pragma once
#ifndef _H_INSDIS
#define _H_INSDIS

#include "instance.h"
#include <cstdio>
#include <cstdlib>

#define MAX_VET_SIZE 32768
#define PI_OPEN 1
#define GAMMA_OPEN 2
#define CYCLE 3
#define END 4
#define NONE 5
#define OUT 6
#define CAP 65533
#define CUP 65534


/**
 * @class InsDis
 *
 * This class calculate the DCJ distance between two genomes
 * with unequal contents.
 * TODO bijection calculation is still not correct
 * TODO turn define into enum
 * TODO decide which function tobe inlined
 * TODO after passing the regression test, remove these debug code
 *
**/
class InsDis : public Instance {
public:
    char *file;         ///< where to read the graph file
    bool is_calculate_bijection;
                        ///< whether to calculate bijection or not 
    char *bij_file;     ///< where to read bijection result
    int **bijection;    ///< data structure to store bijection result
    int *bij_idx;       ///< index of bijection result
    char* dict_file;    ///< because graph is shrinked, we need to have a dictionary
    int *dict;          ///< dictionary to keep 
    int v_size;         ///< num vertices, = 2*gf_size
    int e_size[2];      ///< num edges in CSR
    int **v_idx;        ///< index of CSR
    int **e_idx;        ///< real edge of CSR
    int ***adj;         ///< adjacency format of graph
    int **adj_idx;      ///< index of the adj
    int ***a_tmp;       ///< adj format for branching only
    int **idx_tmp;      ///< index of a_tmp
    int **a;            ///< adj for the real computation
    int *vet_type;      ///< vertex type, see compeau's paper
    int a_size;         ///< size of the adj
    vector<vector <int > > comb;
                        ///< to store combinations for branching
    int avail_id;       ///< the next available vertex id (for matc only)
    bool is_ub;         ///< is computing upper bound or not
    bool is_exem;       ///< is computing exemplar distance or not
    int vet;            ///< which vetex is now branching
    int branch_id;      ///< which branch in the comb
    int num_branch;     ///< number of possible branches for the current vet
    int num_v_after_relabel;
                        ///< how many vertices there are after relabel 
    int encode_size;    ///< how many data to store the encode
    int comb_size;      ///< size of the comb
    int pair_num;       ///< how many pairs there are
    int count;          ///< total count

    /* constructors and destructors */
    InsDis(const char *f, bool is_exe, const char *lf = NULL, const char* fd = NULL, const char* fbi = NULL, const bool *is_cal_bij = NULL); 
    ~InsDis();
    void to_branch(int which_branch);
    void from_branch();
    void compute_bound();
    int get_num_branches();
    int get_encode(int *encode);
    void to_encode(int *encode, int size);
    void print_encode();
    void copy_code(int *encode, int size);
    int get_value();
    int return_private_info();
private:
    void init_graph();
    void scan_graph();
    void read_graph();
    void from_csr_to_adj();
    int compute_indel_dis(); 
    void get_exemplar_adj(int vet, int comb_id);
    void get_matching_adj(int vet, int comb_id);
    void get_deletion_adj(int vet, int comb_id, bool is_exem);
    void get_comb_exem();
    void get_comb_matc();
    bool check_if_in_to(int v, int c);
    void set_vet_type(int size);
    void set_visited(int *vet_type, int pos);
    void graph_vis(char *file);
    void print_v();
    void print_adj();
    void print_tmp();
    void get_num_elem();
    void print_comb();
    void vis_adj();
    void vis_a(int size);
    void vis_tmp();
    void vis_csr();
    bool detect_no_branch(int *encode);
    void read_dict();
    int append_bijection_encode(int* encode, int start);
    void get_bijection_encode(int* encode, int start);
    bool detect_abnornal();
};

#endif

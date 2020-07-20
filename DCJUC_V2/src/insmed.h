#pragma once

#ifndef _H_INSMED
#define _H_INSMED

#include "instance.h"
#include "list.h"
#include "listint.h"
#include "insdis.h"
#include <cstdio>
#include <cstdlib>


#define CAP 65533

/**
 * @class   InsMed      instance used to calculate median genome.
 *
 * TODO     add zero matching
 * TODO     turn order into graph
**/
class InsMed : public Instance{
public:
    enum conn_type { CONN_TYPE_ONE, CONN_TYPE_TWO }; 
    enum comp_type { TYPE_SING, TYPE_MUL }; 
    enum dis_type  { DIS_MODE_EXEM, DIS_MODE_MATC};
    int **v_idx;            ///< v_idx for CSR
    int **e_idx;            ///< e_idx for CSR
    int ***adj;             ///< Adjacency graph 
    int **a_idx;            ///< Adjacency graph index
    int **comp_id;          ///< Component id
    int **comp_type;        ///< Component type
    int *med;               ///< Store median gene order
    char *m_file;           ///< File to read graph
    char *m_tmp_folder;     ///< Folder to store tmp graph
    int branch_selected;    ///< Indicate which branch is selected
    int num_t;              ///< How many number of threads there are 
    int score;              ///< The score we get
    int gf_size;            ///< The size of the gene family
    int v_size;             ///< # vertices
    int e_size[4];          ///< size of edges
    int cycle;              ///< # cycles we already know
    vector<vector<int> > ft;///< Footprints of possible combinations
    vector<vector<int> > signature;   
                            ///< Current status
    List **l_d;             ///< list data structure
    Instance ***ins;        ///< instance data structure
    int dis_mode;           ///< Distance mode, exemplar:1 matching:2
    bool use_heu;           ///< Use heuristic to reduce search space or not
    int thresh;             ///< heuristic threshold
    bool is_branching;      ///< when compute_bound, turn this on or off matters

    InsMed(const char *file, const char* tmp_folder, int tid, const char *lf = NULL, const int *dm = NULL, const bool *uh = NULL, const int *th = NULL);
    InsMed(const InsMed &other);
    void to_branch(int which_branch);
    void from_branch();
    void compute_bound();
    int get_num_branches();
    int get_encode(int* encode);
    void to_encode(int* encode, int size);
    void print_encode();
    void copy_code(int *encode, int size);
    int get_value();
    void compute_num_elem();

private:
    void csr_to_file(char *file, int c1, int c2);
    void init_dis_tools();
    void init_graph();
    void scan_graph();
    void read_graph();
    void from_csr_to_adj(int c);
    void rename_graph_by_car();
    void from_adj_to_csr();
    void init_median();
    int find_next_edge(int *rev_v, int from);
    bool edge_exist(int from, int to);
    int get_deg(int vet);
    int max_deg(int vet);
    void indentify_comps();
    void shrink(int *encode);
    void resume(int *encode);
    void add_possible_branch(int from_one, int to_one, int from_two, int to_two);
    int detect_comps(int *c_type, int c1, int c2);
    void bfs_visit(int *stack, bool *visited, 
            int *c_type, int vet, int type_id, int &s_idx);
    int get_encode_by_component(int *c_type, int *encode, 
            int which_comp, int c1, int c2);
    void marriage_cap_vet();
    void search_cap_edge(int one_vet, int &ret_edge, int &ret_vet);
    bool detect_n_remove_conflict_edge();
    bool check_if_single_cc(int *encode);
    void rename_encode_by_comp(int *comp);
    void print_c_type(int *c_type);
    /* TODO These are debug functions, after development, they should be cleaned up*/
    void print_csr();
    void print_adj();
    void vis_csr();
    void vis_adj();
    void vis_two(int c1, int c2);
    void vis_encode(int *encode);
}; 

#endif

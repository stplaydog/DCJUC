/* Copyright (c) 2015, YMSys and/or its affiliates. All rights reserved.*/

/**
  @file     insmc.cpp

  @brief    Instance for solving maximum clique problem.
 
   TODO Change all neccessary const arguments.
   TODO Change all neccessary const functions. 
   TODO printf to fprintf STDOUT
   TODO Separate debug functions into util

   MODIFIED   (MM/DD/YY)
       zyin    07/10/15 - Test initial graph and assign color using janc.gr
       zyin    07/07/15 - Some basic formatting and the test of assign_C 
**/


#include "insmc.h"
#include "utils.h"
#include <string>
#include <iostream>
#include <cstring>
#include <string>

using namespace std;

#define DEG_VAL 0
#define DEG_IDX 1


/**
 * Constructor 
 *
 * param[in]        f       the input file
 * param[in]        od      the output directory
 * param[in]        lf      the logger file 
 *
 * @return      N/A
 **/
InsMC::InsMC(const char *f, const char *od, const char *lf) :
    Instance(lf)
{
    /* Some file related information */
    in_file   = NULL;
    CC_folder = NULL;

    if(f != NULL)
    {
        in_file = new char[OPTKIT_FILE_SIZE]; 
        snprintf(in_file, OPTKIT_FILE_SIZE, "%s", f);
    }

    if(od != NULL)
    {
        CC_folder = new char[OPTKIT_FILE_SIZE];
        snprintf(CC_folder, OPTKIT_FILE_SIZE, "%s", od);
    }

    /* Some bound related information */
    upper_bound = OPTKIT_NULL;
    lower_bound = OPTKIT_NULL;
    branch_id = OPTKIT_NULL;
    num_branch = OPTKIT_NULL;

    /* init graph and init some related data structures */
    init_graph();
    printf("finished reading graph!\n");

    /* here if we have the folder information just perform the CC algorithm */
    if(CC_folder != NULL)
    {
        separate_CC();
        get_num_triangles();
    }
    else
    {
        /* get initial R, and perform the color assignment function */
        get_full_R(); 
        get_deg_sorted();
        assign_C();

        /* compute upper and lower bound now */
        upper_bound = max_C;
        /* this is important, smaller than this will not be considered */
        lower_bound = 2; 

        printf("initial upper bound %d lower bound %d\n", upper_bound, lower_bound);
        
        /* This code tells us about the stack size. */
        get_num_elem();

        if(logger != NULL)
        {
            fprintf(logger, "====BOUNDS====\n");
            fprintf(logger, "Upper bound: %d Lower bound: %d\n", upper_bound, lower_bound);
        }
    }
}

/**
 * Destructor to free all the things
 *
 * @return      N/A
**/
InsMC::~InsMC()
{
    free(in_file);
    free(v_idx);
    free(e_idx);
    free(Q);
    free(R);
    free(tmp_R);
    free(R_deg);
    free(R_avail);
    free(C_map);
    free(CC_folder);
    fclose(logger);
}

/**
 * One thing to notice is, one branch can generate multiple Qs
 * we compute Rp := R \intersect p
 *
 * @return      N/A
**/
void InsMC::to_branch(int which_branch)
{
    /* Do nothing here */
    branch_id = which_branch;
    Q[Q_size++] = branch_v[branch_id];

    /* Then update upper and lower bound */
    compute_bound();
}

/**
 * Resume from branch
 *
 * @return      N/A
**/
void InsMC::from_branch()
{
    /* No need to do anything */
    branch_id = -1;
    Q_size--;
}

/**
 * Compute the upper and lower bound
 *
 * @return      N/A
**/
void InsMC::compute_bound(){
    /* Upper bound is easy */
    int which_vet = tmp_R[branch_id];
    int color = C_map[which_vet];
    upper_bound = Q_size + color;

    /* Lower bound is itself */
    lower_bound = Q_size;
}

/**
 * Calculate the combinations 
 * and get number of branches
 * this is actually the core algorithm
 * TODO This function is not neccessary.
 *
 * @return      N/A
**/
int InsMC::get_num_branches()
{
    /* number of possible branches is based on R 
       and the adjacent vertices of v */
    num_branch = 0;
    for(int i=0; i<R_size; ++i)
    {
        R_avail[tmp_R[i]] = true;
    }

    /* main body, check every vertex in R 
       and its adjacencies */
    for(int i=0; i<R_size; ++i)
    {
        int v     = tmp_R[i];
        int start = v == 0 ? 0 : v_idx[v-1]; 
        int end   = v_idx[v];

        for(int j=start; j<end; ++j)
        {
            int to = e_idx[j];
            if(R_avail[to] == true)
            {
                branch_v[num_branch++] = v;
                break;
            }
        }
    }

    /* resume */
    for(int i=0; i<R_size; i++)
    {
        R_avail[tmp_R[i]] = false;
    }

    return num_branch;
}

/**
 * From adj to encode
 * for the storage of the graph
 *
 * @param[out]      encode          encode to be written.
 *
 * @return      size of the encode.
**/
int 
InsMC::get_encode(int *encode)
{
    for(int i=0; i<Q_size; ++i)
    {
        encode[i] = Q[i];
    }
    return Q_size;
}

/**
 * Encode is stored in CSR format then transformed into 
 * adj format
 *
 * @param[in]       encode          input of the encode for the graph.
 * @param[in]       size            size of the encode.
 *
 * @return      N/A
**/
void 
InsMC::to_encode(int *encode, int size)
{
    for(int i=0; i<size; ++i) 
    {
        Q[i] = encode[i];
    }

    Q_size = size;

    get_full_R();

    for(int i=0; i<Q_size; ++i)
    {
        join_R(Q[i]);
    }

    get_deg_sorted();

    assign_C();
}

/**
 * Print the content of the encode
 *
 * @return      N/A
**/
void InsMC::print_encode()
{
    printf("encode: ");
    for(int i=0; i<Q_size; i++) printf("%d ", Q[i]);
    printf("\n");
}

/**
 * Copy the 'best' encode into current code
 *
 * @param[in]       encode          input of the encode for the graph.
 * @param[in]       size            size of the encode.
 *
 * @return      N/A
**/
void InsMC::copy_code(int *encode, int size)
{
    for(int i=0; i<size; i++) 
    {
        Q[i] = encode[i];
    }
    Q_size = size;
}

/**
 * Get the best value
 *
 * @return      N/A
**/
int InsMC::get_value()
{
    return upper_bound;
}

/**
 * get number of possible elements in the search list
 *
 * @return      N/A
**/
void
InsMC::get_num_elem(){
    num_count = upper_bound;
}

/**
 * private information is for debugging purposes
 *
 * @return      N/A
**/
int 
InsMC::return_private_info()
{
    if(Q_size<3)
    {
        return 0;
    }
    return 0;
}

/**
 * Allocate memory for basic graph data structures. 
 *
 * @param[in]       v_num       number of vertices
 * @param[in]       e_num       number of edges 
 *
 * @return      N/A
**/
void InsMC::allocate_data_structure(const int v_num, 
        const int e_num)
{
    v_size   = v_num;
    e_size   = e_num;
    v_idx    = new int[v_size + 1];
    e_idx    = new int[e_size + 1];
    Q        = new int[v_size + 1];
    R        = new int[v_size + 1];
    tmp_R    = new int[v_size + 1];
    R_deg    = new int[v_size + 1]; 
    R_avail  = new bool[v_size+1];
    C_map    = new int[v_size+1];
    branch_v = new int[v_size+1];
    for(int i=0; i<v_size; i++) 
    {
        v_idx[i] = 0;
    }
}

/**
 * scan the graph to get some basic information
 *
 * @return      N/A
**/
void
InsMC::init_graph()
{
    /* Declare variables */
    int v_num = OPTKIT_NULL;
    int e_num = OPTKIT_NULL;
    int g_num = OPTKIT_NULL;
    int v_id  = OPTKIT_NULL;
    int v_to  = OPTKIT_NULL;
    int sum   = OPTKIT_NULL;
    int color = OPTKIT_NULL;

    /* Read the head information, 
       and allocate according variables */
    FILE *reader;
    if((reader = fopen(in_file, "r")) == NULL)
    {
        printf("the file %s you input does not exist!\n", in_file);
        exit(1);
    }
    else if(fscanf(reader, "%d %d %d\n", &v_num, &g_num, &e_num)==EOF)
    {
        printf("error here\n");
    }

    int start_read_pos = ftell(reader);

    allocate_data_structure(v_num, e_num);

    /* Scan the real content */
    for(int i=0;i<e_num;i++)
    {
        if(fscanf(reader, "%d %d %d\n",&v_id , &v_to, &color)==EOF)
            printf("error here\n");
        v_idx[v_id] += 1;
    }

    sum=0;
    for(int i=0; i<v_num; i++)
    {
        sum += v_idx[i];
        v_idx[i]=sum;
    }

    /* Populate the content */
    fseek (reader, start_read_pos, SEEK_SET);

    /* this is for the start of the different postitions */
    int *idx = new int[v_size+1];
    for(int i=0; i<v_size; i++)
    {
        idx[i] = i == 0 ? 0 : v_idx[i-1];
    }

    /* Real read part */
    for(int i=0; i<e_num; i++)
    {
        if(fscanf(reader, "%d %d %d\n",&v_id , &v_to, &color)==EOF)
        {
            printf("error here\n");
        }

        int pos = idx[v_id];
        e_idx[pos] = v_to;
        ++idx[v_id];
    }

    /* e_idx also needs to be sorted, for the purpose of join operations */
    for(int i=0; i<v_size; i++)
    {
        int start = i==0?0:v_idx[i-1];
        int end = v_idx[i];
        Utils::q_sort(e_idx, start, end-1);
    }

    /* Write log */
    if(logger != NULL)
    {
        fprintf(logger, "====GRAPH====\n");
        for(int i=0; i<v_size ; i++)
        {
            int start = i==0 ? 0 : v_idx[i-1];
            int end   = v_idx[i];
            for(int j=start; j<end; j++)
            {
                fprintf(logger, "%d %d\n", i, e_idx[j]);
            }
        }
    }

    /* Frees */
    free(idx);

    fclose(reader);
}


/**
 * Detect connected components
 * and output them to file
 * use BFS algorithms
 *
 * @return      N/A
**/
void InsMC::separate_CC()
{
    bool *visited = new bool[v_size+1];
    for(int i=0; i<v_size; ++i) 
    {
        visited[i] = false;
    }
    int num_visited =0;
    int *queue = new int[v_size+1];
    int q_idx = 0;
    int count =0;
    while(true)
    {
        count++;

        //detect a un-visited vertex
        int seed = -1;
        for(int i=0; i<v_size; i++)
        {
            if(visited[i]==false)
            {
                int start = i==0?0:v_idx[i-1];
                int end = v_idx[i];
                int deg = end-start;
                if(deg>=1)
                {
                    seed = i;
                    visited[seed] = true;
                    break;
                }
                else
                {
                    visited[i] = true;
                }
            }
        }

        if(seed==-1)
        {
            break;
        }

        //perform bfs
        char f_buf[100];
        sprintf(f_buf,"%s%d.gr", CC_folder, count);
        FILE *writer = fopen(f_buf, "w");
        q_idx=0;
        queue[q_idx++] = seed;
        int edge_count =0;
        while(q_idx>0)
        {
            int v = queue[--q_idx];
            int start = v==0?0:v_idx[v-1];
            int end = v_idx[v];
            for(int i=start; i<end; i++)
            {
                if(visited[e_idx[i]]==false)
                {
                    visited[e_idx[i]] = true;
                    queue[q_idx++] = e_idx[i];  
                    fprintf(writer, "%d %d\n", v, e_idx[i]);
                    fprintf(writer, "%d %d\n", e_idx[i], v);
                    edge_count += 2;
                }
            }
        }
        fclose(writer);
    }

    printf("there are %d number of CCs in total\n", count);
}


/**
 * at the very begining R is the full list of all vertices
 * when a v is added into Q
 * Q = Q \union v
 * R = R \join \gamma(v)
 * \gamma(v) is the set of vertices that adjacent to v
 *
 * @return      N/A
**/
void
InsMC::join_R(int v)
{
    int m = 0;
    int n = v==0?0:v_idx[v-1];
    int v_end = v_idx[v];
    R_size =0;
    while(m<v_size && n<v_end)
    {
        if(R[m] == e_idx[n])
        {
            R[R_size++] = R[m];
            m++;
            n++;
        }
        else if(R[m] > e_idx[n])
        {
            n++;
        }
        else if(R[m] < e_idx[n])
        {
            m++;
        }
    }
    //keep R sorted
    Utils::q_sort(R, 0, R_size-1);
}

/**
 * this is to get a sorted list of each vertex in the graph
 *
 * @return      N/A
**/
void
InsMC::get_deg_sorted()
{
    /* initialize some info */
    for(int i=0; i<R_size; i++)
    { 
        tmp_R[i] = R[i];
    }
    /* get degrees first */
    for(int i=0; i<R_size; i++)
    {
        int v     = R[i];
        int start = v== 0 ? 0 : v_idx[v-1];
        int end   = v_idx[v];
        R_deg[i]  = (end - start)*-1;
    }
    /* then sort */
    Utils::q_sort_two(R_deg, tmp_R, 0, R_size-1);
}


/**
 * Assign C to each vertex used for the upper bound update
 * using heuristic mentioned in Janes Konc's paper:
 * An improved branch and bound algorithm for the maximum clique problem
 * some improvement might be possible using vertex cover
 *
 * @return      N/A
**/
void InsMC::assign_C()
{
    int total_added=0;
    int now_C = 1;

    for(int i=0; i<R_size; ++i)
    {
        R_avail[tmp_R[i]] = true;
        C_map[tmp_R[i]] = OPTKIT_NULL;
    }
    /* List each tree */
    while(total_added < R_size)
    {
        /* You have to put a first vertex in the C 
           find the first, it's neccessary to have this. */
        int first_pos = find_first_avail_pos();
        int first_vet = tmp_R[first_pos];
        C_map[first_vet] = now_C;
        ++total_added;

        /* For each round, scan the vertex list; 
           and add the vertex that is adjacent to
           current list. */
        for(int i=0; i<R_size; ++i)
        {
            if(i == first_pos) 
                continue;
            /* each vertex will try to see 
               if it's adjacet is already in. */
            int v     = tmp_R[i];
            int start = v == 0 ? 0 : v_idx[v-1];
            int end   = v_idx[v];
            bool is_independent = true;
            for(int j=start; j<end; ++j)
            {
                if(true == R_avail[e_idx[j]] 
                    && C_map[e_idx[j]] == now_C)
                {
                    is_independent = false;
                    break;
                }
            }
            if(true == is_independent 
                    && C_map[v] == OPTKIT_NULL)
            {
                C_map[v] = now_C;
                ++total_added;
            }
        }
        ++now_C;
    }
    for(int i=0; i<R_size; ++i) 
    {
        R_avail[tmp_R[i]] = false;
    }

    max_C = now_C-1;

    if(logger != NULL)
    {
        /* For parallel computing purpose, 
           v need to be sorted. */
        fprintf(logger, "====COLOR====\n");
        int min_v = OPTKIT_NULL;
        for(int i=0; i<R_size; ++i)
        {
            int tmp_min_v = R_size + 1;
            for(int j=0; j<R_size; j++)
            {
                int v = tmp_R[j];
                if(v < tmp_min_v && v > min_v)
                {
                    tmp_min_v   = v;
                }
            }
            min_v = tmp_min_v;
            fprintf(logger, "vertex: %d color: %d\n", min_v, C_map[min_v]);
        }
    }
}

/**
 *
 * @return      N/A
**/
void InsMC::get_num_triangles()
{
    int num_triangle = 0;
    int *pos_larger = (int*)malloc(sizeof(int)*v_size);
    for(int i=0; i<v_size; i++)
    {
        int start = i==0?0:v_idx[i-1];
        int end = v_idx[i];
        for(int j=start; j<end; j++)
        {
            if(e_idx[j]>i)
            {
                pos_larger[i] = j;
                break;
            }
        }
    }
    for(int i=0; i<v_size; i++)
    {
        int start = i==0?0:v_idx[i-1];
        int end = v_idx[i];
        for(int j=start; j<end; j++)
        {
            int to = e_idx[j];
            int to_start = to==0?0:v_idx[to-1];
            int to_end = v_idx[to];
            int m = start+j; 
            int n = pos_larger[to];
            while(m<end && n<to_end)
            {
                if(e_idx[m]==e_idx[n])
                {
                    m++;
                    n++;
                    num_triangle++;
                }
                else if(e_idx[m]<e_idx[n])
                {
                    m++;
                }
                else if(e_idx[m]>e_idx[n])
                {
                    n++;
                }
            }
        }
    }
    printf("total number of triangles %d\n", num_triangle);
}

/**
 * Visualize csr graph
 *
 * @return      N/A
**/
void InsMC::vis_csr()
{
    FILE *writer = fopen("vis/vis_csr.dot", "w");
    fprintf(writer, "graph{\n");
    for(int i=0; i<v_size; i++)
    {
        int start = i==0?0:v_idx[i-1];
        int end = v_idx[i];
        for(int j=start; j<end; j++)
        {
            fprintf(writer, "%d -- %d [color = red]\n", i, e_idx[j]);
        }
    }
    fprintf(writer, "}\n");
}

/**
 * print degree information
 *
 * @return      N/A
**/
void InsMC::print_deg()
{
    for(int i=0; i<R_size; i++)
    {
        printf("%d:%d ", tmp_R[i], R_deg[i]);
    }
    printf("\n");
}

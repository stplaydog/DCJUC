#include "insmed.h"
#include "listint.h"
#include <fstream>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstdio>

extern int g_count;

///////////////////////////////////////////////////////////////////
//
// PUBLIC FUNCTIONS
//
//////////////////////////////////////////////////////////////////

/**
 * this is the constructor to initialize the three input orders
 * and the median order, median order is initialized randomly
 * it is supposed to implement the adequate subgraph algorithm
 * but due to the limit of time, this task is postponded to the future work.
 *
 * @param[in]       file            where to read the graph file
 * @param[in]       tmp_folder      folder to keep tmp files
 * @param[in]       tid             thread number
 * @param[in]       dm              distance mode exemplar:1 matching:2
 * @param[in]       uh              if use heuristic
 * @param[in]       th              threshold for heuristic
 *
 * @return      N/A
**/
InsMed::InsMed(const char *file, const char* tmp_folder, int tid, const char *lf, const int *dm, const bool *uh, const int *th) :
    Instance(lf)
{
    m_file       = NULL;
    m_tmp_folder = NULL;
    num_t = tid;
    /* read order from file */
    if(file != NULL)
    {
        m_file = new char[OPTKIT_FILE_SIZE];
        snprintf(m_file, OPTKIT_FILE_SIZE, "%s", file);
    }
    if(tmp_folder != NULL)
    {
        m_tmp_folder = new char[OPTKIT_FILE_SIZE];
        snprintf(m_tmp_folder, OPTKIT_FILE_SIZE, "%s", tmp_folder);
    }
    init_graph();

    /* decompose the graph by ASs */
    rename_graph_by_car();

    /* change CSR as well*/
    from_adj_to_csr();

    /* now it's time to init median genome */
    init_median();
    marriage_cap_vet();
    detect_n_remove_conflict_edge();
    num_count = 10;

    /* init distance based tools */
    init_dis_tools();

    /* compute the initial bound */
    is_branching = false;
    compute_bound();

    if(dm != NULL)
    {
        dis_mode = *dm;
    }
    if(uh != NULL)
    {
        use_heu = *uh;
    }
    if(th != NULL)
    {
        thresh = *th;
    }
}


/**
 * copy constructor
 *
 * @param[in]       other       the instance tobe copied
 *
**/
InsMed::InsMed(const InsMed &other) :
    Instance(NULL)
{
}

/**
 * branching is using reverse operation
 *
 * @param[in]       which_branch        to branch based on the branch code
 **/
    void
InsMed::to_branch(int which_branch)
{
    int encode[200];
    /* shrink first */
    vector<int> ft_now = ft[which_branch];
    for(int i=0; i<(int)ft_now.size(); i++)
    {
        encode[i] = ft_now[i];
    }
    shrink(encode);
    from_csr_to_adj(3);
    detect_n_remove_conflict_edge();
    is_branching = true;
    branch_selected = which_branch;
    compute_bound();
    from_adj_to_csr();
}

/**
 * use reverse operation as well
 **/
    void
InsMed::from_branch()
{
    int encode[200];
    /* shrink first */
    vector<int> ft_now = ft[branch_selected];
    for(int i=0; i<(int)ft_now.size(); i++)
    {
        encode[i] = ft_now[i];
    }
    resume(encode);
    branch_selected = OPTKIT_NULL;
}

/**
 * for lin kernighan algorithm there is no need to compute bound, just get the score
 **/
    void
InsMed::compute_bound()
{
    g_count = 0;
    int encode[30000];
    int num_comp = 0;
    for(int i=0; i<3; i++)
    {
        /* get the encode for insdis */
        int c1, c2;
        if(0==i)
        {
            c1 = 0;
            c2 = 3;
        }
        else if(1 == i)
        {
            c1 = 1;
            c2 = 3;
        }
        else if(2 == i)
        {
            c1 = 2;
            c2 = 3;
        }

        //printf("to %d th branch color %d!\n", branch_selected, i);
        /* detect number of component first, 
         * this function will also call bfs*/
        int c_type[v_size];
        int num_cc = detect_comps(c_type, c1, c2);
        //print_c_type(c_type);
        for(int k=1; k<=num_cc; k++)
        {
            //printf("computing %d/%d th comp is_branching %s\n", 
            //    k, num_cc, is_branching?"true":"false");

            /* copy only one component to encode  */
            int encode_size = get_encode_by_component(c_type, encode, k, c1, c2);
            rename_encode_by_comp(encode);
            //if(k==1)
            //{
            //    printf("%d %d %d ", encode[0], encode[1], encode[2]);
            //    for(int c=0; c<2; c++)
            //    {
            //        for(int m=0; m<encode[0]; m++)
            //        {
            //            int start = m==0 ? 0 : encode[3+encode[0]*c + m-1];
            //            int end = encode[3+encode[0]*c + m];
            //            start += 3+2*encode[0] + c*encode[1];
            //            end += 3+2*encode[0] + c*encode[1];
            //            for(int n=start; n<end; n++)
            //            {
            //                printf("%d->%d colort %d | ", m, encode[n], c);
            //            }
            //        }
            //    }
            //    printf("\n");
            //}
            if(true == check_if_single_cc(encode))
            {
                num_comp += 1;
            }
            else if(encode[0]==4)
            {
                num_comp += 2;
            }
            else
            {
                /* to encode */
                vis_encode(encode);
                ins[i][0]->to_encode(encode, encode_size);
                /* using bnb to compute distance */
                l_d[i]->reset_list();
                ins[i][0]->compute_bound();
                l_d[i]->upper_bound = ins[i][0]->upper_bound;
                l_d[i]->lower_bound = ins[i][0]->lower_bound;
                l_d[i]->base = ins[i][0]->lower_bound;
                l_d[i]->bnb(ins[i], 0);
                num_comp += (l_d[i]->lower_bound+l_d[i]->upper_bound)/2;
            }
        }
    }
    /* get score */
    score = 3*v_size/2 - num_comp - cycle;
    upper_bound = 3*v_size/2 - num_comp - cycle;
}

/**
 * get number of possible branches
 *
 * @return      number of branches
 **/
    int 
InsMed::get_num_branches()
{
    int num_branches = 0;
    ft.clear();
    if(use_heu == true)
    {
        indentify_comps();
    }
    for(int i=0; i<v_size; i++)
    {
        int start_one = i==0?0:v_idx[3][i-1];
        int end_one = v_idx[3][i];
        for(int j=start_one; j<end_one; j++)
        {
            int from_one = i;
            int to_one = e_idx[3][j];
            for(int m=0; m<v_size; m++)
            {
                int start_two = m==0 ? 0 : v_idx[3][m-1];
                int end_two = v_idx[3][m];
                for(int n=start_two; n<end_two; n++)
                {
                    int from_two = m;
                    int to_two = e_idx[3][n];    
                    if(from_one != to_one 
                            && from_one != from_two 
                            && from_one != to_two 
                            && to_one != from_two 
                            && to_one != to_two 
                            && from_two != to_two
                            && from_one != CAP 
                            && to_one != CAP 
                            && from_two != CAP 
                            && to_two != CAP)
                    { 
                        /* this is the rule to chose whom to exchange */
                        if(use_heu == true)
                        {
                            int score = 0;
                            for(int color=0; color<3; color++)
                            {
                                if(comp_id[color][i] == comp_id[color][m] && comp_id[color][i] != -1)
                                {
                                    if(comp_type[color][comp_id[color][i]]==TYPE_SING)
                                    {
                                        score+=2;
                                    }
                                    else if(comp_type[color][comp_id[color][i]]==TYPE_MUL)
                                    {
                                        score+=1;
                                    }
                                }
                            }
                            if(score > thresh)
                            {
                                num_branches += 2;
                                add_possible_branch(from_one, to_one, from_two, to_two);
                            }
                        }
                        else
                        {
                            num_branches += 2;
                            add_possible_branch(from_one, to_one, from_two, to_two);
                        }
                    }
                }
            }
        }
    }
    return num_branches;
}

/**
 * get encode of current state
 *
 * @param[out]      encode          the array to keep encode
 *
 * @return      size of the encode
 **/
    int 
InsMed::get_encode(int* encode)
{
    encode[0] = v_size;
    encode[1] = e_size[3];
    for(int i=0; i<v_size; i++)
    {
        encode[2+i] = v_idx[3][i];
        //printf("%d ", encode[2+i]);
    }
    for(int i=0; i<e_size[3]; i++)
    {
        encode[2+v_size+i] = e_idx[3][i];
        //printf("%d ", encode[2+v_size+i]);
    }
    //printf("\n");
    return 2+v_size+e_size[3];
}

/**
 * given an encode, transfer current state according to encode
 *
 * @param[in]       encode          input of the encode
 * @param[in]       size            size of the encode
 **/
    void
InsMed::to_encode(int* encode, int size)
{
    copy_code(encode, size);
    compute_bound();
}

/**
 * print the content of the encode
 **/
void
InsMed::print_encode()
{
}

/**
 * copy the encode to the current object
 **/
void
InsMed::copy_code(int *encode, int size)
{
    e_size[3] = encode[1];
    for(int i=0; i<v_size; i++)
    {
        v_idx[3][i] = encode[2+i];
        //printf("%d ", v_idx[3][i]);
    }
    for(int i=0; i<e_size[3]; i++)
    {
        e_idx[3][i] = encode[2+v_size+i];
        //printf("%d ", e_idx[3][i]);
    }
    //printf("\n");
}

/**
 * get the value of the current object
 **/
int 
InsMed::get_value()
{
    return score; 
}

/**
 * average size of the current object
 **/
void
InsMed::compute_num_elem()
{
    num_count = 10;
}


///////////////////////////////////////////////////////////////////
//
// PRIVATE FUNCTIONS
//
///////////////////////////////////////////////////////////////////

/**
 * let's output our CSR into a file
 *
 * @param[in]       file        filename
 * @param[in]       c1          color 1
 * @param[in]       c2          color 2
 *
**/
void InsMed::csr_to_file(char *file, int c1, int c2)
{
    FILE *writer = fopen(file, "w");
    if(!writer)
    {
        printf("File %s does not exits!\n", file);
        exit(1);
    }

    /* write graph to file */
    fprintf(writer, "%d %d %d %d\n", v_size, 2, e_size[c1]+e_size[c2], 0);
    for(int i=0; i<2; i++)
    {
        int c = i==0 ? c1 : c2;
        for(int j=0; j<v_size; j++)
        {
            int start = j==0 ? 0 : v_idx[c][j-1];
            int end = v_idx[c][j];
            for(int k=start; k<end; k++)
            {
                int from = j;
                int to = e_idx[c][k];
                fprintf(writer, "%d %d %d\n", from, to, i+1);
            }
        }
    }
    fclose(writer);
}

/**
 * init some inner data structure for distance computation
 * one list instance and three insdis instances
**/
void InsMed::init_dis_tools()
{
    l_d = (List**)malloc(sizeof(List*)*3);
    ins = (Instance***)malloc(sizeof(Instance**)*3);
    /* now I only need to instantiate once */
    for(int i=0; i<3; i++)
    {
        int c1, c2;
        if(0==i)
        {
            c1 = 0;
            c2 = 1;
        }
        else if(1 == i)
        {
            c1 = 0;
            c2 = 2;
        }
        else if(2 == i)
        {
            c1 = 1;
            c2 = 2;
        }
        /* write CSR to file first */
        char csr_file[100];
        sprintf(csr_file, "%s%d", m_tmp_folder, i);
        csr_to_file(csr_file, c1, c2);
        ins[i] = (Instance**)new InsDis*[2];
        if(dis_mode == DIS_MODE_EXEM)
        {
            for(int j=0; j<2; j++)
            {
                ins[i][j] = new InsDis(csr_file, true, NULL);
            }
        }
        else if(dis_mode == DIS_MODE_MATC)
        {
            for(int j=0; j<2; j++)
            {
                ins[i][j] = new InsDis(csr_file, false, NULL);
            }
        }
        if(ins[i][0]->upper_bound == ins[i][0]->lower_bound)
        {
            ins[i][0]->lower_bound -= 1;
        }
        int buck_size = ins[i][0]->upper_bound + 2;
        int list_size = 1000000;
        int base = ins[i][0]->lower_bound-1;
        int num_tt = 1;
        int num_elem = ins[i][0]->num_count;
        //printf("buck_size %d list_size %d base %d num_tt %d, num_elem %d\n", 
        //    buck_size, list_size, base, num_tt, num_elem);
        l_d[i] = new IntList(buck_size, list_size, base, num_tt, true, num_elem);
    }
}







/**
 * Based on the input of gene order, construct a 
 * CSR graph.
 **/
void InsMed::init_graph()
{
    /* scan the graph to get the basic information */
    scan_graph();
    read_graph();
    for(int c=0; c<3; c++)
    {
        from_csr_to_adj(c);
    }
}

/**
 * scan the graph to get some basic information
 **/
void InsMed::scan_graph(){
    /* declare variables */
    FILE *stream;
    char file_buf[500];
    int v_num;
    int e_num;
    int g_num;
    int v_id=0;
    int v_to=0;
    int color=0;
    int sum=0;

    /* read the head information, and allocate according variables */
    sprintf(file_buf, "%s", m_file);
    if((stream=fopen(file_buf,"r"))==NULL)
    {
        printf("the file %s you input does not exist!\n", file_buf);
        exit(1);
    }
    if(fscanf(stream, "%d %d %d\n", &v_num, &g_num, &e_num)==EOF)
    {
        printf("error here\n");
    }
    v_size = v_num;
    v_idx = new int*[4];
    e_idx = new int*[4];
    comp_id = new int*[4];
    comp_type = new int*[4];
    for(int i=0; i<4; i++)
    {
        v_idx[i] = new int[v_size*4];
        comp_id[i] = new int[v_size*4];
        comp_type[i] = new int[v_size*4];
        for(int j=0; j<v_size*4; j++)
        {
            v_idx[i][j]  = 0;
        }
        e_idx[i] = new int[e_num*4];
        for(int j=0; j<e_num*4; j++)
        {
            e_idx[i][j]  = 0;
        }
        e_size[i] = 0;
    }
    /* scan the real content */
    for(int i=0;i<e_num;i++)
    {
        if(fscanf(stream, "%d %d %d\n",&v_id , &v_to, &color)==EOF)
        {
            printf("error here\n");
        }
        e_size[color-1] += 1;
        v_idx[color-1][v_id] +=1;
    }
    for(int i=0;i<3;i++)
    {
        sum=0;
        for(int j=0; j<v_num; j++)
        {
            sum += v_idx[i][j];
            v_idx[i][j]=sum;
        }
    }
    /* allocate adjs */
    adj = new int**[4];
    a_idx = new int*[4];
    for (int c=0; c<4; c++)
    {
        adj[c] = new int*[v_size*4];
        a_idx[c] = new int[v_size*4];
        for (int i=0; i<v_size*4; i++)
        {
            /* can't have duplication of more than 100 */
            adj[c][i] = new int[100];
        }
    }
    fclose(stream);
}

/**
 * read the graph into local data structure
 **/
void InsMed::read_graph()
{
    FILE *stream;
    char file_buf[500];
    int v_num;
    int e_num;
    int g_num;
    int v_id=0;
    int v_to=0;
    int color=0;
    /* read basic information */
    sprintf(file_buf, "%s", m_file);
    if((stream=fopen(file_buf,"r"))==NULL)
    {
        printf("the file %s you input does not exist!\n", file_buf);
        exit(1);
    }
    if(fscanf(stream, "%d %d %d\n", &v_num, &g_num, &e_num)==EOF)
    {
        printf("error here\n");
    }
    /* this is for the start of the different postitions */
    int **idx = new int*[3];
    for(int i=0; i<3; i++)
    {
        idx[i] = new int[v_size];
        for(int j=0; j<v_size; j++)
        {
            idx[i][j] = j==0?0:v_idx[i][j-1];
        }
    }
    /* real read part */
    for(int i=0; i<e_num; i++)
    {
        if(fscanf(stream, "%d %d %d\n",&v_id , &v_to, &color)==EOF)
        {
            printf("error here\n");
        }
        int pos = idx[color-1][v_id];
        e_idx[color-1][pos] = v_to;
        idx[color-1][v_id]++;
    }
    /* frees */
    for(int i=0; i<3; i++)
    {
        delete [] idx[i];
    }
    delete [] idx;
    fclose(stream);
}

/**
 * transfer csr to adj for the convinience of comuting
 *
 * @param[i]        c               which order to transfer
 **/
void InsMed::from_csr_to_adj(int c)
{
    for (int i=0; i<v_size; i++)
    {
        a_idx[c][i] = 0;
        /* start copying */
        int start = i==0?0:v_idx[c][i-1];
        int end = v_idx[c][i];
        for (int j=start; j<end; j++)
        {
            adj[c][i][a_idx[c][i]++] = e_idx[c][j];
        }
    }
}

/**
 * detect ASs and remove them, add cycle count accordingly
 *
 * TODO refactor code to make it concise
 **/
    void 
InsMed::rename_graph_by_car()
{
    cycle = 0;
    /* detect ASs that three edges of one vertex connected to the same vertex */
    for(int i=0; i<v_size; i++)
    {
        if(a_idx[0][i]==0 || a_idx[1][i]==0 || a_idx[2][i]==0)
        {
            /* no edges */
            continue;
        }
        else if(a_idx[0][i]>1 || a_idx[1][i]>1 || a_idx[2][i]>1)
        {
            /* degree doesn't satisfy condition */
            continue;
        }
        else if(adj[0][i][0] == adj[1][i][0] && adj[1][i][0] == adj[2][i][0])
        {
            /* now we find this type of ASs */
            int to = adj[0][i][0];
            if(to == CAP)
            {
                continue;
            }
            else if(a_idx[0][to]>1 || a_idx[1][to]>1 || a_idx[2][to]>1)
            {
                continue;
            }
            else
            {
                a_idx[0][i] = a_idx[1][i] = a_idx[2][i] =0;
                a_idx[0][to] = a_idx[1][to] = a_idx[2][to] =0;
                cycle += 3;
            }
        }
    }
    /* detect two edges of a vertex connected to the same vertex */
    for(int color=0; color<3; color++)
    {
        /* assign colors */
        int c1=-1, c2=-1, c3=-1;
        if(color==0)
        {
            c1=0; 
            c2=1; 
            c3=2;
        }
        else if(color==1)
        {
            c1=0; 
            c2=2; 
            c3=1;
        }
        else if(color==2)
        {
            c1=1; 
            c2=2; 
            c3=0;
        }

        /* now lets detect */
        for(int i=0; i<v_size; i++)
        {
            /* conditions that doesn't satisfy */
            if(a_idx[c1][i]>1 || a_idx[c2][i]>1 || a_idx[c3][i]>1)
            {
                continue;
            }
            else if(a_idx[c1][i]==0 || a_idx[c2][i]==0)
            {
                continue;
            }
            int to1 = adj[c1][i][0]; 
            int to2 = adj[c2][i][0]; 
            /* can not connected to CAP */
            if(to1 == CAP || to2 == CAP)
            {
                continue;
            }
            if(to1 == to2)
            {
                if(a_idx[c1][to1]>1 || a_idx[c2][to1]>1 || a_idx[c3][to1]>1)
                {
                    continue;
                }
                if((a_idx[c3][i]==0 && a_idx[c3][to1]!=0) ||
                        (a_idx[c3][i]!=0 && a_idx[c3][to1]==0))
                {
                    continue;
                }
                if(a_idx[c3][i]==0 && a_idx[c3][to1]==0)
                {
                    a_idx[c1][i] = 0; 
                    a_idx[c2][i] = 0; 
                    a_idx[c3][i] = 0;
                    a_idx[c1][to1] = 0; 
                    a_idx[c2][to1] = 0; 
                    a_idx[c3][to1] = 0;
                    cycle+=2;
                    continue;
                }
                int neighbor = adj[c3][i][0];
                int to_neighbor = adj[c3][to1][0];
                if(neighbor==CAP || to_neighbor==CAP)
                {
                    continue;
                }
                for(int j=0; j<a_idx[c3][neighbor]; j++)
                {
                    if(adj[c3][neighbor][j] == i)
                    {
                        adj[c3][neighbor][j] = to_neighbor;
                    }
                }
                for(int j=0; j<a_idx[c3][to_neighbor]; j++)
                {
                    if(adj[c3][to_neighbor][j] == to1)
                    {
                        adj[c3][to_neighbor][j] = neighbor;
                    }
                }
                a_idx[c1][i] = 0; 
                a_idx[c2][i] = 0; 
                a_idx[c3][i] = 0;
                a_idx[c1][to1] = 0; 
                a_idx[c2][to1] = 0; 
                a_idx[c3][to1] = 0;
                cycle+=2;
            }
        }
    }
}

/**
 * Turn adj back into CSR
 *
**/
void InsMed::from_adj_to_csr()
{
    /* v_idx */
    for(int i=0; i<3; i++)
    {
        int sum=0;
        for(int j=0; j<v_size; j++)
        {
            sum += a_idx[i][j];
            v_idx[i][j] = sum;
        }
    }
    /* e_idx */
    for(int i=0; i<3; i++)
    {
        int pos =0;
        for(int j=0; j<v_size; j++)
        {
            for(int k=0; k<a_idx[i][j]; k++)
            {
                e_idx[i][pos++] = adj[i][j][k];
            }
        }
    }
    /* change e_size */
    for(int c=0; c<3; c++)
    {
        e_size[c] = v_idx[c][v_size-1];
    }
}

/**
 * this function is used to initialize the median adj in the graph
 * 1) choose the content of the median adj
 * 2) randomly assign the adj
 **/
void 
InsMed::init_median(){
    /* Assigning initial degrees for the CSR */
    for(int i=0; i<v_size; i++)
    {
        v_idx[3][i] = OPTKIT_NULL;
    }
    for(int i=0;i<v_size;i++)
    {
        if(v_idx[3][i] == OPTKIT_NULL)
        {
            int deg = get_deg(i);
            v_idx[3][i] = deg;
        }
    }
    for(int i=0;i<v_size;i++)
    {
        if(v_idx[3][i] == OPTKIT_NULL)
        {
            v_idx[3][i] = 0;
        }
    }

    /* based on degree, compute v_idx, and e_size*/
    int sum = 0;
    for(int i=0; i<v_size; i++)
    {
        sum += v_idx[3][i];
        v_idx[3][i] = sum;
    }
    e_size[3] = sum;

    /* randomly init edges */
    for(int i=0;i<e_size[3]; i++)
    {
        e_idx[3][i] = OPTKIT_NULL;
    }
    int *rev_v = new int[e_size[3]];
    for(int i=0; i<v_size; i++)
    {
        int start = i==0?0:v_idx[3][i-1];
        int end = v_idx[3][i];
        for(int j=start; j<end; j++)
        {
            rev_v[j] = i;
        }
    }
    for(int i=0; i<e_size[3]; i++)
    {
        if(e_idx[3][i] == OPTKIT_NULL)
        {
            int to = find_next_edge(rev_v, i);
            if(to != CAP)
            {
                e_idx[3][i] = rev_v[to];
                e_idx[3][to] = rev_v[i];
            }
            else
            {
                e_idx[3][i] = CAP;
            }
        }
    }
    /* turn from CSR format to adj format */
    from_csr_to_adj(3);
    /* delete temp files */
    delete rev_v;
}


/**
 * find the next available edge
 *
 * @param[in]       rev_v       reverse vertex intex
 * @param[in]       from        vertex 'from'
 *
 * @return      next available edge
 **/
int InsMed::find_next_edge(int *rev_v, int from)
{
    int ret = -1;
    int conflict = 0;
    if(rev_v[from] % 2 == 0)
    {
        conflict = rev_v[from]+1;
    }
    else
    {
        conflict = rev_v[from]-1;
    }
    /* check CAR first */
    int v_from = rev_v[from];
    int start1 = v_from ==0 ? 0 : v_idx[0][v_from-1];
    int end1 = v_idx[0][v_from];
    int start2 = v_from ==0 ? 0 : v_idx[1][v_from-1];
    int end2 = v_idx[1][v_from];
    int start3 = v_from ==0 ? 0 : v_idx[2][v_from-1];
    int end3 = v_idx[2][v_from];
    int vet = -1;
    for(int i=start1; i<end1; i++)
    {
        for(int j=start2; j<end2; j++)
        {
            if(e_idx[0][i] == e_idx[1][j] 
                && !edge_exist(v_from, e_idx[0][i]) 
                && e_idx[0][i] != CAP)
            {
                vet = e_idx[0][i];
                break;
            }
        }
        if(vet != -1)
        {
            break;
        }
    }
    if(vet == -1)
    {
        for(int i=start1; i<end1; i++)
        {
            for(int j=start3; j<end3; j++)
            {
                if(e_idx[0][i] == e_idx[2][j] 
                    && !edge_exist(v_from, e_idx[0][i]) 
                    && e_idx[0][i] != CAP)
                {
                    vet = e_idx[0][i];
                    break;
                }
            }
            if(vet != -1)
            {
                break;
            }
        }
    }
    if(vet == -1)
    {
        for(int i=start2; i<end2; i++)
        {
            for(int j=start3; j<end3; j++)
            {
                if(e_idx[1][i] == e_idx[2][j] 
                    && !edge_exist(v_from, e_idx[0][i])
                    && e_idx[1][i] != CAP)
                {
                    vet = e_idx[1][i];
                    break;
                }
            }
            if(vet != -1)
            {
                break;
            }
        }
    }
    if(vet != -1)
    {
        int start = vet==0 ? 0 : v_idx[3][vet-1];
        int end = v_idx[3][vet];
        for(int i=start; i<end; i++)
        {
            if(e_idx[3][i] == -1)
            {
                ret = i;
                break;
            }
        }
    }
    else if(ret == -1)
    {
        /* then assign randomly */
        for(int i=0; i<e_size[3]; i++)
        {
            /* check if there is already an edge connecting between
             * rev_v[from] and rev_v[i]*/
            int v_to = rev_v[i];
            int v_start = v_from==0 ? 0 : v_idx[3][v_from-1];
            int v_end = v_idx[3][v_from];
            bool is_redundant_vet = false;
            for(int j=v_start; j<v_end; j++)
            {
                if(v_to == e_idx[3][j])
                {
                    is_redundant_vet = true;
                    break;
                }
            }
            if(e_idx[3][i] == -1 
                    && rev_v[i] != conflict 
                    && v_to != v_from
                    && false == is_redundant_vet)
            {
                ret = i;
                break;
            }
        }
    }
    if(ret == -1)
    {
        ret = CAP;
    }
    //printf("-------------return value %d-->%d\n", rev_v[from], rev_v[ret]);
    return ret;
}

/**
 * check if edge already exist for a given vertex
 *
 * @param[in]       from            vertex from 
 * @param[in]       to              vertex to 
 *
 * @return      return true(exist) or false(not exist)
 *
**/
bool InsMed::edge_exist(int from, int to)
{
    bool ret = false;
    int start = from ==0 ? 0 : v_idx[3][from-1];
    int end = v_idx[3][from];
    for(int i = start; i<end; i++)
    {
        if(e_idx[3][i] == to)
        {
            ret = true;
        }
    }
    return ret;
}


/**
 * get degree of a given vertex
 *
 * @param[in]       pos         which pos of the vertex
 *
 * @return degree
 **/
int InsMed::get_deg(int vet)
{
    int start1 = vet==0 ? 0 : v_idx[0][vet-1];
    int end1 = v_idx[0][vet];
    int start2 = vet==0 ? 0 : v_idx[1][vet-1];
    int end2 = v_idx[1][vet];
    int start3 = vet==0 ? 0 : v_idx[2][vet-1];
    int end3 = v_idx[2][vet];
    int d1 = end1 - start1;    
    int d2 = end2 - start2;    
    int d3 = end3 - start3;    
    int ret = -1;
    if(d1==d2 && d1==d3)
    {
        ret = d1;
    }
    else
    {
        if(max_deg(vet)==1)
        {
            ret = 1;
        }
        if(max_deg(vet)>1)
        {
            ret = max_deg(vet)/2;
        }
        if(max_deg(vet)==0)
        {
            ret = 0;
        }
    }
    return ret;
}

/**
 * get the max dgree of a vertex
 *
 * @param[in]       vet         which vertex to calculate
 **/
int InsMed::max_deg(int vet)
{
    int start1 = vet==0 ? 0 : v_idx[0][vet-1];
    int end1 = v_idx[0][vet];
    int start2 = vet==0 ? 0 : v_idx[1][vet-1];
    int end2 = v_idx[1][vet];
    int start3 = vet==0 ? 0 : v_idx[2][vet-1];
    int end3 = v_idx[2][vet];
    int one = end1 - start1;    
    int two = end2 - start2;    
    int three = end3 - start3;    
    int ret = -1;
    if(one>two && one>three)
    {
        ret = one;
    }
    else if(two > three)
    {
        ret = two;
    }
    else
    {
        ret = three;
    }
    return ret;
}

/**
 * reconnect those vertices that is connected to CAP vertex
 * TODO because insdis doesn't have a complete mechanism to
 * process the CAP vertices
 *
 * @note        one vertex cannot connect to the conflict vertex
 *              (which is originated from the same genome).
 *              one vertex cannot connect to itself.
 *
 * **/
void InsMed::marriage_cap_vet()
{
    for(int i=0; i<v_size; i++)
    {
        int start = i==0 ? 0 : v_idx[3][i-1];
        int end = v_idx[3][i];
        for(int j=start; j<end; j++)
        {
            if(e_idx[3][j] == CAP)
            {
                /* search for another edge now */
                int another_vet;
                int another_edge;
                search_cap_edge(i, another_edge, another_vet);
                /* now perform reconnect */
                e_idx[3][another_edge] = i;
                e_idx[3][j] = another_vet;
            }
        }
    }
}

/**
 * find another vertex that can be connected to this vertex which
 * is connected to a CAP vertex.
 * @param[in]       one_vet         the given vertex tobe connected
 * @param[out]      ret_edge        which edge tobe reconnected
 * @param[out]      ret_vet         which vertex tobe reconnected
 *
 **/
void InsMed::search_cap_edge(int one_vet, int &ret_edge, int &ret_vet)
{
    ret_vet = -1;
    ret_edge = -1;

    /* get confliction vertex */
    int conflict = 0;
    if(one_vet % 2 == 0)
    {
        conflict = one_vet+1;
    }
    else
    {
        conflict = one_vet-1;
    }

    /* find the vet and edge */
    for(int i=0; i<v_size; i++)
    {
        if( i== conflict || i==one_vet)
        {
            continue;
        }
        int start = i==0 ? 0 : v_idx[3][i-1];
        int end = v_idx[3][i];
        for(int j=start; j<end; j++)
        {
            if(e_idx[3][j] == CAP)
            {
                /* search for another edge now */
                ret_edge = j;
                ret_vet = i;
                break;
            }
        }
        if(ret_edge != -1)
        {
            break;
        }
    }
}


/**
 * Identify different types of connected component 
 **/
    void 
InsMed::indentify_comps()
{
    bool visited[v_size];
    int queue[v_size];
    int q_idx =0;


    for(int color=0; color<3; color++)
    {
        for(int i=0; i<v_size; i++)
        {
            visited[i] = false;
            comp_id[color][i]  = OPTKIT_NULL;
        }
        int comp_num = 0;
        for(int i=0; i<v_size; i++)
        {
            int start = i== 0 ? 0 : v_idx[color][i-1];
            int end = v_idx[color][i];
            int deg = end - start;
            if( deg == 0 || visited[i] == true)
            {
                visited[i] = true;
                continue;
            }
            q_idx = 0;
            comp_id[color][i] = comp_num;
            comp_type[color][comp_num]= TYPE_SING;
            queue[q_idx++] = i;
            visited[i] = true;
            /* bfs to find connected components */
            while(q_idx > 0)
            {
                /* pop bfs list */
                int v = queue[--q_idx];
                for(int c=0; c<2; c++)
                {
                    int cc = c==0 ? color :3 ;
                    int start = v==0 ? 0 : v_idx[cc][v-1];
                    int end = v_idx[cc][v];
                    int deg = end - start;
                    if(deg>1)
                    {
                        comp_type[cc][comp_num]= TYPE_MUL;
                    }
                    for(int j=start; j<end; j++)
                    {
                        int to = e_idx[cc][j];
                        if(to != OPTKIT_NULL && to != CAP && visited[to]==false)
                        {
                            comp_id[color][to] = comp_num;
                            queue[q_idx++] = to;
                            visited[to] = true;
                        }
                    }
                }
            }
            comp_num++;
        }
    }
}

/**
 * shrink the graph
 *
 * @param[in]       encode      the encode to be shrunk
 *
 **/
void InsMed::shrink(int *encode)
{
    int from_one = encode[0];
    int to_one = encode[1];
    int from_two = encode[2];
    int to_two = encode[3];
    int conn_type = encode[4];
    //printf("from_one %d to_one %d from_two %d to_two %d conn_type %d\n",
    //        from_one, to_one, from_two, to_two, conn_type);
    /* make change */
    int start_from_one = from_one==0?0:v_idx[3][from_one-1];
    int end_from_one = v_idx[3][from_one]; 
    int start_to_one = to_one==0?0:v_idx[3][to_one-1];
    int end_to_one = v_idx[3][to_one];
    int start_from_two = from_two==0?0:v_idx[3][from_two-1];
    int end_from_two = v_idx[3][from_two];
    int start_to_two = to_two==0?0:v_idx[3][to_two-1]; 
    int end_to_two = v_idx[3][to_two];
    if(conn_type == CONN_TYPE_ONE)
    {
        for(int i=start_from_one;i<end_from_one;i++)
        {
            if(e_idx[3][i] == to_one)
            {
                e_idx[3][i] = from_two;
                //printf("%d --> %d\n", from_one, from_two);
            }
        }
        for(int i=start_to_one;i<end_to_one;i++)
        {
            if(e_idx[3][i] == from_one)
            {
                e_idx[3][i] = to_two;
                //printf("%d --> %d\n", to_one, to_two);
            }
        }
        for(int i=start_from_two;i<end_from_two;i++)
        {
            if(e_idx[3][i] == to_two)
            {
                e_idx[3][i] = from_one;
                //printf("%d --> %d\n", from_two, from_one);
            }
        }
        for(int i=start_to_two;i<end_to_two;i++)
        {
            if(e_idx[3][i] == from_two)
            {
                e_idx[3][i] = to_one;
                //printf("%d --> %d\n", to_two, to_one);
            }
        }
    }
    else if(conn_type == CONN_TYPE_TWO)
    {
        for(int i=start_from_one;i<end_from_one;i++)
        {
            if(e_idx[3][i] == to_one)
            {
                e_idx[3][i] = to_two;
            }
        }
        for(int i=start_to_one;i<end_to_one;i++)
        {
            if(e_idx[3][i] == from_one)
            {
                e_idx[3][i] = from_two;
            }
        }
        for(int i=start_from_two;i<end_from_two;i++)
        {
            if(e_idx[3][i] == to_two)
            {
                e_idx[3][i] = to_one;
            }
        }
        for(int i=start_to_two;i<end_to_two;i++)
        {
            if(e_idx[3][i] == from_two)
            {
                e_idx[3][i] = from_one;
            }

        }
    }
    vector<int> a;
    for(int i = 0; i<5; i++)
    {
        a.push_back(encode[i]);
    }
    signature.push_back(a);
}

/**
 * resume the graph from shrinking
 *
 * @param[in]       encode      the encode to be shrunk
**/
void
InsMed::resume(int *encode)
{
    int from_one = encode[0];
    int to_one = encode[1];
    int from_two = encode[2];
    int to_two = encode[3];
    int conn_type = encode[4];
    int start_from_one = from_one==0?0:v_idx[3][from_one-1];
    int end_from_one = v_idx[3][from_one];
    int start_to_one = to_one==0?0:v_idx[3][to_one-1];
    int end_to_one = v_idx[3][to_one];
    int start_from_two = from_two==0?0:v_idx[3][from_two-1];
    int end_from_two = v_idx[3][from_two];
    int start_to_two = to_two==0?0:v_idx[3][to_two-1];
    int end_to_two = v_idx[3][to_two];
    if(conn_type == CONN_TYPE_ONE)
    {
        for(int i=start_from_one;i<end_from_one;i++)
        {
            if(e_idx[3][i] == from_two)
            {
                e_idx[3][i] = to_one;
            }
        }
        for(int i=start_to_one;i<end_to_one;i++)
        {
            if(e_idx[3][i] == to_two)
            {
                e_idx[3][i] = from_one;
            }
        }
        for(int i=start_from_two;i<end_from_two;i++)
        {
            if(e_idx[3][i] == from_one)
            {
                e_idx[3][i] = to_two;
            }
        }
        for(int i=start_to_two;i<end_to_two;i++)
        {
            if(e_idx[3][i] == to_one)
            {
                e_idx[3][i] = from_two;
            }
        }
    } 
    else if(conn_type == CONN_TYPE_TWO)
    {
        for(int i=start_from_one;i<end_from_one;i++)
        {
            if(e_idx[3][i] == to_two)
            {
                e_idx[3][i] = to_one;
            }
        }
        for(int i=start_to_one;i<end_to_one;i++)
        {
            if(e_idx[3][i] == from_two)
            {
                e_idx[3][i] = from_one;
            }
        }
        for(int i=start_from_two;i<end_from_two;i++)
        {
            if(e_idx[3][i] == to_one)
            {
                e_idx[3][i] = to_two;
            }
        }
        for(int i=start_to_two;i<end_to_two;i++)
        {
            if(e_idx[3][i] == from_one)
            {
                e_idx[3][i] = from_two;
            }
        }
    }
    signature.clear();
}

/**
 * add possible branch into current branch list
 *
 * @param[in]       from_one        one vertex from
 * @param[in]       to_one          one vertex to
 * @param[in]       from_two        two vertex from
 * @param[in]       to_two          two vertex to
**/
void InsMed::add_possible_branch(int from_one, int to_one, 
        int from_two, int to_two)
{
    vector<int> ft1;
    ft1.push_back(from_one);
    ft1.push_back(to_one);
    ft1.push_back(from_two);
    ft1.push_back(to_two);
    ft1.push_back(CONN_TYPE_ONE);
    vector<int> ft2;
    ft2.push_back(from_one);
    ft2.push_back(to_one);
    ft2.push_back(from_two);
    ft2.push_back(to_two);
    ft2.push_back(CONN_TYPE_TWO);
    ft.push_back(ft1);
    ft.push_back(ft2);
}

///////////////////////////////////////////////////////////////////
//
// DEBUG FUNCTIONS
//
//////////////////////////////////////////////////////////////////

void InsMed::print_csr()
{
    for(int c=0;c<3; c++)
    {
        for(int i=0; i<v_size; i++)
        {
            int start = i==0?0:v_idx[c][i-1];
            int end = v_idx[c][i];
            for(int j=start; j<end; j++)
            {
                printf("%d %d %d\n", i, e_idx[c][j], c+1);
            }
        }
    }
}

void InsMed::print_adj()
{
    for(int c=0; c<3; c++)
    {
        for(int i=0; i<v_size; i++)
        {
            for(int j=0; j<a_idx[c][i]; j++)
            {
                printf("%d %d %d\n", i, adj[c][i][j], c+1);
            }
        }
    }
}

void InsMed::vis_csr()
{
    FILE *writer = fopen("./med_csr.dot", "w");
    fprintf(writer, "graph{\n");
    for(int c=0;c<3; c++)
    {
        char *color;
        if(c==0)
        {
            color = (char*)"red";
        }
        if(c==1)
        {
            color = (char*)"blue";
        }
        if(c==2)
        {
            color = (char*)"green";
        }
        for(int i=0; i<v_size; i++)
        {
            int start = i==0?0:v_idx[c][i-1];
            int end = v_idx[c][i];
            for(int j=start; j<end; j++)
            {
                fprintf(writer, "%d -- %d [color=%s];\n", i, e_idx[c][j], color);
            }
        }
    }
    fprintf(writer, "}\n");
    fclose(writer);
}
void InsMed::vis_adj()
{
    FILE *writer = fopen("./med_adj.dot", "w");
    fprintf(writer, "graph{\n");
    for(int c=0; c<3; c++)
    {
        char *color=(char*)"";
        if(c==0)
        {
            color = (char*)"red";
        }
        if(c==1)
        {
            color = (char*)"blue";
        }
        if(c==2)
        {
            color = (char*)"green";
        }
        for(int i=0; i<v_size; i++)
        {
            for(int j=0; j<a_idx[c][i]; j++)
            {
                fprintf(writer, "%d -- %d [color=%s];\n", i, adj[c][i][j], color);
            }
        }
    }
    fprintf(writer, "}\n");
    fclose(writer);
}

void InsMed::vis_two(int c1, int c2)
{
    FILE *writer = fopen("./two_csr.dot", "w");
    fprintf(writer, "graph{\n");
    for(int i=0; i<v_size; i++)
    {
        int start = i==0?0:v_idx[c1][i-1];
        int end = v_idx[c1][i];
        for(int j=start; j<end; j++)
        {
            fprintf(writer, "%d -- %d [color=%s];\n", i, e_idx[c1][j], "red");
        }
    }
    for(int i=0; i<v_size; i++)
    {
        int start = i==0?0:v_idx[c2][i-1];
        int end = v_idx[c2][i];
        for(int j=start; j<end; j++)
        {
            fprintf(writer, "%d -- %d [color=%s];\n", i, e_idx[c2][j], "blue");
        }
    }
    fprintf(writer, "}\n");
}

void
InsMed::vis_encode(int *encode)
{
    FILE *writer = fopen("vis_encode.dot", "w");
    fprintf(writer, "graph{\n");
    for(int c=0; c<2; c++)
    {
        const char *c_str = c==0?"red":"blue";
        for(int i=0; i<encode[0]; i++)
        {
            int start = i==0 ? 0 : encode[3+c*encode[0]+i-1];
            int end = encode[3+c*encode[0]+i];
            start += ((3+2*encode[0]) + c*encode[c]);
            end += ((3+2*encode[0]) + c*encode[c]);
            for(int j=start; j<end; j++)
            {
                fprintf(writer, "%d -- %d [color=%s];\n", i, encode[j], c_str);
            }

        }
    }
    fprintf(writer, "}\n");
}

/**
 * Using BFS to detect connected components
 *
 * @param[in]       c_type          the data structure to store component type
 *
 * @return      num possible components
**/
int InsMed::detect_comps(int *c_type, int c1, int c2)
{
    int ret = 0;
    bool visited[v_size];
    int stack[v_size];
    int s_idx;
    int type_id = 0;
    /* if vertex is effective, set it as visited */
    for(int i=0; i<v_size; i++)
    {
        visited[i] = true;
        c_type[i] = -1;
    }
    for(int c=0; c<2; c++)
    {
        int cc = c == 0 ? c1 : c2;
        for(int i=0; i<v_size; i++)
        {
            int start = i==0 ? 0 : v_idx[cc][i-1];
            int end = v_idx[cc][i];
            if(end - start > 0)
            {
                visited[i] = false;
            }
        }
    }
    for(int i=0; i<v_size; i++)
    {
        if(visited[i])
        {
            continue;
        }
        /* let's start BFS for this component now */
        s_idx = 0;
        type_id++;
        bfs_visit(stack, visited, c_type, i, type_id, s_idx);
        while(s_idx > 0)
        {
            int v_now = stack[--s_idx];
            for(int color =0; color<2; color++)
            {
                int c = color == 0 ? c1 : c2;
                int start = v_now ==0 ? 0: v_idx[c][v_now-1];
                int end = v_idx[c][v_now];
                for(int j=start; j<end; j++)
                {
                    if(e_idx[c][j]==CAP)
                    {
                        continue;
                    }
                    else if(visited[e_idx[c][j]])
                    {
                        continue;
                    }
                    bfs_visit(stack, visited, c_type, e_idx[c][j], type_id, s_idx);
                }
            }
        }
    }
    /* return the number of components */
    ret = type_id;
    return ret;
}

/**
 *
 * set variables related to BFS visit, if a vertex is touched
 *
 * @param[in]       stack       stack to keep the BFS frontier
 * @param[in]       visited     the recode of vertex being visited
 * @param[in]       c_type      the connected component type of a vertex
 * @param[in]       vet         which vertex being touched
 * @param[in]       type_id     the connected component id
 * @param[in]       s_idx       stack idx
 *
 * **/
void InsMed::bfs_visit(int *stack, bool *visited, 
        int *c_type, int vet, int type_id, int &s_idx)
{
    if(vet != CAP)
    {
        stack[s_idx++] = vet;
        visited[vet] = true;
        c_type[vet] = type_id;
    }
}

/**
 * Retrive an encode that will be feed to insdis by given a specific component ID
 *
 * @param[in]       c_type          component type of each vertex
 * @param[out]      encode          data structure that keeps the encode to insdis
 * @param[in]       which_comp      which component to be retrieved           
 * @param[in]       c1              color one
 * @param[in]       c2              color two
 *
 * @return      size of the encode
 *
 **/
int InsMed::get_encode_by_component(int *c_type, int *encode, 
        int which_comp, int c1, int c2)
{
    encode[0] = v_size;
    encode[1] = encode[2] = 0;
    int begin = 3;

    for(int c=0; c<2; c++)
    {
        for(int i=0; i<v_size; i++)
        {
            encode[begin+c*v_size+i] = 0;
        }
    }

    /* get degree of vertices and their index */
    for(int c=0; c<2; c++)
    {
        int cc = c == 0 ? c1 : c2;
        int sum = 0;
        for(int i=0; i<v_size; i++)
        {
            if(c_type[i] != which_comp)
            {
                encode[c*v_size + i + begin] = sum;
                continue;
            }
            int start = i==0? 0: v_idx[cc][i-1];
            int end = v_idx[cc][i]; 
            encode[c*v_size + i + begin] = end - start;
            sum += encode[c*v_size + i +begin];
            encode[c*v_size + i + begin] = sum;
            encode[c+1] += (end - start); /* set num edges */
        }
    }
    /* put edges to encode */
    begin += v_size *2;
    for(int c=0; c<2; c++)
    {
        int cc = c == 0 ? c1 : c2;
        int tmp_idx = 0;
        for(int i=0; i<v_size; i++)
        {
            if(c_type[i] != which_comp)
            {
                continue;
            }
            int start = i==0?0:v_idx[cc][i-1];
            int end = v_idx[cc][i]; 
            for(int j=start; j<end; j++)
            {
                encode[begin + encode[c]*c + tmp_idx] = e_idx[cc][j];
                tmp_idx += 1;
            }
        }
    }
    int encode_size = 3 + encode[0]*2 + encode[1] + encode[2]; 
    return encode_size;
}

/**
 * detect if there is any conflict for example:
 * 1) there are two edges of one vertex connecting to the same vertex at the same time
 * 2) ....
 * TODO these are actually bugs, should solve it later
 *
 * @ return true or false
 *
**/
bool InsMed::detect_n_remove_conflict_edge()
{
    bool ret = false;
    /* allocate memory for mark */
    bool **mark = new bool*[4];
    for(int c=0; c<4; c++)
    {
        mark[c] = new bool[e_size[c]];
        for(int j=0; j<e_size[c]; j++)
        {
            mark[c][j] = false;
        }
    }

    /* detect case 1 */
    for(int c=0; c<4; c++)
    {
        for(int i=0; i<v_size; i++)
        {
            int start = i==0 ? 0 : v_idx[c][i-1];
            int end = v_idx[c][i];
            for(int j=start; j<end; j++)
            {
                for(int k=j+1; k<end; k++)
                {
                    if(e_idx[c][j] == e_idx[c][k])
                    {
                        mark[c][k] = true;
                        ret = true;
                    }
                }
            }
        }
    }

    /* allocate memory for temp data structure */
    int **tmp_v_idx = new int*[4];
    int **tmp_e_idx = new int*[4];
    for(int c=0; c<4; c++)
    {
        tmp_v_idx[c] = new int[v_size];
        tmp_e_idx[c] = new int[e_size[c]];
        for(int j=0; j<v_size; j++)
        {
            tmp_v_idx[c][j] = 0;
        }
        for(int j=0; j<e_size[c]; j++)
        {
            tmp_e_idx[c][j] = 0;
        }
    }

    /* get degrees of tmp_v_idx then its real value */
    int tmp_e_size[4];
    for(int c=0; c<4; c++)
    {
        int sum = 0;
        for(int i=0; i<v_size; i++)
        {
            int start = i==0 ? 0 : v_idx[c][i-1];
            int end = v_idx[c][i];
            for(int j=start; j<end; j++)
            {
                if(mark[c][j] == false)
                {
                    tmp_v_idx[c][i] += 1;
                }
            }
            sum += tmp_v_idx[c][i];
            tmp_v_idx[c][i] = sum;
        }
        tmp_e_size[c] = sum;
    }

    /* now copy e_idx to tmp_e_idx */
    for(int c=0; c<4; c++)
    {
        int tmp_pos = 0;
        for(int i=0; i<e_size[c]; i++)
        {
            if(mark[c][i] == false)
            {
                tmp_e_idx[c][tmp_pos] = e_idx[c][i];
                tmp_pos += 1;
            }
        }
    }

    /* now copy everything back */
    for(int c=0; c<4; c++)
    {
        for(int i=0; i<v_size; i++)
        {
            v_idx[c][i] = tmp_v_idx[c][i];
        }
        for(int i=0; i<tmp_e_size[c]; i++)
        {
            e_idx[c][i] = tmp_e_idx[c][i];
        }
        e_size[c] = tmp_e_size[c];
    }

    /* also change the representation of adj */
    for(int c=0; c<4; c++)
    {
        from_csr_to_adj(c);
    }

    /* delete memory for mark */
    for(int c=0; c<4; c++)
    {
        delete [] mark[c];
        delete [] tmp_v_idx[c];
        delete [] tmp_e_idx[c];
    }
    delete [] mark;
    delete [] tmp_v_idx;
    delete [] tmp_e_idx;
    return ret;
}

/**
 * check if a connected comonent represented by encode is single
 * component or not
 *
 * @param[in]       encode          encode to represent a graph
 *
 * @return      true or false
 *
**/
bool InsMed::check_if_single_cc(int *encode)
{
    bool ret = true;
    for(int c=0; c<2; c++)
    {
        for(int m=0; m<encode[0]; m++)
        {
            int start = m==0 ? 0 : encode[3+encode[0]*c + m-1];
            int end = encode[3+encode[0]*c + m];
            start += 3+2*encode[0] + c*encode[1];
            end += 3+2*encode[0] + c*encode[1];
            if(end-start>1)
            {
                ret = false;
            }
        }
    }
    return ret;
}

void InsMed::rename_encode_by_comp(int *encode)
{
    /* allocate the memory for the algorithm */
    int old_encode_size = 3+encode[0]*2+encode[1]+encode[2];
    int tmp_encode[old_encode_size];
    int name_map[v_size];
    bool exist[v_size];
    for(int i=0; i<v_size; i++)
    {
        name_map[i] = -1;
        exist[i] = false;
    }

    int begin = 3;
    /* find which vertex exists */
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<v_size; i++)
        {
            if(i==0  && encode[begin+c*encode[0]+i]>0)
            {
                exist[i] = true;
            }
            else if(encode[begin+c*encode[0]+i]-encode[begin+c*encode[0]+i-1]>0) 
            {
                exist[i] = true;
            }

        }
    }
    /* get the map */
    int total = 0;
    for(int i=0; i<v_size; i++)
    {
        if(exist[i] == true)
        {
            name_map[i] = total;
            total += 1;
        }
    }
    tmp_encode[0] = total;
    tmp_encode[1] = encode[1];
    tmp_encode[2] = encode[2];

    /*rearrange v_idx*/
    for(int c=0; c<2; c++)
    {
        int tmp_idx = 0;
        for(int i=0; i<v_size; i++)
        {
            int start = i==0?0:encode[begin+c*v_size+i-1];
            int end = encode[begin+c*v_size+i];
            if(end-start>0)
            {
                tmp_encode[begin+c*tmp_encode[0]+tmp_idx] = encode[begin+c*v_size+i];
                tmp_idx += 1;
            }
        }
    }
    /*rearrange e_idx*/
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<encode[c+1]; i++)
        {
            if(encode[begin+encode[0]*2+c*encode[c]+i] == CAP)
            {
                tmp_encode[begin+tmp_encode[0]*2+c*tmp_encode[c]+i] = 
                    encode[begin+encode[0]*2+c*encode[c]+i];
                //printf("%d:%d ", begin+tmp_encode[0]*2+c*tmp_encode[c]+i,
                //        encode[begin+encode[0]*2+c*encode[c]+i]);
            }
            else
            {
                tmp_encode[begin+tmp_encode[0]*2+c*tmp_encode[c]+i] = 
                    name_map[encode[begin+encode[0]*2+c*encode[c]+i]];
                //printf("%d:%d ", begin+tmp_encode[0]*2+c*tmp_encode[c]+i,
                //        name_map[encode[begin+encode[0]*2+c*encode[c]+i]]);
            }
        }
        //printf("\n");
    }
    /* copy back */
    for(int i=0; i<begin+tmp_encode[0]*2+tmp_encode[1]+tmp_encode[2]; i++)
    {
        encode[i] = tmp_encode[i];
    }
}

void InsMed::print_c_type(int *c_type)
{
    for(int i=0; i<v_size; i++)
    {
        if(c_type[i] != -1)
        {
            printf("%d:%d ", i, c_type[i]);
        }
    }
    printf("\n");
}

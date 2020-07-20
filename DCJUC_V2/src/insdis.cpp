#include "insdis.h"
#include "utils.h"
#include <string>
#include <iostream>
#include <cstring>
#include <string>

using namespace std;

extern bool control;
extern int g_count;
extern int which_med;

/******************************
 * these are public functions
 * for branching operations
 * ****************************/
void
p_encode(int *encode);
void
vis_encode(int *encode);
string map[7] = {"","pi_open","gamma_open","cycle","end","none","out"};


/**
 * constructor 
 *
 * @param[in]       f           input of the graph file
 * @param[in]       is_exe      indicate if using exemplar distance
 * @param[in]       fd          input of the dictionary file
 * @param[out]      fbi         output of the bijection file
 * @param[in]       is_cal_bij  indicate if calculate bijection or not
 *
 * TODO upper and lower bound are not the smae
 * **/
InsDis::InsDis(const char *f, bool is_exe, const char *lf, const char* fd, const char* fbi, const bool *is_cal_bij) : 
    Instance(lf)
{
    if(f != NULL)
    {
        file = new char[OPTKIT_FILE_SIZE];
        snprintf(file, OPTKIT_FILE_SIZE, "%s", f);
    }
    /* some member variables */
    comb_size = 0;
    v_size = 0;
    /* by default not calculate bijection */
    is_calculate_bijection = false;
    /* some basic info */
    is_exem = is_exe;
    is_ub = true;
    vet = OPTKIT_NULL;
    branch_id = OPTKIT_NULL;
    upper_bound = OPTKIT_NULL;
    lower_bound = OPTKIT_NULL;
    num_branch = OPTKIT_NULL;
    /* init graph and init some related data structures */
    init_graph();
    /* turn csr format into adj format the reason why we need csr is */
    /* that csr can help to compress the graph */
    from_csr_to_adj();
    /* at this time, there is no branch, so branch_id is -1 */
    branch_id = OPTKIT_NULL;
    vet = OPTKIT_NULL;
    compute_bound();
    get_num_elem();
    count=0;

    if(is_cal_bij != NULL)
    {
        is_calculate_bijection = is_cal_bij;

        bij_file = (char*)malloc(sizeof(char)*OPTKIT_FILE_SIZE);
        memcpy(bij_file, fbi, strlen(fbi));

        bijection = (int**)malloc(sizeof(int*)*2);
        for(int i=0; i<2; i++)
        {
            bijection[i] = (int*)malloc(sizeof(int)*v_size*4);
        }
        bij_idx = (int*)malloc(sizeof(int)*2);
        bij_idx[0] = 0;
        bij_idx[1] = 0;

        dict_file = (char*)malloc(sizeof(char)*OPTKIT_FILE_SIZE);
        memcpy(dict_file, fd, strlen(fd));
        dict = (int*)malloc(v_size*4); 

        /* read dictionary files */
        read_dict();
    }
}


/**
 * destructor to free all the things
 * **/
InsDis::~InsDis()
{
    /* free files */
    free(file);
    free(bij_file);
    free(dict_file);

    /* free CSRs */
    free(v_idx[0]);
    free(v_idx[1]);
    free(e_idx[0]);
    free(e_idx[1]);
    free(v_idx);
    free(e_idx);

    /* free ADJs */
    for(int c=0; c<2; c++)
    {
        free(adj_idx[c]);
        free(idx_tmp[c]);
        for(int i=0; i<v_size; i++)
        {
            free(adj[c][i]);
            free(a_tmp[c][i]);
        }
        free(adj[c]);
    }
    free(adj);
    free(adj_idx);
    free(a_tmp);
    free(idx_tmp);

    /* free adj for real computation */
    free(a[0]);
    free(a[1]);
    free(vet_type);

    /* clear vector */
    comb.clear();

    free(dict);
    free(bijection[0]);
    free(bijection[1]);
    free(bijection);
    free(bij_idx);
}


/**
 * set up branch id and calculate upper and lower bound
 *
 * @param[in]       which_branch    the branch_id
 **/
    void 
InsDis::to_branch(int which_branch)
{
    /* do nothing here */
    branch_id = which_branch;
    /* then update upper and lower bound */
    compute_bound();
}

/**
 * resume to original state
 *
**/
void 
InsDis::from_branch()
{
    //no need to do anything
    avail_id = v_size;
}

/**
 * compute the upper and lower bound
**/
    void 
InsDis::compute_bound()
{
    /* get adjs first */
    if(is_exem==true)
    {
        get_exemplar_adj(vet, branch_id);
    }
    else
    {
        get_matching_adj(vet, branch_id);
    }

    /* compute the lower bound */
    a_size = is_exem?v_size:num_v_after_relabel;
    set_vet_type(a_size);
    lower_bound = compute_indel_dis(); 
    
    /* then compute the upper bound */
    get_deletion_adj(vet, branch_id, is_exem);
    set_vet_type(v_size+comb_size);
    upper_bound = compute_indel_dis() - pair_num; 

    /* some patches */
    if(upper_bound<lower_bound) 
    {
        upper_bound = lower_bound;
    }
    if(upper_bound -1 == lower_bound)
    {
        int encode[30000];
        get_encode(encode);
        bool is_no_branch = detect_no_branch(encode);
        if(is_no_branch == true)
        {
            /* round lower bound to upper bound */
            lower_bound = upper_bound;
        }
    }
}

/**
 * calculate the combinations 
 * and get number of branches
**/
int 
InsDis::get_num_branches()
{
    /* first find the first available branch */
    for(int i=0; i<v_size; i++)
    {
        printf("deg %d:(%d-%d)\n", i, adj_idx[0][i], adj_idx[1][i]);
        if(adj_idx[0][i]>1 || adj_idx[1][i]>1)
        {
            vet = i;
            break;
        }
    }

    /* then compute the possible combinations */
    if(is_exem == true)
    {
        get_comb_exem();
    }
    else
    {
        get_comb_matc();
    }
    num_branch = comb.size();
    return num_branch;
}

/**
 * from adj to encode, for the storage of the graph
 *
 * @param[in]       encode          translate current state into encode
**/
int 
InsDis::get_encode(int *encode)
{
    /* turn adj into CSR first */
    /* fill into tmp adj */
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<v_size+comb_size; i++) 
        {
            idx_tmp[c][i] =0;
        }
    }
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<v_size+comb_size; i++)
        {
            /* simply copy, if not vet and its neighbors */
            if(check_if_in_to(i, c)==false && i != vet)
            {
                /* for each of this vet's neighbor */
                for(int j=0; j<adj_idx[c][i]; j++)
                {
                    /* not vet's neighbor */
                    if(check_if_in_to(adj[c][i][j], c)==false)
                    {
                        int from = i;
                        int to = adj[c][i][j];
                        int pos = idx_tmp[c][from]++;
                        a_tmp[c][from][pos] = to;
                    }
                }
            }
            else
            {
                /* including vet and its to vertices */
                idx_tmp[c][i] = 0; 
            }
        }
    }
    /* cope the rest */
    /* do if exemplar */
    if(is_exem==true && vet != OPTKIT_NULL)
    {
        for(int c=0; c<2; c++)
        {
            int from = vet;
            /* comb[b][0] is the first, comb[b][1] is the second */
            int to = adj[c][from][comb[branch_id][c]];
            a_tmp[c][from][0] = to; 
            a_tmp[c][to][0] = from; 
            idx_tmp[c][from]  = 1;
            idx_tmp[c][to]  = 1;
            for(int j=0; j<adj_idx[c][to]; j++)
            {
                if(adj[c][to][j] != from)
                {
                    int pos = idx_tmp[c][to]++;
                    int toto = adj[c][to][j];
                    int toto_pos = idx_tmp[c][toto]++;
                    a_tmp[c][to][pos] = toto;
                    a_tmp[c][toto][toto_pos] = to;
                }
            }
            /* deletions only, and copy to tmp */
            for(int i=0; i<adj_idx[c][vet]; i++)
            {
                if(adj[c][vet][i] != to)
                {
                    int another = adj[c][vet][i];
                    for(int j=0, k=0; j<adj_idx[c][another]; j++)
                    {
                        if(adj[c][another][j] != from)
                        {
                            int toto = adj[c][another][j];
                            int toto_pos = idx_tmp[c][toto]++;
                            a_tmp[c][another][k++] = toto;
                            a_tmp[c][toto][toto_pos] = another;
                        }
                    }
                    idx_tmp[c][another] = adj_idx[c][another]-1;
                } 
            }
        }
    }
    /* do if matching */
    else if(is_exem == false && vet != OPTKIT_NULL)
    {
        /* color 1 */
        int avail = avail_id;
        int from = vet;
        int to = adj[0][vet][0];
        a_tmp[0][from][0] = to;
        a_tmp[0][to][0] = from;

        /* mapping */
        if(is_calculate_bijection)
        {
            bijection[0][bij_idx[0]++] = from; 
            bijection[0][bij_idx[0]++] = to; 
        }

        idx_tmp[0][vet] = 1;
        idx_tmp[0][adj[0][vet][0]]=1;
        for(int i=1; i<adj_idx[0][vet]; i++)
        {
            int from = avail+i-1;    
            int to = adj[0][vet][i];

            /* mapping */
            if(is_calculate_bijection)
            {
                bijection[0][bij_idx[0]++] = vet; 
                bijection[0][bij_idx[0]++] = to; 
            }

            a_tmp[0][from][0] = to;
            a_tmp[0][to][idx_tmp[0][to]++] = from;
            idx_tmp[0][from] = 1;
        }
        /* color 2 */
        for(unsigned int j=0; j<comb[branch_id].size(); j++)
        {
            int from, to;
            int id = comb[branch_id][j];    
            if(id==0)
            {
                from = vet;
            }
            else
            {
                from = avail + id -1;
            }
            to = adj[1][vet][j];

            /* mapping */
            if(is_calculate_bijection)
            {
                bijection[1][bij_idx[1]++] = vet; 
                bijection[1][bij_idx[1]++] = to; 
            }

            a_tmp[1][from][0] = to;
            a_tmp[1][to][idx_tmp[1][to]] = from;
            idx_tmp[1][from] = 1;
            idx_tmp[1][to] += 1;
        }
        //cope to
        for(int c=0; c<2; c++)
        {
            for(int i=0; i<adj_idx[c][vet]; i++)
            {
                int to = adj[c][vet][i];
                for(int j=0; j<adj_idx[c][to]; j++)
                {
                    if(adj[c][to][j] != vet)
                    {
                        /* if toto is also a adj of vet */
                        if(check_if_in_to(adj[c][to][j], c) && to < adj[c][to][j])
                        {
                            continue; // just skip
                        }
                        /* otherwise */
                        else
                        {
                            int toto = adj[c][to][j];
                            int to_i = idx_tmp[c][to]++;
                            int toto_i = idx_tmp[c][toto]++;
                            a_tmp[c][to][to_i] = toto;
                            a_tmp[c][toto][toto_i] = to;
                        }
                    }
                }
            }
        }
    }
    /* get v_idx first */
    /* color 0 */
    int new_v_size; 
    if(is_exem == true)
    {
        /* this is dependent on which algorithm you are using */
        new_v_size = v_size; 
    }
    else
    {
        /* this is not correct */
        new_v_size = v_size + comb_size; 
    }
    encode[0] = new_v_size;
    /* encode v_idx */
    int start = 3;
    for(int c=0; c<2; c++){
        int sum =0;
        for(int i=0; i<new_v_size; i++){
            sum += idx_tmp[c][i]; 
            encode[start+i] = sum;
        }
        encode[c+1] = sum;
        start += new_v_size;
    }
    /* burn in e_idx then */
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<new_v_size; i++)
        {
            for(int j=0; j<idx_tmp[c][i]; j++)
            {
                encode[start++] = a_tmp[c][i][j];
            }
        }
    }
    /* append the bijection information at the end of the code */
    if(is_calculate_bijection)
    {
        start = append_bijection_encode(encode, start);
    }
    return start;
}

/**
 * encode is stored in CSR format then transformed into adj format
 *
 * @param[in]       encode      code to be interpreted 
 * @param[in]       size        size of the encode 
**/
void 
InsDis::to_encode(int *encode, int size)
{
    /* the first entry of encode is e_size[0] the second is e_size[1] */
    v_size = encode[0];
    e_size[0] = encode[1];
    e_size[1] = encode[2];
    int g_start = 0;
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<v_size; i++)
        {
            adj_idx[c][i] = 0;
            for(int j=0; j<v_size; j++)
            {
                adj[c][i][j] = -1;
            }
        }
    }
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<encode[0]; i++)
        {
            int start = i==0?0:encode[3+c*encode[0]+i-1];
            int end = encode[3+c*encode[0]+i];
            start += ((3+2*encode[0]) + c*encode[c]);
            end += ((3+2*encode[0]) + c*encode[c]);
            for(int j=start; j<end; j++)
            {
                adj[c][i][adj_idx[c][i]++] = encode[j];
            }
            g_start = end;
        }
    }
    avail_id = v_size;
    comb_size = 0;
    /* get the bijection information */
    if(is_calculate_bijection)
    {
        get_bijection_encode(encode, g_start);
    }
}

/**
 * print the content of the encode, and print the bijection information
**/
void 
InsDis::print_encode()
{
    int encode[30000];
    get_encode(encode);
    printf("num_v %d num_e1 %d num_e2 %d\n", encode[0], encode[1], encode[2]);
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<encode[0]; i++)
        {
            int start = i==0?0:encode[3+c*encode[0]+i-1];
            int end = encode[3+c*encode[0]+i];
            start += ((3+2*encode[0]) + c*encode[c+1]);
            end += ((3+2*encode[0]) + c*encode[c+1]);
            for(int j=start; j<end; j++)
            {
                printf("%d:%d ", i, encode[j]);
            }

        }
        printf("\n");
    }
    /* at the end append the bijection info */
    if(is_calculate_bijection)
    {
        printf("come here! bij_file %s\n", bij_file);
        FILE* writer = fopen(bij_file, "a");
        if(writer == NULL)
        {
            printf("the file %s does not exists!\n", bij_file);
            exit(1);
        }
        for(int i=0; i<bij_idx[i]; i+=2)
        {
            fprintf(writer, "%d %d:%d %d\n", 
                    dict[bijection[0][i]], dict[bijection[0][i+1]],
                    dict[bijection[1][i]], dict[bijection[1][i+1]]);
            printf("%d %d:%d %d\n", 
                    dict[bijection[0][i]], dict[bijection[0][i+1]],
                    dict[bijection[1][i]], dict[bijection[1][i+1]]);
        }
        fclose(writer);
    }
}

/**
 * copy code from another instance to this instance
 *
 * @param[in]       encode      the encode tobe copied
 * @param[in]       size        size of the encode
 **/
    void 
InsDis::copy_code(int *encode, int size)
{
    to_encode(encode, size);
}

/**
 * return some private inforation that programer needs 
 *
 * @return      any information that programmer interested in
 **/
    int 
InsDis::get_value()
{
    return upper_bound;
}


/******************************
 * these are private functions
 * for graph operations
 * ****************************/

/**
 * transfer csr to adj for the convinience of comuting
 **/
    void 
InsDis::from_csr_to_adj()
{
    int extra_num_v = 0;
    for (int c=0; c<2; c++)
    {
        for (int i=0; i<v_size; i++)
        {
            adj_idx[c][i] = 0;
            /* start copying */
            int start = i==0?0:v_idx[c][i-1];
            int end = v_idx[c][i];
            extra_num_v += (end-start)==0?0:(end-start)-1;
            for (int j=start; j<end; j++)
            {
                adj[c][i][adj_idx[c][i]++] = e_idx[c][j];
            }
        }
    }
    num_v_after_relabel = v_size + extra_num_v;
}


/**
 * use deletion to calculate lower bound
 **/
    void 
InsDis::get_deletion_adj(int vet, int comb_id, bool is_exem)
{
    a_size = is_exem==true?v_size:num_v_after_relabel;
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<a_size; i++)
        {
            a[c][i] = OPTKIT_NULL;
        }
    }
    for(int i=0; i<v_size; i++)
    {
        if(i==vet)
        {
            /* no problem with this */
            if(is_exem==true)
            {
                for(int c=0; c<2; c++)
                {
                    int from = i;
                    int to = adj[c][from][comb[comb_id][c]];
                    /* another patch */
                    if(a[c][to] != OPTKIT_NULL) 
                    {
                        int toto = a[c][to];
                        if(toto != CUP)
                        {
                            a[c][toto] = OPTKIT_NULL;
                        }
                    }
                    a[c][from] = to; 
                    a[c][to] = from; 
                }
            }
            else
            {
                /* color 0, expand directly */
                int avail = avail_id;
                int from = i;
                int to = adj[0][from][0];
                a[0][from] = to;
                a[0][to] = from;
                for(int j=1; j<adj_idx[0][i]; j++)
                {
                    int from = avail++; 
                    int to = adj[0][i][j];
                    a[0][from] = to;
                    a[0][to] = from;
                }
                /* color 1, expand according to mapping */
                for(unsigned int j=0; j<comb[branch_id].size(); j++)
                {
                    int id = comb[branch_id][j];    
                    if(id==0)
                    {
                        from = i;
                    }
                    else
                    {
                        from = avail_id + id -1;
                    } 
                    to = adj[1][i][j];
                    a[1][from] = to;
                    a[1][to] = from;
                }
            }
        }
        else
        {
            for(int c=0; c<2; c++)
            {
                int c_another = c==0?1:0;
                /* just delete */
                if((adj_idx[c][i] > 1 || adj_idx[c_another][i] > 1) 
                        && a[c][i]==OPTKIT_NULL)
                {
                    a[c][i] = OPTKIT_NULL;
                }
                /* just copy */
                else if(adj_idx[c][i] == 1 && a[c][i]==OPTKIT_NULL)
                { 
                    /* because it might be zero */
                    int from = i;
                    int to = adj[c][from][0];
                    if(to != CAP && adj_idx[c][to]==1 && adj_idx[c_another][to]==1)
                    {
                        a[c][from] = to; 
                        a[c][to] = from; 
                    }
                    else
                    {
                        a[c][i] = CUP;
                    }
                }
            }
        }
    }
}

/**
 * turning into a adj that is computable of the vet into possible combination
 * to compute exemplar distance
 *
 * @param[in]       vet         which vet tobe considered
 * @param[in]       comb_id     indicate which combination 
**/
void 
InsDis::get_exemplar_adj(int vet, int comb_id)
{
    for(int c=0; c<2; c++)
    {
        /* before execution init a and a_size and vet_type */
        for(int i=0; i<v_size; i++)
        {
            a[c][i] = OPTKIT_NULL;
        }
        for(int i=0; i<v_size; i++)
        {
            int from, to;
            /* correct! good job */
            if(i==vet)
            {
                from = i;
                to = adj[c][from][comb[comb_id][c]];
                if(to == OPTKIT_NULL)
                {
                    printf("vet %d from %d comb[comb_id][c] %d comb_id %d \n", vet, from, comb[comb_id][c], comb_id);
                    vis_adj();
                    exit(1);
                }
                /* already assigned, should roll back this operation */
                if(to != CAP && a[c][to] != OPTKIT_NULL ) 
                {
                    int toto = a[c][to];
                    a[c][toto] = OPTKIT_NULL;
                }
                a[c][from] = to; 
                a[c][to] = from; 
                /* calculate the bijection */
                if(is_calculate_bijection)
                {
                    bijection[c][bij_idx[c]++] = from; 
                    bijection[c][bij_idx[c]++] = to; 
                }
            }
            else if(a[c][i]==OPTKIT_NULL && adj_idx[c][i]>0 
                && check_if_in_to(i, c)==false)
            {
                from = i;
                if(adj_idx[c][from] < 0)
                {
                    continue;
                }
                /* random position */
                /* it does'nt matter if degree is 1 or more than 1 */
                int pos = rand()%adj_idx[c][from];
                to = adj[c][from][pos];
                if(CAP == to)
                {
                    ;
                }
                else if(a[c][to] != OPTKIT_NULL)
                {
                    ;
                }
                else if(to == vet)
                {
                    ;
                }
                else
                {
                    a[c][from] = to; 
                    a[c][to] = from; 
                }
            }
        }
    }
}

/**
 * turning into a adj that is computable of the vet into possible combination
 * to compute matching distance
 *
 * @param[in]       vet         which vet tobe considered
 * @param[in]       comb_id     indicate which combination 
 *
**/
void
InsDis::get_matching_adj(int vet, int comb_id)
{
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<num_v_after_relabel; i++)
        {
            a[c][i] = OPTKIT_NULL;
        }
    }
    int avail = avail_id;
    for(int i=0; i<v_size; i++){
        if(i==vet ){
            /* color 0, expand directly */
            int from = i;
            int to = adj[0][from][0];
            a[0][from] = to;
            a[0][to] = from;
            for(int j=1; j<adj_idx[0][i]; j++)
            {
                int from = avail+j-1; 
                int to = adj[0][i][j];
                a[0][from] = to;
                a[0][to] = from;
            }
            /* color 1, expand according to mapping */
            for(unsigned int j=0; j<comb[branch_id].size(); j++)
            {
                int id = comb[branch_id][j];    
                if(id==0)
                {
                    from = i;
                }
                else
                {
                    from = avail + id -1;
                }
                to = adj[1][i][j];
                a[1][from] = to;
                a[1][to] = from;
            }
            avail += adj_idx[0][i]>adj_idx[1][i]?
                (adj_idx[0][i] - 1):
                (adj_idx[1][i] - 1);
        }
        /* for the rest direct matching */
        else
        {
            int ava[2];
            for(int c=0; c<2; c++)
            {
                ava[c] = avail;
                if(check_if_in_to(i, c)==true)
                {
                    continue;
                }
                if(a[c][i] == OPTKIT_NULL && adj_idx[c][i]>0)
                {
                    for(int j=0; j<adj_idx[c][i]; j++)
                    {
                        int from  = i;
                        int to = adj[c][from][j];
                        /* in case that two duplicated vertices */
                        /* are connected, only one stretch is needed */
                        if(check_if_in_to(to, c)==false && 
                            a[c][to] == OPTKIT_NULL && 
                            (adj_idx[c][to]==1 || to>from))
                        {
                            a[c][from] = to;
                            a[c][to] = from;
                            break;
                        }
                    }
                }
                if(adj_idx[c][i] > 1)
                {
                    /* if degree is equal to 1 this step will be skipped */
                    for(int j=0; j<adj_idx[c][i]; j++)
                    {
                        if(adj[c][i][j]==a[c][i])
                        {
                            continue;
                        }
                        int from = ava[c];
                        int to = adj[c][i][j];
                        if(a[c][to] ==OPTKIT_NULL && a[c][from]==OPTKIT_NULL)
                        {
                            a[c][from] = to;
                            a[c][to] = from;
                            ava[c] += 1;
                        }
                    }
                }
            }
            avail = ava[0]>ava[1]?ava[0]:ava[1];
        }
    }
    /* can we add some logic here to update the missing edges */
    /* for example, vertices in 'to', but itself is irregular vertex? */

    int ava[2];
    for(int c=0; c<2; c++)
    {
        ava[c] = avail;
        for(int i=0; i<v_size; i++)
        {
            if(check_if_in_to(i, c)==true && adj_idx[c][i]>1)
            {
                for(int j=0; j<adj_idx[c][i]; j++)
                {
                    if(adj[c][i][j] != vet)
                    {
                        int from = ava[c]++;
                        int to = adj[c][i][j];
                        if(a[c][from]==OPTKIT_NULL && a[c][to]==OPTKIT_NULL)
                        {
                            a[c][from] = to;
                            a[c][to] = from;
                        }
                    }
                }
            }
        }
    }
}

/**
 * get all exemplar combinations
 * it doesn't matter if deg1 != deg2
 * exemplar comb is easy to compute
**/
void 
InsDis::get_comb_exem()
{
    int deg1 = adj_idx[0][vet];
    int deg2 = adj_idx[1][vet];
    comb.clear();
    for(int i=0; i<deg1; i++)
    {
        for(int j=0; j<deg2; j++)
        {
            vector<int> tmp;
            tmp.push_back(i);
            tmp.push_back(j);
            comb.push_back(tmp);
        }
    }
    for(int c=0; c<2; c++)
    {
        for(int i=v_size; i<v_size+comb_size; i++) 
        {
            adj_idx[c][i] = 0;
        }
    }
}

/**
 * get all matching combinations
 * by using BFS
 * in this algorithm we should distinguish between larger and smaller deg edges
 * **/
void 
InsDis::get_comb_matc(){
    /* some preprocess steps */
    int deg1 = v_idx[0][vet]-(vet==0?v_idx[0][0]:v_idx[0][vet-1]);
    int deg2 = v_idx[1][vet]-(vet==0?v_idx[1][0]:v_idx[1][vet-1]);
    int smaller = deg1<deg2?deg1:deg2;
    int larger = deg1<deg2?deg2:deg1;
    /* it really matters which is smaller */
    comb.clear();
    vector<vector<int > > queue;
    for(int i=0; i<larger; i++)
    {
        vector<int> tmp;
        tmp.push_back(i);
        queue.push_back(tmp);
    }
    /* start bfs */
    while(queue.size()!=0)
    {
        vector<int> item = queue.back();
        queue.pop_back();
        /* check final */
        if((signed)item.size()==smaller)
        {
            comb.push_back(item);
        }
        else
        {
            for(int i=0; i<larger; i++)
            {
                bool find = false;
                /* check if item already in */
                for(int j=0; j<(signed)item.size(); j++)
                {
                    if(i==item[j])
                    {
                        find = true;
                        break;
                    }
                }
                if(find==false)
                {
                    item.push_back(i);
                    queue.push_back(item);
                }
            }
        }
    }
    /* post process to get id map of each */
    if(deg1==larger)
    {
        ;
    }
    else
    {
        for(unsigned int i=0; i<comb.size(); i++)
        {
            int tmp[larger];
            int avail = smaller;
            for(int j=0; j<larger; j++) 
            {
                tmp[j] = -1;
            }
            for(unsigned int j=0; j<comb[i].size(); j++)
            {
                tmp[comb[i][j]] = j;
            }
            for(int j=0; j<larger; j++) 
            {
                if(tmp[j] == -1) 
                {
                    tmp[j] = avail++;
                }
            }
            comb[i].clear();
            for(int j=0; j<larger; j++) 
            {
                comb[i].push_back(tmp[j]);
            }
        }
    } 
    comb_size = larger-1;
    for(int c=0; c<2; c++)
    {
        for(int i=v_size; i<v_size+comb_size; i++) 
        {
            adj_idx[c][i] = 0;
        }
    }
}

/**
 * Initiate the Ins_Dis first
**/
void
InsDis::init_graph()
{
    /* scan the graph to get the basic information */
    scan_graph();
    read_graph();
    avail_id = v_size;
}


/**
 * read the dictionary file into local storage
 *
 * @param[in]       dict        data structure to keep the info
**/
void 
InsDis::read_dict()
{
    FILE *reader = fopen(dict_file, "r");
    if(reader == NULL)
    {
        printf("error file %s does not exists! \n", dict_file);
        exit(1);
    }
    int num_line = v_size;
    
    /* read the content */
    int new_v = -1;
    int old_v = -1;
    for(int i=0; i<num_line; i++)
    {
        fscanf(reader, "%d %d\n", &new_v, &old_v);
        dict[new_v] = old_v;
    }

    /* close file */
    fclose(reader);
}

/**
 * scan the graph to get some basic information
**/
void
InsDis::scan_graph(){
    /* declare variables */
    FILE *stream;
    char file_buf[OPTKIT_FILE_SIZE];
    int v_num;
    int e_num;
    int g_num;
    int v_id=0;
    int v_to=0;
    int color=0;
    int i,j,sum=0;

    /* read the head information, and allocate according variables */
    sprintf(file_buf, "%s", file);
    if((stream=fopen(file_buf,"r"))==NULL)
    {
        printf("the file %s you input does not exist!\n", file_buf);
        exit(1);
    }
    if(fscanf(stream, "%d %d %d %d\n", &v_num, &g_num, &e_num, &pair_num)==EOF)
    {
        printf("error here\n");
    }
    v_size = v_num;
    v_idx = (int**)malloc(sizeof(int*)*2);
    e_idx = (int**)malloc(sizeof(int*)*2);
    for(i=0;i<2;i++)
    {
        v_idx[i] = (int*)malloc(sizeof(int)*v_size*4);
        for(int j=0; j<v_size*4; j++)
        {
            v_idx[i][j]  = 0;
        }
        e_idx[i] = (int*)malloc(sizeof(int)*e_num*4);
        for(int j=0; j<e_num*4; j++)
        {
            e_idx[i][j]  = 0;
        }
        e_size[i] = 0;
    }
    /* scan the real content */
    for(i=0;i<e_num;i++)
    {
        if(fscanf(stream, "%d %d %d\n",&v_id , &v_to, &color)==EOF)
        {
            printf("error here\n");
        }
        e_size[color-1] += 1;
        v_idx[color-1][v_id] +=1;
    }
    for(i=0;i<2;i++)
    {
        sum=0;
        for(j=0;j<v_num;j++)
        {
            sum += v_idx[i][j];
            v_idx[i][j]=sum;
        }
    }
    /* allocate adjs */
    adj = (int***)malloc(sizeof(int**)*2);
    adj_idx = (int**)malloc(sizeof(int*)*2);
    a_tmp = (int***)malloc(sizeof(int**)*2);
    idx_tmp = (int**)malloc(sizeof(int*)*2);
    a = (int**)malloc(sizeof(int*)*2);
    vet_type = (int*)malloc(sizeof(int)*(v_size)*4);
    for (int c=0; c<2; c++)
    {
        adj[c] = (int**)malloc(sizeof(int*)*(v_size)*4);
        adj_idx[c] = (int*)malloc(sizeof(int)*(v_size)*4);
        a_tmp[c] = (int**)malloc(sizeof(int*)*(v_size)*4);
        idx_tmp[c] = (int*)malloc(sizeof(int)*(v_size)*4);
        a[c] = (int*)malloc(sizeof(int)*(v_size)*4);
        for (int i=0; i<(v_size)*4; i++)
        {
            /* can't have duplication of more than 100 */
            adj[c][i] = (int*)malloc(sizeof(int)*100);
            a_tmp[c][i] = (int*)malloc(sizeof(int)*100);
            a[c][i] = -1;
        }
    }
    fclose(stream);
}

/**
 * read the graph into local data structure
**/
void
InsDis::read_graph()
{
    FILE *stream;
    char file_buf[OPTKIT_FILE_SIZE];
    int v_num;
    int e_num;
    int g_num;
    int v_id=0;
    int v_to=0;
    int color=0;
    int i,j;
    /* read basic information */
    sprintf(file_buf, "%s", file);
    if((stream=fopen(file_buf,"r"))==NULL)
    {
        printf("the file %s you input does not exist!\n", file_buf);
        exit(1);
    }
    if(fscanf(stream, "%d %d %d %d\n", &v_num, &g_num, &e_num, &pair_num)==EOF)
    {
        printf("error here\n");
    }
    /* this is for the start of the different postitions */
    int **idx = (int**)malloc(sizeof(int*)*2);
    for(i=0;i<2;i++)
    {
        idx[i]=(int*)malloc(sizeof(int)*v_size);
        for(j=0;j<v_size;j++)
        {
            idx[i][j] = j==0?0:v_idx[i][j-1];
        }
    }
    /* real read part */
    for(i=0;i<e_num;i++)
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
    for(i=0; i<2; i++)
    {
        free(idx[i]);
    }
    free(idx);
    fclose(stream);
}

/**
 * also need to have the function to compute
 * the types of different veritces
 *
 * @param[in]       size        size of the adj
**/
void 
InsDis::set_vet_type(int size)
{
    for(int i=0;i<size; i++)
    {
        if(a[1][i]!=OPTKIT_NULL && a[1][i]!=CAP 
            && a[0][i]==OPTKIT_NULL && a[0][i]!=CUP 
            && a[1][i] != CUP)
        {
            vet_type[i]=PI_OPEN;
        }
        else if(a[0][i]!=OPTKIT_NULL && a[0][i]!=CAP 
            && a[1][i]==OPTKIT_NULL && a[0][i]!=CUP 
            && a[1][i] != CUP)
        {
            vet_type[i]=GAMMA_OPEN;
        }
        else if(a[0][i]!=OPTKIT_NULL && a[1][i]!=OPTKIT_NULL && 
                a[0][i]!=CAP && a[1][i]!=CAP &&
                a[0][i]!=CUP && a[1][i] != CUP)
        {
            vet_type[i]=CYCLE;
        }
        else if(a[0][i]==CAP || a[1][i]==CAP)
        {
            vet_type[i]=END;
        }
        else if(a[0][i]==CUP || a[1][i]==CUP)
        {
            vet_type[i]=OUT;
        }
        else
        {
            vet_type[i]=NONE;    
        }
    }
}

/**
 * implementation of set visited
 *
 * @param[in]       vet_type        vet_type array
 * @param[in]       pos             which pos tobe set
**/
void 
InsDis::set_visited(int *vet_type, int pos)
{
    vet_type[pos] = NONE;
}

/**
 * compute distance based on number of 
 * different CCs in the graph
**/
int  
InsDis::compute_indel_dis(){ 
    int c1=0,c2=1;
    int left = 0;
    int right = 0;

    int num_c=0;
    int num_even=0;
    int num_pi_pi=0;
    int num_gamma_gamma=0;
    int num_pi_gamma=0;
    //int num_reg_even=0;
    int num_pi_odd=0;
    int num_pi_even=0;
    int num_gamma_odd=0;
    int num_gamma_even=0;
    int delta=0;
    int num_out_two=0;
    int num_out_one=0;

    /* this is dependent on which algorithm you are using */
    if(is_exem == true)
    {
        a_size = v_size; 
    }
    else
    {
        a_size = num_v_after_relabel; 
    }
    int i;
    /* open path */
    for(i=0;i<a_size;i++)
    {
        if(vet_type[i]!=NONE && vet_type[i]!=CYCLE 
            && vet_type[i]!= END && vet_type[i] != OUT)
        {
            if(vet_type[i]==PI_OPEN)
            {
                left = i;
                set_visited(vet_type, left);
                while(1)
                {
                    right = a[c2][left];
                    if(vet_type[right]==PI_OPEN)
                    {
                        num_pi_pi++;
                        set_visited(vet_type, right);
                        break;
                    }
                    if(vet_type[right]==END)
                    {
                        num_pi_odd++;
                        set_visited(vet_type, right);
                        break;
                    }
                    if(vet_type[right]==OUT)
                    {
                        num_out_one++;
                        set_visited(vet_type, right);
                        break;
                    }
                    set_visited(vet_type, right);

                    left=a[c1][right];
                    if(vet_type[left]==GAMMA_OPEN)
                    {
                        num_pi_gamma++;
                        set_visited(vet_type, left);
                        break;
                    }
                    if(vet_type[left]==END)
                    {
                        num_pi_even++;
                        set_visited(vet_type, left);
                        break;
                    }
                    if(vet_type[left]==OUT)
                    {
                        num_out_one++;
                        set_visited(vet_type, left);
                        break;
                    }
                    set_visited(vet_type, left);
                }
            }
            else if(vet_type[i]==GAMMA_OPEN)
            {
                left = i;
                set_visited(vet_type, left);
                while(1)
                {
                    right = a[c1][left];
                    if(vet_type[right]==GAMMA_OPEN)
                    {
                        num_gamma_gamma++;
                        set_visited(vet_type, right);
                        break;
                    }
                    if(vet_type[right]==END)
                    {
                        num_gamma_odd++;
                        set_visited(vet_type, right);
                        break;
                    }
                    if(vet_type[right]==OUT)
                    {
                        num_out_one++;
                        set_visited(vet_type, right);
                        break;
                    }
                    vet_type[right]=NONE;
                    set_visited(vet_type, right);
                    left=a[c2][right];
                    if(vet_type[left]==PI_OPEN)
                    {
                        num_pi_gamma++;
                        set_visited(vet_type, left);
                        break;
                    }
                    if(vet_type[left]==END)
                    {
                        num_gamma_even++;
                        set_visited(vet_type, left);
                        break;
                    }
                    if(vet_type[left]==OUT)
                    {
                        num_out_one++;
                        set_visited(vet_type, left);
                        break;
                    }
                    vet_type[left] = NONE;
                    set_visited(vet_type, left);
                }
            }
        }
    }
    /* calculate clean paths */
    int tmp_c1 = 0;
    int tmp_c2 = 0;
    for(i=0;i<a_size;i++)
    {
        if(vet_type[i]!=NONE && 
            vet_type[i]!=CYCLE && 
            vet_type[i] != OUT)
        {
            left = i;
            set_visited(vet_type, left);
            if(a[c1][left]==CAP && a[c2][left]!=CAP)
            {
                tmp_c1=1;
                tmp_c2=0;
            }
            else if(a[c1][left]!=CAP && a[c2][left]==CAP)
            {
                tmp_c1=0;
                tmp_c2=1;
            }
            else if(a[c1][left]==CAP && a[c2][left]==CAP)
            {
                num_even++;
                continue;
            }
            else if (a[c1][left]==CUP || a[c2][left]==CUP)
            {
                num_out_one++;
                continue;
            }
            while(1)
            {
                right = a[tmp_c1][left];
                if(right == OPTKIT_NULL)
                {
                    break;
                }
                if(vet_type[right]==END)
                {
                    set_visited(vet_type, right);
                    break;
                }
                if(vet_type[right]==OUT)
                {
                    num_out_one++;
                    set_visited(vet_type, right);
                    break;
                }
                set_visited(vet_type, right);

                left=a[tmp_c2][right];
                if(vet_type[left]==END)
                {
                    num_even++;
                    set_visited(vet_type, left);
                    break;
                }
                if(vet_type[left]==OUT)
                {
                    num_out_one++;
                    set_visited(vet_type, left);
                    break;
                }
                set_visited(vet_type, left);
            }
        }
    }
    /* calculate out vertices of two */
    tmp_c1 = 0;
    tmp_c2 = 0;
    for(i=0;i<a_size;i++)
    {
        if(vet_type[i] == OUT)
        {
            left = i;
            set_visited(vet_type, left);
            if(a[c1][left]==CUP && a[c2][left]!=CUP)
            {
                tmp_c1=1;
                tmp_c2=0;
            }
            else if(a[c1][left]!=CUP && a[c2][left]==CUP)
            {
                tmp_c1=0;
                tmp_c2=1;
            }
            else if(a[c1][left]==CUP && a[c2][left]==CUP)
            {
                num_out_two++;
                continue;
            }
            while(1)
            {
                right = a[tmp_c1][left];
                if(right == OPTKIT_NULL)
                {
                    break;
                }
                if(vet_type[right]==OUT || vet_type[right]==NONE)
                {
                    num_out_two++;
                    set_visited(vet_type, right);
                    break;
                }
                set_visited(vet_type, right);

                left=a[tmp_c2][right];
                if(vet_type[left]==OUT || vet_type[left]==NONE)
                {
                    num_out_two++;
                    set_visited(vet_type, left);
                    break;
                }
                set_visited(vet_type, left);
            }
        }
    }
    /* calculate cycles */
    for(i=0;i<a_size;i++)
    {
        if(vet_type[i]==CYCLE)
        {
            int start = i;
            left = i;
            set_visited(vet_type, left);
            while(1)
            {
                right = a[c1][left];
                if(right==start)
                {
                    num_c++;
                    set_visited(vet_type, right);
                    break;
                }
                vet_type[right]=NONE;

                left=a[c2][right];
                if(left==start)
                {
                    num_c++;
                    set_visited(vet_type, left);
                    break;
                }
                set_visited(vet_type, left);
            }
        }
    }

    /* calculate delta */
    if((num_pi_gamma%2==0)&&
            ((num_pi_odd>num_pi_even && num_gamma_odd>num_gamma_even)||
             (num_pi_odd<num_pi_even && num_gamma_odd<num_gamma_even)))
    {
        delta=1;
    }
    else
    {
        delta=0;
    }

    /* calculate result */
    int result = 
        ((num_c+num_pi_pi+num_gamma_gamma+(int)floor((double)num_pi_gamma/2))
         +(num_even+delta+
             num_pi_odd>num_pi_even?num_pi_even:num_pi_odd+
             num_gamma_odd>num_gamma_even?num_gamma_even:num_gamma_odd)/2);
    result += num_out_one;
    result += num_out_two;
    return result;    
}

/**
 * reverse check
 *
 * @param[in]       v       which vertex
 * @param[in]       c       which color
 *
**/
bool 
InsDis::check_if_in_to(int v, int c){
    bool result = false;
    if(vet==OPTKIT_NULL)
    {
        return result;
    }
    for(int i=0; i<adj_idx[c][vet]; i++)
    {
        if(v==adj[c][vet][i])
        {
            result=true;
            break;
        }
    }
    return result;
}

/**
 * Visualize the graph
 *
 * @param[in]       file        file tobe written
 **/
void 
InsDis::graph_vis(char *file){
    int j,from;
    FILE *writer = fopen(file, "w");
    fprintf(writer, "graph G{\n");
    for(from=0;from<v_size; from++)
    {
        int start = from==0?0:v_idx[0][from-1];
        int end = v_idx[0][from];
        for(j=start;j<end;j++)
        {
            int to = e_idx[0][j];
            fprintf(writer, "%d -- %d [color=red];\n", from, to);
        }
    }
    for(from=0;from<v_size; from++)
    {
        int start = from==0?0:v_idx[1][from-1];
        int end = v_idx[1][from];
        for(j=start;j<end;j++)
        {
            int to = e_idx[1][j];
            fprintf(writer, "%d -- %d [color=blue];\n", from, to);
        }
    }
    fprintf(writer, "}\n");
    fclose(writer);
}

/**
 * print v
 *
**/
void
InsDis::print_v()
{
    int size = is_exem==true?v_size:num_v_after_relabel;
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<size; i++)
        {
            printf("%d:%d ", i, a[c][i]);
        }
        printf("\n");
    }
    for(int i=0; i<size; i++)
    {
        cout<<map[vet_type[i]]<<" ";
    }
    printf("\n");
}

/**
 * print adj 
 *
**/
void
InsDis::print_adj()
{
    int size = is_exem==true?v_size:num_v_after_relabel;
    string map[7] = {"","pi_open","gamma_open","cycle","end","none","out"};
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<size; i++)
        {
            for(int j=0; j<adj_idx[c][i]; j++)
            {
                printf("%d:%d ", i, adj[c][i][j]);
            }
        }
        printf("\n");
    }
}

/**
 * print tmp 
 *
**/
void
InsDis::print_tmp()
{
    int size = is_exem==true?v_size:num_v_after_relabel;
    string map[7] = {"","pi_open","gamma_open","cycle","end","none","out"};
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<size; i++)
        {
            for(int j=0; j<idx_tmp[c][i]; j++)
            {
                printf("%d:%d ", i, a_tmp[c][i][j]);
            }
        }
        printf("\n");
    }
}

/**
 * get number of elements
**/
void
InsDis::get_num_elem()
{
    num_count = 3;
    num_count += (is_exem?v_size:num_v_after_relabel)*2;
    num_count += (e_size[0] + e_size[1]);
}

/**
 * print different combinations
**/
void
InsDis::print_comb()
{
    for(unsigned int i=0; i<comb.size(); i++)
    {
        for(unsigned int j=0; j<comb[i].size(); j++)
        {
            printf("%d ", comb[i][j]);
        }
        printf("\n");
    }
}

/**
 * visualize adjacency graph
**/
void
InsDis::vis_adj(){
    FILE *writer = fopen("vis_adj.dot", "w");
    fprintf(writer, "graph{\n");
    for(int c=0; c<2; c++)
    {
        const char *c_str = c==0?"red":"blue";
        for(int i=0; i<v_size; i++)
        {
            for(int j=0; j<adj_idx[c][i]; j++)
            {
                fprintf(writer, "%d -- %d [color=%s];\n", i, adj[c][i][j], c_str);
            }
        }
    }
    fprintf(writer, "}\n");
    fclose(writer);
}

/**
 * visualize tmp graph
**/
void
InsDis::vis_tmp()
{
    int size = is_exem?v_size:v_size+comb_size;
    FILE *writer = fopen("vis_tmp.dot", "w");
    fprintf(writer, "graph{\n");
    for(int c=0; c<2; c++)
    {
        const char *c_str = c==0?"red":"blue";
        for(int i=0; i<size; i++)
        {
            for(int j=0; j<idx_tmp[c][i]; j++)
            {
                fprintf(writer, "%d -- %d [color=%s];\n", i, a_tmp[c][i][j], c_str);
            }
        }
    }
    fprintf(writer, "}\n");
    fclose(writer);
}

/**
 * visualize adjacency graph
**/
void
InsDis::vis_a(int size)
{
    FILE *writer = fopen("vis_a.dot", "w");
    fprintf(writer, "graph{\n");
    for(int c=0; c<2; c++)
    {
        const char *c_str = c==0?"red":"blue";
        for(int i=0; i<size; i++)
        {
            if(a[c][i] != OPTKIT_NULL)
            {
                fprintf(writer, "%d -- %d [color=%s];\n", i, a[c][i], c_str);
            }
        }
    }
    fprintf(writer, "}\n");
    fclose(writer);
}

void
p_encode(int *encode)
{
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<encode[0]; i++)
        {
            int start = i==0?0:encode[3+c*encode[0]+i-1];
            int end = encode[3+c*encode[0]+i];
            start += ((3+2*encode[0]) + c*encode[c]);
            end += ((3+2*encode[0]) + c*encode[c]);
            for(int j=start; j<end; j++)
            {
                printf("%d->%d:%d ", j, i, encode[j]);
            }
        }
        printf("\n");
    }
}

/**
 * visualize encode
**/
void
vis_encode(int *encode)
{
    FILE *writer = fopen("vis_encode.dot", "w");
    fprintf(writer, "graph{\n");
    for(int c=0; c<2; c++)
    {
        const char *c_str = c==0?"red":"blue";
        for(int i=0; i<encode[0]; i++)
        {
            int start = i==0?0:encode[3+c*encode[0]+i-1];
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

void
InsDis::vis_csr()
{
    FILE *writer = fopen("vis_csr.dot", "w");
    fprintf(writer, "graph{\n");
    for(int c=0; c<2; c++)
    {
        const char *c_str = c==0?"red":"blue";
        for(int i=0; i<v_size; i++)
        {
            int start = i ==0?0:v_idx[c][i-1];
            int end = v_idx[c][i];
            for(int j=start; j<end; j++)
            {
                fprintf(writer, "%d -- %d [color=%s];\n", i, e_idx[c][j], c_str);
            }
        }
    }
    fprintf(writer, "}\n");
    fclose(writer);
}


int 
InsDis::return_private_info()
{
    return vet;
}

bool
InsDis::detect_no_branch(int *encode)
{
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<encode[0]; i++)
        {
            int start = i==0?0:encode[3+c*encode[0]+i-1];
            int end = encode[3+c*encode[0]+i];
            start += ((3+2*encode[0]) + c*encode[c]);
            end += ((3+2*encode[0]) + c*encode[c]);
            if((end-start) > 1)
            {
                return false;
            }
        }
    }
    return true;
}

/**
 * append the bijection mapping info to the end of the encode
 * 1) the first int is the size of the bijection mapping
 * 2) then append the real bijection information
 *
 * param[out]        encode      the encode tobe written out
 * param[in]         start       start position of the encode
 *
 * @return      total length
**/
int InsDis::append_bijection_encode(int* encode, int start)
{
    /* copy the bijection length first */
    encode[start] = bij_idx[0];
    encode[start+1] = bij_idx[1];
    start += 2;
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<bij_idx[i]; i++)
        {
            encode[start+i] = bijection[c][i];
        }
    }
    return start + 2 + bij_idx[0] + bij_idx[1];
}

/**
 * copy the bijection info into local from encode
 *
 * param[in]         encode      the encode tobe extract from
 * param[in]         start       start position of the encode
 **/
void InsDis::get_bijection_encode(int* encode, int start)
{
    bij_idx[0] = encode[start];
    bij_idx[1] = encode[start+1];
    start += 2;
    /* read into bijection */
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<bij_idx[i]; i++)
        {
            bijection[c][i] = encode[start*c+i];
        }
    }
}

bool InsDis::detect_abnornal()
{
    bool ret = false;
    for(int c=0; c<2; c++)
    {
        for(int i=0; i<v_size;i++)
        {
            for(int j=0; j<adj_idx[c][i]; j++)
            {
                int from = i;
                int to = adj[c][i][j];
                bool found = false;
                if(to == CAP)
                {
                    continue;
                }
                for(int k=0; k<adj_idx[c][to];k++)
                {
                    if(adj[c][to][k] == from)
                    {
                        found = true;
                    }
                }
                if(found == false && to != CAP)
                {
                    ret = true;
                }
            }
        }
    }
    return ret;
}

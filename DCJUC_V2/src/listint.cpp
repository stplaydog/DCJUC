#include "listint.h"
#include "instance.h"

#ifndef LIST_SIZE
#define LIST_SIZE 10000000
#endif

/**
 * The constructor.
 *
 * @param[in]       b_size
 * @param[in]       list_size
 * @param[in]       b
 * @param[in]       num_t
 * @param[in]       is_ub
 * @param[in]       num_elem
 * @param[in]       is_enumerate
 *
 * @return      N/A
 *
**/
IntList::IntList(int32_t b_size, int32_t list_size, int32_t b, int32_t num_t, bool is_ub, int num_elem, bool *is_enumerate) :
    List(b_size, list_size, b, num_t, is_ub, is_enumerate)
{ 
    content = new int*[b_size];

    for(int i=0; i<b_size; i++)
    {
        content[i] = new int[LIST_SIZE*num_elem];
    }

    c_idx = new int*[b_size];
    for(int i=0; i<b_size; i++)
    {
        c_idx[i] = new int[LIST_SIZE];
    }

    c_pos = new int[b_size];    
}

/**
 * The destructor.
 *
 * @return      N/A
**/
IntList::~IntList()
{
    for(int i=0; i<buck_size; ++i)
    {
        free(content[i]);
        free(c_idx[i]);
    }
    free(content);
    free(c_idx);
    free(c_pos);
    /* free original stuff */
    for(int i=0; i<num_threads; ++i)
    {
        free(count[i]);
        free(idx_sum[i]);
        free(end[i]);
        free(start[i]);
        free(eliminate[i]);
        free(nei_score[i]);
    }
    free(count);
    free(idx_sum);
    free(end);
    free(start);
    free(eliminate);
    free(nei_score);
    free(content_size);
    free(num);
    free(elim_idx);
}

/**
 * 
 * @return      N/A
**/
void IntList::prepare_parallel_list(int buck_id)
{
    /* get start position for each thread */
    for (int i=0; i<buck_size; ++i)
    {
        int sum = i== buck_id ? 0 : num[i];
        for(int j=0; j<num_threads; j++)
        {
            start[j][i] = sum;
            sum += count[j][i];
            end[j][i] = sum;
        }
    }

    /* compute array position */
    for(int i=0; i<buck_size; ++i)
    {
        if(i == buck_id)
        {
            continue;
        }

        int sum = c_idx[i][num[i] == 0 ? 0 : num[i] - 1];
        sum += idx_sum[0][i];
        for(int j=1; j<num_threads; j++)
        {
            int pre_pos = start[j][i] == 0 ? 0 : start[j][i] - 1;
            c_idx[i][pre_pos] = pre_pos == 0 ? 0 : sum;
            sum += idx_sum[j][i];
        }
    }
    /* cope with zero bucket */
    if(is_ub == true)
    {
        c_idx[0][0] = 0;
        int sum = idx_sum[0][buck_id];
        for(int j=1; j<num_threads; j++)
        {
            int pre_pos = start[j][buck_id] == 0 ? 0 : start[j][buck_id] - 1;
            c_idx[0][pre_pos] = pre_pos == 0 ? 0 : sum;
            sum += idx_sum[j][buck_id];
        }
    }
    else
    {
        c_idx[buck_size-1][0] = 0;
        int sum = idx_sum[0][buck_id];
        for(int j=1; j<num_threads; ++j)
        {
            int pre_pos = start[j][buck_id] == 0 ? 0 : start[j][buck_id] - 1;
            c_idx[buck_size-1][pre_pos] = pre_pos == 0 ? 0 : sum;
            sum += idx_sum[j][buck_id];
        }
    }
}

/**
 * Reset the number in the list, which means cleaning.
 *
 * @return      N/A
**/
int IntList::copy_zero(int buck_id)
{
    int ret = OPTKIT_NULL;
    if(is_ub == true)
    {
        num[0] = end[num_threads-1][buck_id];
        /* copy back */
        int copy_size = num[0] == 0 ? 0 : c_idx[0][num[0] - 1]; 
#pragma omp parallel for num_threads(num_threads) 
        for(int i=0; i<copy_size; ++i)
        {
            content[buck_id][i] = content[0][i];
        }

#pragma omp parallel for num_threads(num_threads) 
        for(int i=0; i<num[0]; ++i)
        {
            c_idx[buck_id][i] = c_idx[0][i];
            c_idx[0][i] = 0;
        }
        num[buck_id] = num[0];
        num[0]=0;
        ret = copy_size;
    }
    else
    {
        num[buck_size-1] = end[num_threads-1][buck_id];
        /* copy back */
        int copy_size = num[buck_size-1]==0?0:c_idx[buck_size-1][num[buck_size-1]-1]; 
#pragma omp parallel for num_threads(num_threads) 
        for(int i=0; i<copy_size; ++i)
        {
            content[buck_id][i] = content[buck_size-1][i];
        }

#pragma omp parallel for num_threads(num_threads) 
        for(int i=0; i<num[buck_size-1]; i++)
        {
            c_idx[buck_id][i] = c_idx[buck_size-1][i];
            c_idx[buck_size-1][i] = 0;
        }

        num[buck_id] = num[buck_size-1];
        num[buck_size-1] = 0;
        ret = copy_size;
    }

    return ret;
}

/**
 * Reset the number in the list, which means cleaning.
 *
 * @return      N/A
**/
void IntList::copy(List *other, int buck_id, int size)
{
    IntList *l = (IntList*)other;
    int start = num[buck_id]-size;    
    int offset = c_idx[buck_id][start-1];

    /* copy c_idx */
    for(int i=start, j=0; i<num[buck_id]; i++, j++)
    {
        l->c_idx[buck_id][j] = c_idx[buck_id][i]-offset;
    }

    /* copy content */
    int pos_start = c_idx[buck_id][start-1];
    int pos_end = c_idx[buck_id][num[buck_id]-1];
    for(int i=pos_start, j=0; i<pos_end; i++, j++)
    {
        l->content[buck_id][j] = content[buck_id][i];
    }

    l->num[buck_id] = size;
    num[buck_id] -= size;
}


/**
 * Reset the number in the list, which means cleaning.
 *
 * @return      N/A
**/
void IntList::print_bucket(int buck_id)
{
    printf("buck %d num %lld:", buck_id, num[buck_id]);
    for (int i=0; i<num[buck_id]; i++)
    {
        printf("|");
        int start = i==0?0:c_idx[buck_id][i-1];
        int end = c_idx[buck_id][i];
        for(int j=start; j<end; j++)
            printf("%d ", content[buck_id][j]);
    }
    printf("\n");
}


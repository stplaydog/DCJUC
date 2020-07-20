#pragma once

#include "list.h"

#ifndef _H_LISTMC
#define _H_LISTMC

using namespace std;

/**
 *  @class      IntList
 *
 *  @brief      Maximum clique list.
 *
 *  @details    This is the list for solving maximum clique problem.
**/
class IntList : public List
{
public:
    int **content;  ///< The real content that stores the encode info.
    int **c_idx;    ///< Iontent index indicate start and end position.
    int *c_pos;     ///< Indicate how many position have been filled.

    IntList(int32_t b_size, int32_t list_size, int32_t b, int32_t num_t, bool is_ub, int num_elem, bool *is_enumerate=NULL);
    ~IntList();

    /**
     * Add one element into the list.
     *
     * @param[in]       pos         which position to insert
     * @param[in]       buck_id     which bucket to insert
     * @param[in]       num_code    size of the encode
     * @param[in]       encode      the encode to be stored
     * @param[in]       ins_id      which instance
     *
     * @return      N/A
    **/
    void add(int64_t pos, int32_t buck_id, int num_code, int *encode, int ins_id)
    {
        int stt = pos == 0 ? 0 : c_idx[buck_id][pos-1];
        for(int i=0; i<num_code; i++)
        {
            content[buck_id][stt++] = encode[i];
        }

        c_idx[buck_id][pos] =stt;
    }

    /**
     * Add one element into the list.
     *
     * @param[in]       pos         which position to insert
     * @param[in]       buck_id     which bucket to insert
     * @param[in]       encode      the encode to be stored
     * @param[in]       size        size of the encode
     * @param[in]       ins_id      which instance
     *
     * @return      N/A
    **/
    void get(int64_t pos, int32_t buck_id, int *encode, int *number, int ins_id)
    {
        int start = pos == 0 ? 0 : c_idx[buck_id][pos-1];
        int end   = c_idx[buck_id][pos];
        for(int i=0, j=start; j<end; i++,j++)
        {
            encode[i] = content[buck_id][j];
        }

        *number = (end-start);
    }

    /**
     * Reset the number in the list, which means cleaning.
     *
     * @return      N/A
    **/
    void reset_num()
    {
        for(int i=lower_bound; i<upper_bound; i++)
        {
            if(start[0][i-base] < end[num_threads-1][i-base])
            {
                num[i-base] = end[num_threads-1][i-base];
            }
        }
    }

    void prepare_parallel_list(int buck_id);
    int  copy_zero(int buck_id);
    void copy(List *other, int buck_id, int size);
    void print_bucket(int buck_id);
};

#endif

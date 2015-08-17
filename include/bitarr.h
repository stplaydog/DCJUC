#pragma once
#include <stdlib.h>

#define BITS_PER_CHAR 8

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))
#define SET_BIT(var,pos) ((var) |= (1<<(pos)))
#define UNSET_BIT(var,pos) ((var) &= ~(1<<(pos)))


struct bitarr
{
        int size;
        int num;
        char *arr;
};
typedef struct bitarr tbarr;
typedef struct bitarr *pbarr;

void
create_bit_arr(pbarr arr, int size);
void
empty_bit_arr(pbarr arr);
void
add_bit_arr_pos(pbarr arr, int pos);
void
remove_bit_arr_pos(pbarr arr, int pos);


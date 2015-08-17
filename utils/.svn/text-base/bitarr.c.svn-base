#include "bitarr.h"

void 
create_bit_arr(pbarr arr, int size){
        int bit_arr_size = size/BITS_PER_CHAR +1;
        arr->arr = (char*)malloc(sizeof(char)*bit_arr_size);
        int i, j;
#pragma omp parallel for num_threads(num_threads) private(i, j)
        for(i=0;i<bit_arr_size;i++){
                for(j=0;j<BITS_PER_CHAR;j++){
                        UNSET_BIT(arr->arr[i], j);
                }
        }
        arr->size = bit_arr_size;
        arr->num = 0;
}

void    
empty_bit_arr(pbarr arr){
        int i, j;
#pragma omp parallel for num_threads(num_threads) private(i, j)
        for(i=0; i<arr->size; i++)
                for(j=0; j<BITS_PER_CHAR; j++)
                        UNSET_BIT(arr->arr[i], j);
        arr->num = 0;
}

void 
add_bit_arr_pos(pbarr arr, int pos){
        int pos1, pos2;
        pos1 = pos/BITS_PER_CHAR;
        pos2 = pos%BITS_PER_CHAR;
        //arr->arr[pos1] |= 1<<pos2; // 1 is TRUE
        SET_BIT(arr->arr[pos1], pos2);
        arr->num++;
}

void 
remove_bit_arr_pos(pbarr arr, int pos){
        int pos1, pos2;
        pos1 = pos/BITS_PER_CHAR;
        pos2 = pos%BITS_PER_CHAR;
        //arr->arr[pos1] |= 0<<pos2; // 0 is FALSE
        UNSET_BIT(arr->arr[pos1], pos2);
        arr->num--;
}

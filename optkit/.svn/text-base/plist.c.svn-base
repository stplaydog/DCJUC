#include "plist.h"

void 
create_plist(pbnbl l,
    void (*fpc)(void*, int), void (*fpp)(void*, int), 
    int buck_size, int base){
    int i, j, k;
    l->buck_size = buck_size;
    l->base = base;
    l->child_size = (int*)malloc(sizeof(int)*buck_size);
    l->parent_size = (int*)malloc(sizeof(int)*buck_size);
    l->child_num = (int*)malloc(sizeof(int)*buck_size);
    l->parent_num = (int*)malloc(sizeof(int)*buck_size);
    l->child = (void**)malloc(sizeof(void*)*buck_size);   	 
    l->parent = (void**)malloc(sizeof(void*)*buck_size);   	 
    l->avail = (int**)malloc(sizeof(int*)*buck_size);
    for(i=0; i<buck_size; i++){
        l->avail[i] = (int*)malloc(sizeof(int)*LIST_SIZE);
	l->child_size[i] = LIST_SIZE;
	l->parent_size[i] = LIST_SIZE;
	l->child_num[i] = 0;
	l->parent_num[i] = 0;
	(*fpc)(l->child[i], LIST_SIZE);
        (*fpp)(l->parent[i], LIST_SIZE);
    }
}

void 
free_plist(pbnbl l, 
    void (*fpc)(void*, int), void (*fpp)(void*, int)){
    int i;
    for(i=0; i<l->buck_size; i++){
        (*fpc)(l->child[i], l->child_size[i]);
        (*fpp)(l->parent[i], l->parent_size[i]);
    }
    free(l->child_size);
    free(l->parent_size);
    free(l->child_num);
    free(l->parent_num);
    free(l->parent);
    free(l->child);
    free(l);
}

void
add_parent(pbnbl l, void *parent, int pos_p, 
    void (*a_p)(void*, void*, int)){
    (*a_p)(l, parent, pos_p);
}

void
add_child(pbnbl l, void *child, int pos_c, int pos_p, 
    void (*a_c)(void*, void*, int, int)){
    (*a_c)(l, child, pos_c, pos_p);
}

void
get_list(pbnbl l, void *elem, int buck_id, int pos_c,
   void (*g_l)(void*, void*, int, int)){
    (*g_l)(l, elem, buck_id, pos_c);
}

void
double_child(pbnbl l, int buck_id,
	void (*d_c)(void*, int)){
	(*d_c)(l, buck_id);
}

void
double_parent(pbnbl l, int buck_id,
	void (*d_p)(void*, int)){
	(*d_p)(l, buck_id);
}

#pragma once

#include <stdio.h>
#include <stdlib.h>

#define LIST_SIZE 1000

void (*c_c)(void*, int);
void (*c_p)(void*, int);
void (*f_c)(void*, int);
void (*f_p)(void*, int);
void (*a_c)(void*, void*, int);
void (*a_p)(void*, void*, int);
void (*g_l)(void*, void*, int);

struct bnb_list
{
	int buck_size;
	int base;
	int *child_size;
	int *parent_size;
	int *child_num;
	int *parent_num;
	int **avail;
	void **child;	
	void **parent;	
};
typedef struct bnb_list tbnbl;
typedef struct bnb_list *pbnbl;

void 
create_plist(pbnbl l,
    void (*c_c)(void*, int), void (*c_p)(void*, int), 
    int buck_size, int base);
void 
free_plist(pbnbl l, 
    void (*f_c)(void*, int), void (*f_p)(void*, int));
void
add_parent(pbnbl l, void *parent, int pos_p, 
    void (*a_p)(void*, void*, int));
void
add_child(pbnbl l, void *child, int pos_c, int pos_p, 
    void (*a_c)(void*, void*, int, int));
void
get_list(pbnbl l, void *elem, int buck_id, int pos_c,
    void (*g_l)(void*, void*, int, int));

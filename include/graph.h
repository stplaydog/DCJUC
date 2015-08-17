#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <immintrin.h>

#define MAX_GENE_SIZE 100000
#define TRUE 0
#define FALSE 1
#define NUL -1
#define VISITED -1
#define NUM_AS 16
#define CAP 65533
#define TYPE_MUL 0
#define TYPE_SING 1

//struct time_stamp{
//	int time;
//	int stamp;
//};
//typedef struct time_stamp ts;
//typedef struct time_stamp *pts;

#define ERROR_PRINT() {printf("Error on (%d) line in (%s) file\n", __LINE__, __FILE__); exit(123);}
struct gr
{
	int v_size;
	int **v_idx; //4*sth
	int **v_deg;
	int *v_check;

	int **e_idx;
	int *e_size;
	int **e_check;
	int *zero;
	int cycle;
	int distance[3];
};
typedef struct gr tg;
typedef struct gr *pg;


struct gr_adj
{
	int ***adj;
	int **idx;
	int *e_size;
	int v_size;
	int c_size;
};
typedef struct gr_adj tga;
typedef struct gr_adj *pga;

void init_gr(pg g, char *file);
void rename_graph_by_car(pg g);
void graph_vis(pg g, char *file);
void graph_vis_two(pg g, char *file, int c1, int c2);
void prepare_dis(pg g);
void prepare_graph_for_distance(pg g1, pg g2);
void free_graph(pg g);
void 
free_graph_two(pg g);
void 
indentify_comps(int** comp_id, int **comp_type, pg g);
void from_gr_to_adj(pg g, pga a);
